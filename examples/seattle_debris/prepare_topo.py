from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.windows import Window, from_bounds
from shapely import wkt
from shapely.geometry import box
from clawpack.geoclaw import topotools

from setrun import X1, X2, Y1, Y2

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / 'SeattleBuildings'

def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_topo(path, x, y, south_first):
    topo = topotools.Topography()
    topo.set_xyZ(x, y, south_first)
    # .tt3 file rows are written north to south
    topo.write(str(path), topo_type=3, Z_format='%.3f', grid_registration='lower')
    restored = topotools.Topography(str(path), topo_type=3)
    np.testing.assert_allclose(restored.x, x, rtol=0, atol=1e-10)
    np.testing.assert_allclose(restored.y, y, rtol=0, atol=1e-10)
    np.testing.assert_allclose(restored.Z, south_first, rtol=0, atol=0.00051)

def main():
    tif = SOURCE / 'seattle_19arcsec_navd88_wgs84.tif'
    inventory = SOURCE / 'Seattle_Inventory.csv'
    with rasterio.open(tif) as src:
        if src.crs.to_epsg() != 4326 or src.transform.a <= 0 or src.transform.e >= 0 or src.transform.b or src.transform.d: raise ValueError('wrong terrain format')
        # add one cell around the model domain.
        region = from_bounds(X1, Y1, X2, Y2, src.transform)
        col0, row0 = int(np.floor(region.col_off))-1, int(np.floor(region.row_off))-1
        col1 = int(np.ceil(region.col_off+region.width))+1
        row1 = int(np.ceil(region.row_off+region.height))+1
        if col0 < 0 or row0 < 0 or col1 > src.width or row1 > src.height: raise ValueError('terrain does not cover the simulation domain')
        window = Window(col0, row0, col1-col0, row1-row0)
        transform = src.window_transform(window)
        raw = src.read(1, window=window, masked=True)
        if np.ma.getmaskarray(raw).any() or not np.isfinite(raw).all(): raise ValueError('some missing elevations in the simulation domain')
        bare = np.asarray(raw, dtype=float) - 2.484

    # use cell centers, then put latitude in south-to-north order
    x = transform.c + (np.arange(bare.shape[1]) + 0.5) * transform.a
    y = (transform.f + (np.arange(bare.shape[0]) + 0.5) * transform.e)[::-1]
    domain = box(x[0], y[0], x[-1], y[-1])
    buildings = [wkt.loads(value) for value in pd.read_csv(inventory).geometry]
    buildings = [p for p in buildings if p.intersects(domain)]
    if any(not p.is_valid for p in buildings): raise ValueError('Invalid building geometry')
    footprint = rasterize([(p, 1) for p in buildings], out_shape=bare.shape, transform=transform, fill=0)
    # flip the raster and footprint arrays to match the south-to-north latitude order
    bare = bare[::-1].copy()
    footprint = footprint[::-1]
    # force pier buildings above water even where the seabed is deep
    combined = np.where(footprint, np.maximum(bare, 0.) + 20., bare)
    write_topo(HERE / 'seattle_19arcsec_mhw.tt3', x, y, bare)
    terrain = HERE / 'seattle_buildings_mhw.tt3'
    write_topo(terrain, x, y, combined)
    digest = file_hash(terrain)
    (HERE / 'terrain_provenance.json').write_text(json.dumps({'combined_sha256': digest})+'\n')
    print(f'Wrote terrain with {len(buildings)} buildings; SHA256: {digest}')

if __name__ == '__main__':
    main()
