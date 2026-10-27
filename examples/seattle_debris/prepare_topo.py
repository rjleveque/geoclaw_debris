from pathlib import Path
import hashlib
import json

import netCDF4
import numpy as np
import pandas as pd
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely import wkt
from shapely.geometry import box
from clawpack.geoclaw import topotools

from setrun import X1, X2, Y1, Y2

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / 'SeattleBuildings'
MHW_NAME = 'puget_sound_mhw_19_n47x75_w122x50_2023v1.nc'


def read_mhw_region(source):
    # read the model domain and one extra cell to avoid edge effects in geoclaw
    with netCDF4.Dataset(source) as ds:
        lon = np.asarray(ds['lon'][:], dtype=float)
        lat = np.asarray(ds['lat'][:], dtype=float)
        if not (np.all(np.diff(lon) > 0) and np.all(np.diff(lat) > 0)): raise ValueError('MHW longitude and latitude must increase')
        col0 = np.searchsorted(lon, X1, side='left') - 1
        col1 = np.searchsorted(lon, X2, side='right') + 1
        row0 = np.searchsorted(lat, Y1, side='left') - 1
        row1 = np.searchsorted(lat, Y2, side='right') + 1
        if min(col0, row0) < 0 or col1 > len(lon) or row1 > len(lat): raise ValueError('MHW terrain does not cover the simulation domain')
        x = lon[col0:col1]
        y = lat[row0:row1]
        raw = ds['z'][row0:row1, col0:col1]
        if np.ma.getmaskarray(raw).any() or not np.isfinite(raw).all(): raise ValueError('some missing elevations in the simulation domain')
        return x, y, np.asarray(raw, dtype=float)

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
    source = str(HERE / MHW_NAME)
    inventory = SOURCE / 'Seattle_Inventory.csv'
    x, y, bare = read_mhw_region(source)
    transform = from_origin(x[0] - (x[1]-x[0])/2, y[-1] + (y[1]-y[0])/2, x[1]-x[0], y[1]-y[0])
    domain = box(x[0], y[0], x[-1], y[-1])
    buildings = [wkt.loads(value) for value in pd.read_csv(inventory).geometry]
    buildings = [p for p in buildings if p.intersects(domain)]
    if any(not p.is_valid for p in buildings): raise ValueError('Invalid building geometry')
    footprint = rasterize([(p, 1) for p in buildings], out_shape=bare.shape, transform=transform, fill=0)[::-1]
    # force pier buildings above water even where the seabed is deep
    combined = np.where(footprint, np.maximum(bare, 0.) + 20., bare)
    write_topo(HERE / 'seattle_19arcsec_mhw.tt3', x, y, bare)
    terrain = HERE / 'seattle_buildings_mhw.tt3'
    write_topo(terrain, x, y, combined)
    digest = file_hash(terrain)
    # used by combine_fgout_nc.py
    (HERE / 'terrain_provenance.json').write_text(json.dumps({'source': source, 'vertical_datum': 'MHW', 'combined_sha256': digest}, indent=2) + '\n')
    print(f'wrote MHW terrain from {source} with {len(buildings)} buildings; SHA256: {digest}')

if __name__ == '__main__':
    main()
