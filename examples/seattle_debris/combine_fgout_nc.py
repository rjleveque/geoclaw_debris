import json
from pathlib import Path
import time

import netCDF4
from clawpack.geoclaw import fgout_tools
from prepare_topo import file_hash


def combine():
    outdir = Path('_output').resolve()
    fgno = 1
    provenance = json.loads((outdir / 'terrain_provenance.json').read_text())
    terrain = Path(__file__).resolve().parent / 'seattle_buildings_mhw.tt3'
    if file_hash(terrain) != provenance['combined_sha256']: raise ValueError('provenance does not match the current terrain')
    grid = fgout_tools.FGoutGrid(fgno, str(outdir), 'binary32')
    grid.read_fgout_grids_data()
    for n in range(1, grid.nout+1):
        if not (outdir / f'fgout{fgno:04d}.b{n:04d}').exists(): raise ValueError(f'incomplete: missing frame {n}/{grid.nout}')
    first = grid.read_frame(1)
    destination = Path('pike_place_buffered_fgout.nc')
    temporary = destination.with_suffix('.partial.nc')
    with netCDF4.Dataset(temporary, 'w') as ds:
        ds.description = 'uniform 3 m Pike Place waterfront flow on corrected terrain'
        ds.history = f'created {time.ctime()} from {outdir}'
        ds.terrain_sha256 = provenance['combined_sha256']
        ds.fgno = fgno
        for name, values, units in [('lon', first.x, 'degrees_east'), ('lat', first.y, 'degrees_north')]:
            ds.createDimension(name, len(values))
            var = ds.createVariable(name, 'f8', (name,))
            var[:] = values
            var.units = units
        ds.createDimension('time', grid.nout)
        ds.createVariable('time', 'f8', ('time',)).units = 'seconds'
        for name in ['h', 'u', 'v', 'B']:
            var = ds.createVariable(name, 'f4', ('time', 'lon', 'lat'), zlib=True, complevel=1)
            var.units = 'm/s' if name in ['u', 'v'] else 'meters'
        for n in range(1, grid.nout+1):
            frame = first if n == 1 else grid.read_frame(n)
            ds['time'][n-1] = frame.t
            for name in ['h', 'u', 'v', 'B']: ds[name][n-1] = getattr(frame, name)
    temporary.replace(destination)
    print(f'created {destination} from fgout grid {fgno}')

if __name__ == '__main__':
    combine()
