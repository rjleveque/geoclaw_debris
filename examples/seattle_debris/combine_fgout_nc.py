# combine fgout binary frames into a single netcdf file

from clawpack.geoclaw import fgout_tools

fgout_grid = fgout_tools.FGoutGrid(1, '_output', 'binary32')
fgout_grid.read_fgout_grids_data()

fgout_frames = []
for n in range(1, fgout_grid.nout + 1):
    fgout_frames.append(fgout_grid.read_frame(n))

fgout_tools.write_netcdf(fgout_frames, fname_nc='elliott_bay_fgout.nc',
                         qois=['h', 'u', 'v', 'B'], datatype='f4',
                         include_B0=False, include_Bfinal=False,
                         verbose=True)
