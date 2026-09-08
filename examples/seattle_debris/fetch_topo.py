##download high-res topography/bathymetry for elliot bay from the NOAA puget sound 1/3 arc-second DEM using Clawpack's topotools


from clawpack.geoclaw import topotools

# elliot bay and seattle waterfront
extent = [-122.40, -122.32, 47.58, 47.64]
topo = topotools.read_netcdf('puget_sound', extent=extent, verbose=True)
fname = 'elliott_bay.tt3'
topo.write(fname, topo_type=3)
print(f'Saved {fname}: shape={topo.Z.shape}, '
      f'Z range [{topo.Z.min():.1f}, {topo.Z.max():.1f}] m')
