# 3m surface perturbation on the west side, tapering to 0 eastward (so that the water is moving, not just a still lake)

import numpy as np
from clawpack.geoclaw import topotools

x1, x2 = -122.40, -122.32
y1, y2 = 47.58, 47.64
cellsize = 0.0004
x = np.arange(x1, x2 + cellsize/2, cellsize)
y = np.arange(y1, y2 + cellsize/2, cellsize)
X, Y = np.meshgrid(x, y) # make the grid we're gonna be working with

x_ramp_start = -122.38
x_ramp_end = -122.36
eta_pert = np.where(X<x_ramp_start, 3.0, np.where(X>x_ramp_end, 0.0, 3.0 * (x_ramp_end-X)/(x_ramp_end-x_ramp_start))) # makes a linear ramp of depth between the start and end points, from 3m deep to 0m deep

topo = topotools.Topography()
topo.set_xyZ(x, y, eta_pert)
topo.write('qinit.tt1', topo_type=1)
print(f'wrote qinit.tt1: {len(x)}x{len(y)}, 'f'range [{eta_pert.min():.1f}, {eta_pert.max():.1f}] m')
