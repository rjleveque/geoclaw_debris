# 3m surface perturbation on the west side, tapering to 0 eastward (so that the water is moving, not just a still lake)

import numpy as np
from clawpack.geoclaw import topotools
from setrun import X1, X2, Y1, Y2
#keep setrun and qinit xs and ys synced
x1, x2 = X1, X2
y1, y2 = Y1, Y2
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
