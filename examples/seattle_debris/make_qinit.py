import numpy as np
from clawpack.geoclaw import topotools
from setrun import X1, X2, Y1, Y2, SURGE_HEIGHT_M
#keep setrun and qinit xs and ys synced
x1, x2 = X1, X2
y1, y2 = Y1, Y2
cellsize = 0.0001
x = np.arange(x1, x2 + cellsize/2, cellsize)
y = np.arange(y1, y2 + cellsize/2, cellsize)
X, Y = np.meshgrid(x, y) # make the grid we're gonna be working with

width = X2 - X1
x_ramp_start = X1 + 0.15 * width
x_ramp_end = X1 + 0.55 * width
eta_pert = np.where(X < x_ramp_start, SURGE_HEIGHT_M, np.where(X > x_ramp_end, 0.0, SURGE_HEIGHT_M * (x_ramp_end-X)/(x_ramp_end-x_ramp_start)))

topo = topotools.Topography()
topo.set_xyZ(x, y, eta_pert)
topo.write('qinit.tt1', topo_type=1)
print(f'wrote qinit.tt1: {len(x)}x{len(y)}, 'f'range [{eta_pert.min():.1f}, {eta_pert.max():.1f}] m')
