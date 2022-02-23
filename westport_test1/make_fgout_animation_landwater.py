"""
Make an mp4 animation of fgout grid results. 
This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.
The tuple update_artists contains the list of Artists that must be changed
in update.  Modify this as needed.

This version also introduces some debris particles to be tracked.

Note that imshow is used here, and to apply imshow to an array from fgout,
you must transpose and flipud, and then specify the extent as 
fgout.extent_edges, e.g. to plot the topography:
   imshow(flipud(fgout.B.T), extent=fgout1.extent_edges, cmap=...)

Alternatively, could use plottools.pcolorcells, e.g.
   plottools.pcolorcells(fgout.X,fgout.Y,fgout.B, cmap=...)
but that doesn't work well when a transparency is specified for viewing plots
on top of a background image.  The edges drawn by pcolormesh come out darker
than the faces, even when trying to adjust facecolor and alpha.  Seems to be
a known bug with the underlying matplotlib routine pcolormesh.

"""

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from clawpack.geoclaw import fgout_tools
from matplotlib import animation, colors

sys.path.insert(0,'/Users/rjl/git/geoclaw_debris/common_python')
import debris_tools
    
fgno = 11  # which fgout grid

outdir = '_output'
format = 'binary'  # format of fgout grid output

# list of frames of fgout solution to use in animation:
fgframes = range(200,230)

figsize = (10,8)

bgimage = None  # if None, then color plots of B0 will be used.

if 1:
    # provide a background image 
    graphics_dir = '/Users/rjl/git/WestportMaritime/graphics/'
    bgimage = imread(graphics_dir+'fgout11CT.png')
    bgimage_extent = [-124.16, -124.08, 46.885, 46.92]  # corners of bgimage

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 

# initial topography (used to determine land / water points):
fgout0 = fgout_grid.read_frame(1)  # assumed to be at time 0 (no subsidence)
B0_fcn = fgout_tools.make_fgout_fcn_xy(fgout0, 'B')

# -----------------
# Debris particles

# Deterime time t of first fgout frame, to initialize particles
fgout1 = fgout_grid.read_frame(fgframes[0])
t0 = fgout1.t
fgout_extent = fgout1.extent_edges

# make a B(x,y) function that can be used to determine land vs. water:
#B_fcn_xy = fgout_tools.make_fgout_fcn_xy(fgout1, 'B')

# Initialize debris_paths dictionary and set
# initial debris particle locations (x0,y0) at time t0.
# Require a list of dbnos and each array 
#     debris_paths[dbno]
# in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

debris_paths = {}
dbnos_water = []
dbnos_land = []
grounding_depth = {}
drag_factor = {}

# set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
u0 = 0.
v0 = 0.
x1 = fgout_extent[0] + 0.001
x2 = fgout_extent[1] - 0.0005
y1 = fgout_extent[2] + 0.001
y2 = fgout_extent[3] - 0.0005

dxd = 0.002
dyd = dxd * cos(47*pi/180)
xd = arange(x1, x2, dxd)
yd = arange(y1, y2, dyd)
mxd = len(xd)
myd = len(yd)

for i in range(len(xd)):
    x0 = xd[i]
    for j in range(len(yd)):
        y0 = yd[j]
        dbno = myd*i + j
        db = array([[t0, x0, y0, u0, v0]])
        debris_paths[dbno] = db
        
        if B0_fcn(x0,y0) < 0:
            dbnos_water.append(dbno)
            grounding_depth[dbno] = 0.
            drag_factor[dbno] = None
        elif B0_fcn(x0,y0) > 0:
            dbnos_land.append(dbno)
            grounding_depth[dbno] = 1.
            drag_factor[dbno] = 0.02

dbnos = dbnos_water + dbnos_land
print('Created %i initial debris particles' % len(dbnos))

# Compute debris path for each particle by using all the fgout frames
# in the list fgframes (first frame should be one used to set t0 above):

debris_paths = debris_tools.make_debris_paths(fgframes, fgout_grid, 
                            debris_paths, dbnos, drag_factor, grounding_depth)


# done computing debris paths
# ------------------------------

# ----------
# Plotting:

# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

fgout1 = fgout_grid.read_frame(fgframes[0])

plot_extent = fgout1.extent_edges

ylat = fgout1.Y.mean()  # for aspect ratio of plots

fig,ax = subplots(figsize=figsize)

if bgimage is not None:
    ax.imshow(bgimage,extent=bgimage_extent)
else:
    # plot initial B0 topo  (created separately for this test data):
    imshow(flipud(fgout0.B.T), extent=fgout0.extent_edges,
                          cmap=geoplot.land_colors)
    clim(-5,5)
    
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])

# plot some quantity of interest (here the depth h):

# set transparency alpha (1=opaque)
if bgimage is None:
    a = 1
else:
    a = 0.4
    
bounds_depth = array([1e-6,1,2,3])
cmap_depth = colors.ListedColormap([[.3,.3,1,a],[0,0,1,a],[.5,0,.5,a]])

# Set color to transparent where no water:
cmap_depth.set_under(color=[1,1,1,0])
norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)

qoi_plot = imshow(flipud(fgout1.h.T), extent=fgout1.extent_edges,
                  cmap=cmap_depth, norm=norm_depth)
                  
cb = colorbar(qoi_plot, extend='max', shrink=0.7)
cb.set_label('meters')
title_text = title('Depth h at %s after quake' % fgout1.t_hms)

# add debris:
            
xd_water,yd_water = debris_tools.get_debris_xy(fgout1.t, debris_paths, dbnos_water)
xd_land,yd_land = debris_tools.get_debris_xy(fgout1.t, debris_paths, dbnos_land)

points_water, = ax.plot(xd_water, yd_water, color='b',linestyle='',
                        marker='.',markersize=4)
points_land, = ax.plot(xd_land, yd_land, color='r',linestyle='',
                        marker='.',markersize=4)


ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])


# The artists that will be updated for subsequent frames:
update_artists = (qoi_plot, title_text, points_water, points_land)
        
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    print('Updating plot at time %s' % fgout.t_hms)    
    
    # unpack update_artists (must agree with definition above):
    qoi_plot, title_text, points_water, points_land = update_artists
        
    title_text.set_text('Depth h at %s after quake' % fgout.t_hms)
    qoi_plot.set_data(flipud(fgout.h.T))  # for imshow
        
    # update debris points:
    xd_water,yd_water = debris_tools.get_debris_xy(fgout.t, debris_paths,
                                                   dbnos_water)
    xd_land,yd_land = debris_tools.get_debris_xy(fgout.t, debris_paths,
                                                   dbnos_land)
    points_water.set_data(xd_water,yd_water)
    points_land.set_data(xd_land,yd_land)
    
    update_artists = (qoi_plot, title_text, points_water, points_land)
    return update_artists

def plot_fgframe(fgframeno):
    """
    Convenience function for plotting one frame interactively.
    But if you use this function in IPython and then try to make the animation,
    it may get into an infinite loop (not sure why).  Close the figure to abort.
    """
    update(fgframeno, *update_artists)
                

def make_anim():
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=fgframes, 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':

    anim = make_anim()
    
    # Output files:
    name = 'fgout11_landwater_anim'

    fname_mp4 = name + '.mp4'
    fname_html = None # name + '.html'
    
    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        fname_html = name + '.html'
        animation_tools.make_html(anim, file_name=fname_html, title=name)
    
    
    
    
