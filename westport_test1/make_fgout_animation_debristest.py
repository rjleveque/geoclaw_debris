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
from matplotlib import animation, colors
from datetime import timedelta 


if 1:
    from clawpack.geoclaw import fgout_tools
    sys.path.insert(0,'../common_python')
    import debris_tools
    graphics_dir = '/Users/rjl/git/WestportMaritime/graphics/'
else:
    # local versions for self-contained directory:
    import fgout_tools, debris_tools
    graphics_dir = './'
    
fgno = 11  # which fgout grid

outdir = '_output_S-A-Whole'
format = 'binary'  # format of fgout grid output

# list of frames of fgout solution to use in animation:
fgframes = range(200,230)

figsize = (10,8)

# background image?

#bgimage = None  # if None, then color plots of fgout.B will be used.

bgimage = imread(graphics_dir+'fgout11CT.png')
bgimage_extent = [-124.16, -124.08, 46.885, 46.92]  # corners of bgimage

# For plotting debris particles:
color = 'r'
marker = 's'
markersize = 4

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 

# -----------------
# Debris particles

# Deterime time t of first fgout frame, to initialize particles
fgout1 = fgout_grid.read_frame(fgframes[0])
t0 = fgout1.t

# Initialize debris_paths dictionary and set
# initial debris particle locations (x0,y0) at time t0.
# Require a list of dbnos and each array 
#     debris_paths[dbno]
# in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

debris_paths = {}
dbnos = []
grounding_depth = {}
drag_factor = {}

# set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
u0 = 0.
v0 = 0.

# Choose initial locations, in this case points along a line at constant longitude
xd = -124.125
yd = arange(46.89, 46.915, .002)

for j in range(len(yd)):
    x0 = xd
    y0 = yd[j]
    
    # a massless tracer particle:
    dbno = j
    db = array([[t0, x0, y0, u0, v0]])
    debris_paths[dbno] = db
    grounding_depth[dbno] = 0.
    drag_factor[dbno] = None
    dbnos.append(dbno)

    # make another particle with a drag factor, for comparison
    # (same starting location)
    dbno = j + 100
    db = array([[t0, x0, y0, u0, v0]])
    debris_paths[dbno] = db
    grounding_depth[dbno] = 0.
    drag_factor[dbno] = 0.02
    dbnos.append(dbno)

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
    imshow(flipud(fgout1.B.T), extent=fgout1.extent_edges,
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
title_text = title('Depth h at %s after quake' % timedelta(seconds=fgout1.t))

# add debris:

xd,yd = debris_tools.get_debris_xy(fgout1.t, debris_paths, dbnos)
# plot points.  Note that plot returns a list with a single element, 
# extract this as points via unpacking as points,
points, = ax.plot(xd, yd, color=color, linestyle='', marker=marker, 
                     markersize=markersize)                    

ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])


# The artists that will be updated for subsequent frames:
update_artists = (qoi_plot, title_text, points)
        
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    print('Updating plot at time %s' % timedelta(seconds=fgout.t))    
    
    # unpack update_artists (must agree with definition above):
    qoi_plot, title_text, points = update_artists
        
    title_text.set_text('Depth h at %s after quake' % timedelta(seconds=fgout.t))
    qoi_plot.set_data(flipud(fgout.h.T))  # for imshow
        
    # update debris points:
    xd,yd = debris_tools.get_debris_xy(fgout.t, debris_paths, dbnos)
    points.set_data(xd,yd)
    
    update_artists = (qoi_plot, title_text, points)
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
    name = 'fgout11_debris_anim'

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
    
    
    
    
