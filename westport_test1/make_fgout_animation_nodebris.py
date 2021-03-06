"""
Make an mp4 animation of fgout grid results. 
This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.
The tuple update_artists contains the list of Artists that must be changed
in update.  Modify this as needed.

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
qoi = 'speed'  # what to plot, 'h', 'speed' or 'eta'

outdir = '_output_S-A-Whole'
format = 'binary'  # format of fgout grid output

fgframes = range(200,230)  # frames of fgout solution to use in animation

figsize = (10,8)

# background image?

bgimage = None  # if None, then color plots of fgout.B will be used.

bgimage = imread(graphics_dir+'fgout11CT.png')
bgimage_extent = [-124.16, -124.08, 46.885, 46.92]  # corners of bgimage

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 


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


# set trasparency alpha for imshow plots of qoi (1=opaque)
if bgimage is None:
    a = 1
else:
    a = 0.4

if qoi == 'eta':
    eta = ma.masked_where(fgout1.h<0.001, fgout1.eta)

    qoi_plot = imshow(flipud(eta.T), extent=fgout1.extent_edges,
                      cmap=geoplot.tsunami_colormap)
    
    clim(-5,5)
    cb = colorbar(qoi_plot, extend='both', shrink=0.7)
    cb.set_label('meters')
    title_text = title('Surface eta at %s after quake' \
                % timedelta(seconds=fgout1.t))

    
elif qoi == 'h':
        
    bounds_depth = array([1e-6,1,2,3])
    cmap_depth = colors.ListedColormap([[.3,.3,1,a],[0,0,1,a],[.5,0,.5,a]])
    
    # Set color to transparent where no water:
    cmap_depth.set_under(color=[1,1,1,0])
    norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)

    qoi_plot = imshow(flipud(fgout1.h.T), extent=fgout1.extent_edges,
                      cmap=cmap_depth, norm=norm_depth)
                      
    cb = colorbar(qoi_plot, extend='max', shrink=0.7)
    cb.set_label('meters')
    title_text = title('Depth h at %s after quake' \
                % timedelta(seconds=fgout1.t))

elif qoi == 'speed':
        
    bounds_speed = array([1e-6,2,4,6])
    cmap_speed = colors.ListedColormap([[.3,.3,1,a],[0,0,1,a],[.5,0,.5,a]])
    
    # Set color to transparent where no water:
    cmap_speed.set_under(color=[1,1,1,0])
    norm_speed = colors.BoundaryNorm(bounds_speed, cmap_speed.N)

    qoi_plot = imshow(flipud(fgout1.s.T), extent=fgout1.extent_edges,
                      cmap=cmap_speed, norm=norm_speed)
                      
    cb = colorbar(qoi_plot, extend='max', shrink=0.7)
    cb.set_label('meters/sec')
    title_text = title('Speed s at %s after quake' % timedelta(seconds=fgout1.t))


ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])


# The artists that will be updated for subsequent frames:
update_artists = (qoi_plot, title_text)
        
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    print('Updating plot at time %s' % timedelta(seconds=fgout.t))    
    
    # unpack update_artists (must agree with definition above):
    qoi_plot, title_text = update_artists
        
    if qoi == 'eta':
        title_text.set_text('Surface eta at %s after quake' % timedelta(seconds=fgout.t))
        eta = ma.masked_where(fgout.h<0.001, fgout.eta)
        qoi_plot.set_data(flipud(eta.T))  # for imshow
        
    elif qoi == 'h':
        title_text.set_text('Depth h at %s after quake' % timedelta(seconds=fgout.t))
        qoi_plot.set_data(flipud(fgout.h.T))  # for imshow

    elif qoi == 'speed':
        title_text.set_text('Speed s at %s after quake' % timedelta(seconds=fgout.t))
        qoi_plot.set_data(flipud(fgout.s.T))  # for imshow
        
    update_artists = (qoi_plot, title_text)
    return update_artists

def plot_fgframe(fgframeno):
    """
    Convenience function for plotting one frame.
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
    name = 'fgout11_anim'

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
    
    
    
    
