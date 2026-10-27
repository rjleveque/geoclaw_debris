import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from clawpack.geoclaw import topotools


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('flow', type=Path)
    parser.add_argument('--bare-topo', type=Path, default=Path(__file__).with_name('seattle_19arcsec_mhw.tt3'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('street_flooding.png'))
    parser.add_argument('--threshold', type=float, default=0.2)
    parser.add_argument('--view', type=float, nargs=4, default=(-122.356, -122.339, 47.606, 47.616))
    args = parser.parse_args()
    if args.threshold <= 0: parser.error('--threshold must be positive')

    with Dataset(args.flow) as ds:
        lon, lat, times = [np.asarray(ds[name][:]) for name in ('lon', 'lat', 'time')]
        depth = np.ma.filled(ds['h'][:], np.nan)
        terrain = np.ma.filled(ds['B'][0], np.nan)
    topo = topotools.Topography(str(args.bare_topo), topo_type=3)
    x, y = np.meshgrid(lon, lat, indexing='ij')
    bare = RegularGridInterpolator((topo.x, topo.y), topo.Z.T)(np.column_stack((x.ravel(), y.ravel()))).reshape(x.shape)
    if not np.isfinite(depth).all() or not np.isfinite(terrain).all(): raise ValueError('flow contains missing elevations or depths')

    building = terrain - bare >= 10
    open_land = (bare >= 0) & ~building & (depth[0] < args.threshold)
    wet = (depth >= args.threshold).any(axis=0) & open_land
    arrival = np.where(wet, times[(depth >= args.threshold).argmax(axis=0)], np.nan)
    maximum = np.where(wet, depth.max(axis=0), np.nan)

    west, east, south, north = args.view
    ix = np.flatnonzero((lon >= west) & (lon <= east))
    iy = np.flatnonzero((lat >= south) & (lat <= north))
    if len(ix) < 2 or len(iy) < 2: parser.error('--view does not overlap the flow grid')
    view = np.ix_(ix, iy)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    panels = ((maximum, 'max water depth', 'Blues', 'm'), (arrival, f'first arrival at {args.threshold:g} m depth', 'magma_r', 'seconds'))
    for ax, (values, title, cmap, units) in zip(axes, panels):
        ax.set_facecolor('#eee9df')
        image = ax.pcolormesh(lon[ix], lat[iy], np.ma.masked_invalid(values[view].T), cmap=cmap, shading='nearest', vmin=0)
        b = building[view].T
        ax.pcolormesh(lon[ix], lat[iy], np.ma.masked_where(~b, b), cmap=ListedColormap(['#464d56']), shading='nearest')
        fig.colorbar(image, ax=ax, label=units, shrink=0.72)
        ax.set(xlim=(west, east), ylim=(south, north), xlabel='Longitude', ylabel='Latitude', title=title)
        ax.set_aspect(1 / np.cos(np.deg2rad(lat[iy].mean())))
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    axes[0].legend(handles=[Patch(facecolor='#464d56', label='building cells'), Patch(facecolor='#eee9df', label='unflooded or excluded')], loc='upper left', fontsize=8)
    flooded = np.isfinite(maximum[view]).sum()
    dry = open_land[view].sum()
    fig.suptitle(f'seattle waterfront: {flooded} of {dry} initially dry open-ground cells flooded in this view')
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(f'{args.output}: {flooded}/{dry} initially dry open-ground cells flooded in view')

if __name__ == '__main__':
    main()
