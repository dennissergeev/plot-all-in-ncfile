from __future__ import print_function, division, absolute_import
import argparse
import cartopy.crs as ccrs
import cf_units
import datetime
import glob
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import numpy as np
import os
import sys

import plot_util

import um_utils as umu

DESCR = """
Plot snapshots from subsets of variables from a UM output.
"""

parser = argparse.ArgumentParser(os.path.basename(__file__),
                                 description=DESCR,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("input_file", nargs='*', help="File to extract data from (netCDF or PP)")
parser.add_argument("-v", "--variables", nargs="*",
                    help="Variables to plot; if not supplied,\n"
                         "will plot the first variable in file")

parser.add_argument("-lev", "--level", default=0, type=float,
                    help="Level to plot\n"
                         "(if the variable has a vertical dimension)")

parser.add_argument("-levname", "--level_name", type=str,
                    help="Pressure level to plot\n"
                         "(if the variable's vertical dimension is pressure)")

parser.add_argument("-p", "--pr2f", type=str,
                    help="Print to file; if not supplied,\n"
                         "show on display")

parser.add_argument("--unrotatewind", default=True, type=bool,
                    help="If the file contains u- and v-winds\n"
                         "on a rotated pole coordinate system,\n"
                         "rotate them to a 'normal' coordinate system")

parser.add_argument("--use_interface", default='iris', type=str,
                    help="Interface to use (only iris is operational)")

if __name__ == '__main__':
    args = parser.parse_args()
    use = args.use_interface

    #print(args)

    # Read the given data file
    fn_in = args.input_file
    try:
        ds = iris.load(fn_in)
    except IOError:
        print('Could not open file {}'.format(fn_in))
        sys.exit(1)

    umu.replace_unknown_names(ds, use=use)

    print('File contains:')
    print(ds)

    if args.unrotatewind:
        umu.unrotate_wind(ds, replace=True)

    if args.variables:
        varnames = args.variables
    else:
        varnames = [ds[0].name()]

    print(ds)
    for ivarname in varnames:
        print('Plotting {}...'.format(ivarname))
        icube = umu.get_cube(ds, ivarname)

        zcoord = None
        if len(icube.coords(axis='z')) > 0:
            if args.level_name:
                for icoord in icube.coords(axis='z'):
                    if args.level_name.lower() in icoord.name().lower():
                        zcoord = icoord
                if zcoord is None:
                    raise ValueError('no z-coord matches {}'.format(args.level_name))
            else:
                zcoord = icube.coords(axis='z')[0]
            ilev = np.argmin(abs(args.level - zcoord.points))
        plot_data = icube[ilev,...]
        # Plot this timeslice
        fig = plt.figure()
        ax = None
        #ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        var_kwargs = dict()
        ax, plot = plot_util.geo_plot(plot_data, ax=ax, **var_kwargs)
        cb = plot_util.add_colorbar(plot, ax=ax, orientation='vertical')
        plot_util.label_ax(ax, cb, icube)

        #tcoord = cube.coord('time')
        #um_dt = cf_units.num2date(tcoord.points, tcoord.units.name, tcoord.units.calendar)
    plt.show()
