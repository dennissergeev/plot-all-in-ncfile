from __future__ import print_function, division
import argparse
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
#
#import phys_meteo as met
#import um_utils as umu
#import var_utils as var

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

parser.add_argument("-p", "--pr2f", type=str,
                    help="Print to file; if not supplied,\n"
                         "show on display")

if __name__ == '__main__':
    args = parser.parse_args()
    use = 'iris'

    print(args.input_file)
    print(args.variables)

    # Read the given data file
    fn_in = args.input_file
    try:
        ds = iris.load(fn_in)
    except RuntimeError:
        print('Could not open file {}'.format(fn_in))
        sys.exit(1)

    print('File contains:')
    print(ds)

    if not args.variables:
        qplt.contourf(ds[0][0,0,...])
    plt.show()

    #flist = sorted(glob.glob(sys.argv[1]))
    #plid = sys.argv[2]

    #try:
    #    jlev = float(sys.argv[3])
    #except:
    #    jlev = 950.

    lon0 = 11
    lat0 = 75
    mapkw = dict(lon1=lon0-16,lon2=lon0+40,lat1=lat0-9,lat2=lat0+7,tick_incr=[5.,1.],resolution='i',fill=True)
    cmap = 'viridis'
    cbkw = dict(orientation='vertical', shrink=1)

    flist = []

    for ifile in flist:
        print(ifile)
        f = iris.load(ifile)
        umu.replace_unknown_names(f, use='iris')
        umu.unrotate_wind(f, replace=True)

        u_cube = umu.get_cube(f, 'transformed_x_wind')
        v_cube = umu.get_cube(f, 'transformed_y_wind')
        w_cube = umu.get_cube(f, 'tendency') # CHECK!!!!!!!!!!!!!!!!!!
        pv_cube = umu.get_cube(f, 'potential vorticity')
        zg_cube = umu.get_cube(f, 'geopotential_height')
        lon2d, lat2d = [i.points for i in u_cube.aux_coords[-2:]]

        tcoord = f[0].coord('time')
        fcst_dt = cf_units.num2date(tcoord.points, tcoord.units.name, tcoord.units.calendar)
        for n, it in enumerate(tcoord):
            idt = fcst_dt[n]
            print(idt)

            dims = 'zyx'
            if len(tcoord.points) > 1:
                u = u_cube[n, ...]
                v = v_cube[n, ...]
                w = w_cube[n, ...]
                pv = pv_cube[n, ...]
                zg = zg_cube[n, ...]
            else:
                u = u_cube  
                v = v_cube  
                w = w_cube  
                zg = zg_cube
                pv = pv_cube
    
            mcoords = umu.get_model_real_coords(w, dims=dims)
            plev = mcoords[dims.find('z')].points
            ilev = np.argmin(abs(plev[:] - jlev))    
            #
            # Additional diagnostics of wind
            #
            uwind, info3d = tools.prep_data(u.data, dims)
            vwind, _ = tools.prep_data(v.data, dims)
            W2D = standard.WindHorizontal(uwind, vwind, um_res_m, um_res_m)
            vo = W2D.vort_z()
            vo = tools.recover_data(vo, info3d)
    
            #
            # Plotting
            #
            fig, ax = plt.subplots()
            bm = mymap.make_map(ax=ax, **mapkw)
            xx, yy = bm(lon2d, lat2d)
            
            if plid == 'pv':
                c1 = bm.contourf(xx, yy, pv.data[ilev,:,:]*1e6, np.arange(-5,5.5,0.5), cmap=plt.cm.BrBG_r, extend='both')
                c1.cmap.set_under('k')
                c1.cmap.set_over('r')
                cb1 = plt.colorbar(c1, ax=ax)
                cb1.ax.text(0.5, -0.05, r'${0}$'.format('PV units'), va='top', ha='center',size=pp.fntsz_c)
                ax.set_title('Potential vorticity ($x10^6$) at {}hPa'.format(jlev)+\
                             '\n'+idt.strftime('%b-%d %H:%M'),size=pp.fntsz_t)
                
                subsubdir = 'um_pv'
                imgname = 'um_pv_{0}hpa_'.format(jlev) + idt.strftime('%y%m%d-%H%M')
