"""
Plotting utility functions.

"""
import cartopy
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

import warnings

_PLOTTYPE_ARGS = {
    'pcolormesh': dict(linewidth='0'),
    'pcolor': dict(linewidth='0'),
    'contourf': dict(),
}

def _determine_cmap_params(plot_data, vmin=None, vmax=None, cmap=None,
                           center=None, robust=False, extend=None,
                           cnorm=None):
    """
    Use some heuristics to set good defaults for colorbar and range.

    Adapted from xray.
    https://github.com/xray/xray/blob/master/xray/plot/utils.py

    Parameters
    ----------
    plot_data: Numpy array

    Returns
    -------
    cmap_params : dict
    """
    calc_data = plot_data.copy()
    if vmin is None:
        vmin = np.nanmin(calc_data)

    if vmax is None:
        vmax = np.nanmax(calc_data)

    # Simple heuristics for whether these data should  have a divergent map
    divergent = ((vmin < 0) and (vmax > 0)) or center is not None

    # Now set center to 0 so math below makes sense
    if center is None:
        center = 0

    # A divergent map should be symmetric around the center value
    if divergent:
        vlim = max(abs(vmin - center), abs(vmax - center))
        vmin, vmax = -vlim, vlim

    # Now add in the centering value and set the limits
    vmin += center
    vmax += center

    # Choose default colormaps if not provided
    if cmap is None:
        if divergent:
            cmap = "RdBu_r"
        else:
            cmap = "viridis"

    # Allow viridis before matplotlib 1.5
    if mpl.__version__ < 1.5:
        cmap == "Spectral"

    return dict(vmin=vmin, vmax=vmax, cmap=cmap)


def add_colorbar(mappable, fig=None, ax=None, thickness=0.025,
                 shrink=0.1, pad=0.05, orientation='horizontal'):
    """ Add a colorbar into an existing axis or figure. Need to pass
    either an Axis or Figure element to the appropriate keyword
    argument. Should elegantly handle multi-axes figures.

    Parameters
    ----------
    mappable : mappable
        The element set with data to tailor to the colorbar
    fig : Figure
    ax: Axis
    thickness: float
        The width/height of the colorbar in fractional figure area,
        given either vertical/horizontal orientation.
    shrink: float
        Fraction of the width/height of the figure to leave blank
    pad : float
        Padding between bottom/right subplot edge and the colorbar
    orientation : str
        The orientation of the colorbar

    """
    if (fig is None) and (ax is None):
        raise ValueError("Must pass either 'fig' or 'ax'")
    elif fig is None:
        # Plot on Axis
        cb = plt.colorbar(mappable, ax=ax, pad=pad, orientation=orientation)
    else:
        # Plot onto Figure's set of axes
        axes = fig.get_axes()

        # Get coordinates for making the colorbar
        ul = axes[0]
        lr = axes[-1]
        top = ul.get_position().get_points()[1][1]
        bot = lr.get_position().get_points()[0][1]
        right = lr.get_position().get_points()[1][0]
        left = ul.get_position().get_points()[0][0]

        # Calculate colorbar positioning and geometry
        if orientation ==  'vertical':
            cb_left = right + pad
            cb_width = thickness
            cb_bottom = bot + shrink
            cb_height = (top - shrink) - cb_bottom
        elif orientation == 'horizontal':
            cb_left = left + shrink
            cb_width = (right - shrink) - cb_left
            cb_height = thickness
            cb_bottom = (bot - pad) - cb_height
        else:
            raise ValueError("Uknown orientation '%s'" % orientation)

        cax = fig.add_axes([cb_left, cb_bottom,
                            cb_width, cb_height])
        cb = fig.colorbar(mappable, cax=cax, orientation=orientation)

    return cb

def label_ax(ax, cb, cube):
    """ Label the colorbar with the best info available """
    title = cube.name()
    ax.set_title(title)

    cb_label = '{}'.format(cube.units)
    cb.set_label(cb_label)

def geo_plot(cube, ax=None, method='contourf',
             projection='PlateCarree', **kwargs):
    """ Create a global plot of a given variable.

    Parameters:
    -----------
    cube : iris.cube.Cube
        The cube to be plotted.
    ax : axis
        An existing matplotlib axis instance, else one will be created.
    method : str
        String to use for looking up name of plotting function
    projection : str or tuple
        Name of the cartopy projection to use and any args
        necessary for initializing it passed as a dictionary;
        see func:`make_geoaxes` for more information
    **kwargs : dict
        Any additional keyword arguments to pass to the plotter,
        including colormap params. If 'vmin' is not in this
        set of optional keyword arguments, the plot colormap will be
        automatically inferred.
    """

    # Extract longitude and latitude
    xcoord = cube.coords(axis='x')[0:1]
    ycoord = cube.coords(axis='y')[0:1]

    for ixc, iyc in zip(xcoord, ycoord):
        if  'rotated' in ixc.coord_system.grid_mapping_name.lower() \
        and 'rotated' in iyc.coord_system.grid_mapping_name.lower():

            lon2d, lat2d = np.meshgrid(ixc.points, iyc.points)
            nplon = ixc.coord_system.grid_north_pole_longitude
            nplat = ixc.coord_system.grid_north_pole_latitude
            # Unrotate coordinates
            lon, lat = iris.analysis.cartography.unrotate_pole(lon2d, lat2d,
                                                               nplon, nplat)
        else:
            lon, lat = ixc.points, iyc.points

    clon = 0.5*(lon.min()+lon.max())
    clat = 0.5*(lat.min()+lat.max())

    # Set up plotting function
    if method in _PLOTTYPE_ARGS:
        extra_args = _PLOTTYPE_ARGS[method].copy()
    else:
        raise ValueError("Don't know how to deal with '%s' method" % method)
    extra_args.update(**kwargs)

    # Alias a plot function based on the requested method and the
    # datatype being plotted
    plot_func = plt.__dict__[method]

    # `transform` should be the ORIGINAL coordinate system -
    # which is always a simple lat-lon coordinate system in CESM
    # output
    extra_args['transform'] = ccrs.PlateCarree()

    # Was an axis passed to plot on?
    new_axis = ax is None

    if new_axis: # Create a new cartopy axis object for plotting
        # When dealing with regions close to poles,
        # use Lambert Conformal projection
        if (lat > 45).all():
            projection = ('LambertConformal',
                          dict(central_longitude=clon,
                               central_latitude=clat))

        if isinstance(projection, (list, tuple)):
            if len(projection) != 2:
                raise ValueError("Expected 'projection' to only have 2 values")
            projection, proj_kwargs = projection[0], projection[1]
        else:
            proj_kwargs = {}

        # hack to look up the name of the projection in the cartopy
        # reference system namespace; makes life a bit easier, so you
        # can just pass a string with the name of the projection wanted.
        proj = ccrs.__dict__[projection](**proj_kwargs)
        ax = plt.axes(projection=proj)
    else: # Set current axis to one passed as argument
        if not hasattr(ax, 'projection'):
            raise ValueError("Expected `ax` to be a GeoAxes instance")
      #  plt.sca(ax) # Why?

    # Set up map
    ax.coastlines('110m')
    ax.add_feature(cartopy.feature.LAND, facecolor='k', alpha=0.25)
    #ax.add_feature(cartopy.feature.OCEAN)

    try:
        gl = ax.gridlines(crs=extra_args['transform'], draw_labels=True,
                          linewidth=0.5, color='grey', alpha=0.8)
        #LON_TICKS = [ -180, -90, 0, 90, 180 ]
        #LAT_TICKS = [ -90, -60, -30, 0, 30, 60, 90 ]
        gl.xlabels_top   = False
        gl.ylabels_right = False
        #gl.xlines = grid
        #gl.ylines = grid
        #gl.xlocator = mticker.FixedLocator(LON_TICKS)
        #gl.ylocator = mticker.FixedLocator(LAT_TICKS)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    except TypeError:
        warnings.warn("Could not label the given map projection.")

    # Infer colormap settings if not provided
    if not ('vmin' in kwargs):
        #warnings.warn("Re-inferring color parameters...")
        cmap_kws = _determine_cmap_params(cube.data)
        extra_args.update(cmap_kws)

    gp = plot_func(lon, lat, cube.data, **extra_args)

    return ax, gp
