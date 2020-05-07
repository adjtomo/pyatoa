#!/usr/bin/env python3
"""
Basemap map making tools

Functions used to produce maps using Basemap that have a standard look across
the Pyatoa workflow. Map creation functions are located outside this toolbox
"""
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from scipy.interpolate import griddata

from pyatoa.utils.srcrcv import gcd_and_baz
from pyatoa.utils.calculate import myround


def legend():
    """
    Create a legend on the current axes with astandard look
    """
    leg = plt.legend(loc="lower right")
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(1)


def place_scalebar(m, map_corners, **kwargs):
    """
    Put the scale bar in a corner at a reasonable distance from each edge

    Handy reminder for moving the scalebar around:
        latitude is up, down
        longitude is right, left

    :type m: Basemap
    :param m: basemap object
    :type map_corners: dict of floats
    :param map_corners: [lat_bot,lat_top,lon_left,lon_right]
    :type loc: str
    :param loc: location of scalebar, 'upper-right' or 'lower-right'
    """
    loc = kwargs.get("loc", "upper-right")

    mc = map_corners
    if loc == "upper-right":
        latscale = mc['lat_min'] + (mc['lat_max'] - mc['lat_min']) * 0.94
        lonscale = mc['lon_min'] + (mc['lon_max'] - mc['lon_min']) * 0.875
    if loc == "lower-right":
        latscale = mc['lat_min'] + (mc['lat_max'] - mc['lat_min']) * 0.04
        lonscale = mc['lon_min'] + (mc['lon_max'] - mc['lon_min']) * 0.9
    m.drawmapscale(lonscale, latscale, lonscale, latscale, 100,
                   yoffset=0.01 * (m.ymax-m.ymin), zorder=5000, linewidth=2,
                   fontsize=13
                   )


def build_colormap(array, cmap=cm.jet_r):
    """
    Build a custom range colormap. Round values before.

    :type array: numpy.array
    :param array: array to build colormap from
    :type cmap: matplotlib colormap
    :param cmap: colormap to create a custom range for
    :rtype colormap: matplotlib.cm.ScalarMappable
    :return colormap: custom colormap
    """
    vmax = myround(np.nanmax(array), base=1, choice='up')
    vmin = myround(np.nanmin(array), base=1, choice='down')
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    colormap = cm.ScalarMappable(norm=norm, cmap=cmap)

    return colormap


def event_beachball(m, event, fm_type="focal_mechanism", **kwargs):
    """
    Plot event beachball for a given moment tensor attribute from event object.

    Note:
    if strike_dip_rake chosen, nodal plane 1 is used for focal mechanism,
    assuming that is the preferred plane

    :type m: Basemap
    :param m: basemap object
    :type fm_type: str
    :param fm_type: focal mech. type, 'focal_mechanism' or 'strike_dip_rake'
    :type event: obspy.core.event.Event
    :param event: event object which should contain focal mechanism
    """
    width = kwargs.get("width", 2.6E4)
    facecolor = kwargs.get("facecolor", 'r')
    linewidth = kwargs.get("linewidth", 1)
    zorder = kwargs.get("zorder", 1000)

    eventx, eventy = m(event.preferred_origin().longitude,
                       event.preferred_origin().latitude
                       )

    # No focal mechanism? Just plot a ploint
    if not hasattr(event, 'focal_mechanisms'):
        m.scatter(eventx, eventy, marker="o", color=facecolor,
                  edgecolor="k", s=105, zorder=zorder, linewidth=1.75)
    else:
        if fm_type == "focal_mechanism":
            fm = event.focal_mechanisms[0].moment_tensor.tensor or \
                 event.preferred_focal_mechanism().moment_tensor.tensor
            beach_input = [fm['m_rr'], fm['m_tt'], fm['m_pp'],
                           fm['m_rt'], fm['m_rp'], fm['m_tp']
                           ]
        elif fm_type == "strike_dip_rake":
            nod_plane = event.focal_mechanisms[0].nodal_planes or \
                        event.preferred_focal_mechanism().nodal_planes
            # try determine the preferred nodal plane, default to 1
            try:
                sdr = nod_plane[f"nodal_plane_{nod_plane.preferred_plane}"]
            except AttributeError:
                sdr = nod_plane.nodal_plane_1
            beach_input = [sdr.strike, sdr.dip, sdr.rake]

        b = beach(beach_input, xy=(eventx, eventy), width=width,
                  linewidth=linewidth, facecolor=facecolor)
        b.set_zorder(zorder)
        ax = plt.gca()
        ax.add_collection(b)


def plot_stations(m, inv, event=None, annotate_names=False,
                  color_by_network=False):
    """
    Fill map with stations based on station availability and network

    :type m: Basemap
    :param m: basemap object
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing relevant network and stations
    :type event: obspy.core.event.Event
    :param event: event object which should contain focal mechanism
    :type annotate_names: bool
    :param annotate_names: whether or not station names should placed nearby
    :type color_by_network: bool
    :param color_by_network: decided the coloring of different networks
    """
    # Sort the networks by station name
    code, num_sta = [], []
    for net in inv:
        code.append(net.code)
        num_sta.append(len(net))
    num_sta, network_codes = zip(*sorted(zip(num_sta, code)))

    # A simple list of colors to color by network.
    # NOTE: Assuming we never plot more networks than length of available colors
    if color_by_network:
        available_colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
    else:
        available_colors = ['None'] * len(network_codes)

    # Plot stations
    for i, net_code in enumerate(network_codes):
        net = inv.select(net_code)[0]
        # For legend
        m.scatter(0, 0, marker='v', color=available_colors[i],
                  linestyle='-', s=80, linewidth=1.25, zorder=1,
                  label=code)
        for sta in net:
            # If an event is given, check that the station is active at time
            if event and not sta.is_active(time=event.preferred_origin().time):
                continue
            x, y = m(sta.longitude, sta.latitude)
            m.scatter(x, y, marker='v', color=available_colors[i],
                      edgecolor='k', linestyle='-', s=80, linewidth=1.25,
                      zorder=1001
                      )
            # Annotate the station name next to the station
            if annotate_names:
                plt.annotate(f"{net.code}.{sta.code}", xy=(x, y), xytext=(x, y),
                             zorder=6, fontsize=7, bbox=dict(facecolor='w',
                                                             edgecolor='k',
                                                             boxstyle='round')
                             )


def plot_stations_simple(m, lats, lons):
    """
    Fill map with stations based on latitude longitude values, nothing fancy

    :type m: Basemap
    :param m: basemap object
    :type lats: np.ndarray
    :param lats: array of latitude values
    :type lons: np.ndarray
    :param lons: array of longitude values
    """
    x, y = m(lons, lats)
    m.scatter(x, y, marker='v', color='None', edgecolor='k', linestyle='-',
              s=95, linewidth=1.75, zorder=100
              )


def annotate_srcrcv_information(m, event, inv):
    """
    Annotate event receiver information into bottom right corner of the map

    :type m: Basemap
    :param m: basemap object
    :type event: obspy.core.event.Event
    :param event: event object
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing relevant network and stations
    """
    event.origins[0].time.precision = 0
    gcdist, baz = gcd_and_baz(event, inv[0][0])
    event_id = event.resource_id.id.split('/')[1]

    plt.annotate(s=(f"{event_id} / {inv[0].code}.{inv[0][0].code}\n"
                    f"{event.preferred_origin().time}\n"
                    f"{event.preferred_magnitude().magnitude_type}="
                    f"{event.preferred_magnitude().mag:.2f}\n"
                    f"Depth(km)={event.preferred_origin().depth*1E-3:.2f}\n"
                    f"Dist(km)={gcdist:.2f}\n"
                    f"BAz(deg)={baz:.2f}"),
                 xy=(m.xmin + (m.xmax-m.xmin) * 0.675,
                     m.ymin + (m.ymax-m.ymin) * 0.035), multialignment='right',
                 fontsize=10
                 )


def connect_source_receiver(m, event, sta, **kwargs):
    """
    Draw a dashed line connecting the station and receiver, highlight station
    with a colored marker. Useful for visualizing a source-receiver pair.

    :type m: Basemap
    :param m: basemap object
    :type event: obspy.core.event.Event
    :param event: event object
    :type sta: obspy.core.inventory.Inventory
    :param sta: inventory containing relevant network and stations
    """
    linestyle = kwargs.get("linestyle", "--")
    linewidth = kwargs.get("linewidth", 1.1)
    linecolor = kwargs.get("color", "k")
    marker = kwargs.get("marker", "v")
    markercolor = kwargs.get("markercolor", "r")
    zorder = kwargs.get("zorder", 100)

    # Get coordinates in map extent
    event_x, event_y = m(event.preferred_origin().longitude,
                         event.preferred_origin().latitude)
    station_x, station_y = m(sta.longitude, sta.latitude)

    # Plot line, station, event
    m.plot([event_x, station_x], [event_y, station_y], linestyle, linewidth,
           c=linecolor, zorder=zorder-10)
    m.scatter(station_x, station_y, marker=marker, color=markercolor,
              edgecolor='k', linestyle='-', s=75, zorder=zorder)
    m.scatter(event_x, event_y, marker="o", color=markercolor, edgecolor="k",
              s=105, zorder=zorder, linewidth=1.75)


def interpolate_and_contour(m, x, y, z, len_xi, len_yi, fill=True,
                            cbar_label='', colormap='viridis',
                            interpolation_method='cubic', marker='o'):
    """
    Interpolate over scatter points on a regular grid, and create a contour plot

    Station locations are irregular over the map and it's sometimes difficult to
    visualize data based soley on the colors of the station markers. This
    function will interpolate the station values (e.g. misfit) over the area
    covered by the stations, and create a contour plot that can be overlaid onto
    the map object

    :type m: Basemap
    :param m: basemap object
    :type x: list
    :param x: x coordinates
    :type y: list
    :param y: y coordinates
    :type z: list
    :param z: values to interpolate
    :type len_xi: int
    :param len_xi: x-axis length of the regular grid to interpolate onto
    :type len_yi: int
    :param len_yi: y-axis length of the regular grid to interpolate onto
    :type fill: bool
    :param fill: contourf or countour
    :type cbar_label: str
    :param cbar_label: custom label for the colorbar
    :type colormap: str
    :param colormap: matplotlib colorscale to be shared between contour, scatter
    :type interpolation_method: str
    :param interpolation_method: 'linear', 'nearest', 'cubic', passed to
        scipy.interpolate.griddata
    :type marker: str
    :param marker: marker type for the station scatterplot
    """
    # Create a grid of data to interpolate data onto
    xi, yi = np.mgrid[min(x):max(y):len_xi,
                      min(x):max(y):len_yi]

    # np meshgrid can also be used, but it gives a jagged edge to contourf
    # xi, yi = np.meshgrid(np.linspace(min(x), max(y), len_xi),
    #                      np.linspace(min(x), max(y), len_yi)
    #                      )

    # Interpolate the data and plot as a contour overlay
    zi = griddata((x, y), z, (xi, yi), method=interpolation_method)

    # Use contourf, for filled contours
    if fill:
        m.contourf(xi, yi, zi, vmin=0, zorder=100, alpha=0.6, cmap=colormap)
        # Add a colorbar
        max_value = myround(np.nan_to_num(zi).max(), base=10, choice='up')
        cmap = plt.cm.ScalarMappable(cmap=colormap)
        cmap.set_array(zi)
        cmap.set_clim(0., max_value)
        cbar = plt.colorbar(cmap, boundaries=np.arange(0, max_value+10, 10),
                            shrink=0.8, extend='min')
        cbar.set_label(cbar_label, rotation=270)

    # Use contour for only lines, places values inline with contour
    else:
        cs = m.contour(xi, yi, zi, vmin=0, zorder=100, alpha=0.6, cmap=colormap)
        # Label the contour values
        plt.clabel(cs, cs.levels, fontsize=8, inline=True, fmt='%.1E')

    # Plot the points that were used in the interpolation
    m.scatter(x, y, c=z, alpha=1., edgecolor='k', linestyle='-', s=100,
              cmap=colormap, linewidth=2, zorder=101, marker=marker
              )


def initiate_basemap(map_corners=None, scalebar=True, **kwargs):
    """
    Set up the basemap object in the same way each time

    :type map_corners: dict of floats
    :param map_corners: {lat_min,lat_max,lon_min,lon_max}
    :type scalebar: bool
    :param scalebar: add a scalebar to the map
    :rtype m: Basemap
    :return m: basemap object
    """
    continent_color = kwargs.get("contininent_color", "w")
    lake_color = kwargs.get("lake_color", "w")
    coastline_zorder = kwargs.get("coastline_zorder", 5)
    coastline_linewidth = kwargs.get("coastline_linewidth", 2.0)
    fill_color = kwargs.get("fill_color", "w")
    scalebar_location = kwargs.get("scalebar_location", "upper-right")
    latlon_linewidth = kwargs.get("latlon_linewidth", 0.)   

    # Initiate map and draw in style. Stereographic projection if regional
    # corners are given, otherwise world map
    if map_corners:
        m = Basemap(projection='stere', resolution='h', rsphere=6371200,
                    lat_0=(map_corners['lat_min'] + map_corners['lat_max'])/2,
                    lon_0=(map_corners['lon_min'] + map_corners['lon_max'])/2,
                    llcrnrlat=map_corners['lat_min'],
                    urcrnrlat=map_corners['lat_max'],
                    llcrnrlon=map_corners['lon_min'],
                    urcrnrlon=map_corners['lon_max'],
                    )
        m.drawparallels(np.arange(int(map_corners['lat_min']),
                                  int(map_corners['lat_max']), 1),
                        labels=[1, 0, 0, 0], linewidth=latlon_linewidth,
                        )
        m.drawmeridians(np.arange(int(map_corners['lon_min']),
                                  int(map_corners['lon_max']) + 1, 1),
                        labels=[0, 0, 0, 1], linewidth=latlon_linewidth,
                        )
    else:
        m = Basemap(projection='cyl', resolution='c', llcrnrlat=-90,
                    urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180,
                    )

    m.drawcoastlines(linewidth=coastline_linewidth, zorder=coastline_zorder)
    m.fillcontinents(color=continent_color, lake_color=lake_color)
    m.drawmapboundary(fill_color=fill_color)

    if scalebar and map_corners:
        place_scalebar(m, map_corners, loc=scalebar_location)

    return m

