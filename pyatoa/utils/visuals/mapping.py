#!/usr/bin/env python3
"""
Map making functionalities

Produces a basemap of target region with beachball representing event,
all available stations with data for the given origin time, relevant station
highlighted, connecting line between station and event, and important
information annotated (e.g. misift information, distance, BAz etc.)
"""
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from scipy.interpolate import griddata

from pyatoa.utils.operations.source_receiver import gcd_and_baz
from pyatoa.utils.operations.calculations import myround
from pyatoa.utils.asdf import extractions
from pyatoa.utils.visuals import map_extras


def legend():
    leg = plt.legend(loc="lower right")
    leg.get_frame().set_edgecolor('k')
    leg.get_frame().set_linewidth(1)


def standalone_colorbar(bounds, steps):
    """
    generate a colorbar as a new figure
    TO DO: clean this up and potentially add it into basemap
    """
    fig = plt.figure(figsize=(0.5, 8))
    mpl.rcParams['font.size'] = 25
    mpl.rcParams['axes.linewidth'] = 4

    ax = fig.add_axes([0.2,0.05,0.4,0.75])
    # bounds = range(0,90,10)
    cmap = mpl.cm.get_cmap('jet_r')
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, extend='max',
                                   orientation='vertical')
    cb.set_label("Event Depth (km)")
    cb.ax.invert_yaxis()
    plt.show()


def place_scalebar(m, map_corners, loc="upper-right"):
    """
    Put the scale bar in the corner at a reasonable distance from each edge
    Handy reminder:
        latitude is up, down
        longitude is right, left
        ;)

    :type m: Basemap
    :param m: basemap object
    :type map_corners: dict of floats
    :param map_corners: [lat_bot,lat_top,lon_left,lon_right]
    :type loc: str
    :param loc: location of scalebar, 'upper-right' or 'lower-right'
    """
    mc = map_corners
    if loc == "upper-right":
        latscale = mc['lat_min'] + (mc['lat_max'] - mc['lat_min']) * 0.94
        lonscale = mc['lon_min'] + (mc['lon_max'] - mc['lon_min']) * 0.875
    if loc == "lower-right":
        latscale = mc['lat_min'] + (mc['lat_max'] - mc['lat_min']) * 0.24
        lonscale = mc['lon_min'] + (mc['lon_max'] - mc['lon_min']) * 0.9
    m.drawmapscale(lonscale, latscale, lonscale, latscale, 100,
                   yoffset=0.01 * (m.ymax-m.ymin), zorder=5000, linewidth=2,
                   fontsize=13
                   )


def build_colormap(array):
    """
    Build a custom range colormap, hardcoded colormap. Round values before.

    :type array: numpy.array
    :param array: array to build colormap from
    :rtype colormap: matplotlib.cm.ScalarMappable
    :return colormap: custom colormap
    """
    vmax = myround(np.nanmax(array), base=1, choice='up')
    vmin = myround(np.nanmin(array), base=1, choice='down')
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.jet_r
    colormap = cm.ScalarMappable(norm=norm, cmap=cmap)

    return colormap


def scale_magnitude(magitude):
    """
    Short function to standardize magnitude scaling on plots
    :param magitude: float
    :return:
    """
    raise NotImplementedError


def event_beachball(m, event, fm_type="focal_mechanism", **kwargs):
    """
    Plot event beachball for a given geonet moment tensor list,
    read in the from the GeoNet moment tensor csv file.

    :type m: Basemap
    :param m: basemap object
    :type fm_type: str
    :param fm_type: focal mechanisms type,
        'focal_mechanism' or 'strike_dip_rake'
    :type event: obspy.core.event.Event
    :param event: event object which should contain focal mechanism
    """
    width = kwargs.get("width", 2.6E4)
    facecolor = kwargs.get("facecolor", 'r')
    zorder = kwargs.get("zorder", 1000)

    eventx, eventy = m(event.preferred_origin().longitude,
                       event.preferred_origin().latitude
                       )

    # No focal mechanism? Just plot a ploint
    if not hasattr(event, 'focal_mechanisms'):
        m.scatter(eventx, eventy, marker="o", color=facecolor,
                  edgecolor="k", s=105, zorder=zorder, linewidth=1.75)

    if fm_type == "focal_mechanism":
        beach_input = [
            event.focal_mechanisms[0].moment_tensor.tensor['m_rr'],
            event.focal_mechanisms[0].moment_tensor.tensor['m_tt'],
            event.focal_mechanisms[0].moment_tensor.tensor['m_pp'],
            event.focal_mechanisms[0].moment_tensor.tensor['m_rt'],
            event.focal_mechanisms[0].moment_tensor.tensor['m_rp'],
            event.focal_mechanisms[0].moment_tensor.tensor['m_tp']
        ]
    elif fm_type == "strike_dip_rake":
        beach_input = [
            event.focal_mechanisms[0].nodal_planes.nodal_plane_1.strike,
            event.focal_mechanisms[0].nodal_planes.nodal_plane_1.dip,
            event.focal_mechanisms[0].nodal_planes.nodal_plane_1.rake
        ]
    b = beach(beach_input, xy=(eventx, eventy), width=width, linewidth=1,
              facecolor=facecolor)
    b.set_zorder(zorder)
    ax = plt.gca()
    ax.add_collection(b)


def plot_stations(m, inv, event=None, **kwargs):
    """
    Fill map with stations based on station availability and network
    TO DO: fix this up and comment it nicely

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
    annotate_names = kwargs.get("annotate_names", False)
    color_by_network = kwargs.get("color_by_network", False)

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
        net = inv.select(net_code)
        # For legend
        m.scatter(0, 0, marker='v', color=available_colors[i],
                  linestyle='-', s=80, linewidth=1.25, zorder=1,
                  label=code
                  )
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
                plt.annotate("{n}.{s}".format(net.code, sta.code),
                             xy=(x,y), xytext=(x,y), zorder=6, fontsize=7,
                             bbox=dict(facecolor='w', edgecolor='k',
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
    Annotate event receiver information into hard coded map area

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

    plt.annotate(
        s=("{id} / {net}.{sta}\n"
           "{date}\n"
           "{type}={mag:.2f}\n"
           "Depth(km)={depth:.2f}\n"
           "Dist(km)={dist:.2f}\n"
           "BAz(deg)={baz:.2f}").format(
            id=event_id, net=inv[0].code, sta=inv[0][0].code,
            date=event.preferred_origin().time,
            depth=event.preferred_origin().depth*1E-3,
            type=event.preferred_magnitude().magnitude_type,
            mag=event.preferred_magnitude().mag, dist=gcdist, baz=baz
        ),
        xy=(m.xmin + (m.xmax-m.xmin) * 0.675,
            m.ymin + (m.ymax-m.ymin) * 0.035),
        multialignment='right', fontsize=10
    )


def connect_source_receiver(m, event, sta, **kwargs):
    """
    draw a dashed line connecting the station and receiver, highlight station

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

    event_x, event_y = m(event.preferred_origin().longitude,
                         event.preferred_origin().latitude)
    station_x, station_y = m(sta.longitude, sta.latitude)
    m.plot([event_x, station_x], [event_y, station_y], linestyle, linewidth,
           c=linecolor, zorder=zorder-10)
    m.scatter(station_x, station_y, marker=marker, color=markercolor,
              edgecolor='k', linestyle='-', s=75, zorder=zorder)
    m.scatter(event_x, event_y, marker="o", color=markercolor, edgecolor="k",
              s=105, zorder=zorder, linewidth=1.75)


def interpolate_and_contour(m, x, y, z, len_xi, len_yi, fill=True,
                            colormap='viridis', interpolation_method='cubic'):
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
    :type colormap: str
    :param colormap: matplotlib colorscale to be shared between contour, scatter
    :type interpolation_method: str
    :param interpolation_method: 'linear', 'nearest', 'cubic', passed to
        scipy.interpolate.griddata
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
        cs = m.contourf(xi, yi, zi, vmin=0, zorder=100, alpha=0.6,
                        cmap=colormap)
        # Add a colorbar
        cbar = m.colorbar(cs, location='right', pad="5%")
        cbar.set_label("misfit")
    # Use contour for only lines, places values inline with contour
    else:
        cs = m.contour(xi, yi, zi, vmin=0, zorder=100, alpha=0.6, cmap=colormap)
        # Label the contour values
        plt.clabel(cs, cs.levels, fontsize=8, inline=True, fmt='%.1E')

    # Plot the points that were used in the interpolation
    m.scatter(x, y, c=z, alpha=0.75, edgecolor='k', linestyle='-', s=100,
              cmap=colormap, linewidth=2, zorder=101
              )


def initiate_basemap(map_corners, scalebar=True, **kwargs):
    """
    set up the basemap object in the same way each time

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
    coastline_linewidth = kwargs.get("coastline_linewidth", 1.0)

    # Initiate map and draw in style
    m = Basemap(projection='stere', resolution='h', rsphere=6371200,
                lat_0=(map_corners['lat_min'] + map_corners['lat_max'])/2,
                lon_0=(map_corners['lon_min'] + map_corners['lon_max'])/2,
                llcrnrlat=map_corners['lat_min'],
                urcrnrlat=map_corners['lat_max'],
                llcrnrlon=map_corners['lon_min'],
                urcrnrlon=map_corners['lon_max'],
                )
    m.drawcoastlines(linewidth=coastline_linewidth, zorder=coastline_zorder)
    m.fillcontinents(color=continent_color, lake_color=lake_color)
    m.drawparallels(np.arange(int(map_corners['lat_min']),
                              int(map_corners['lat_max']), 1),
                    labels=[1, 0, 0, 0], linewidth=0.0)
    m.drawmeridians(np.arange(int(map_corners['lon_min']),
                              int(map_corners['lon_max'])+1, 1),
                    labels=[0, 0, 0, 1], linewidth=0.0)

    if scalebar:
        place_scalebar(m, map_corners)

    return m


def standalone_map(map_corners, inv, catalog, annotate_names=False,
                   show_nz_faults=False, color_by_network=False,
                   figsize=(10,9.4), dpi=100, show=True, save=None):
    """
    To be used in a standalone mapmaker. Plots a catalog and an inventory
    to show all events and all stations for e.g. a given tomographic inversino

    :type map_corners: dict of floats
    :param map_corners: [lat_bot,lat_top,lon_left,lon_right]
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing relevant network and stations
    :type catalog: obspy.core.event.catalog.Catalog
    :param catalog: a catalog of events to plot
    :type map_corners: dict of floats
    :param map_corners: {lat_min,lat_max,lon_min,lon_max}
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing relevant network and stations
    :type annotate_names: bool
    :param annotate_names: annotate station names next to markers
    :type show_nz_faults: bool
    :param show_nz_faults: call hardcoded fault plotting scripts (TO DO change)
    :type color_by_network: bool
    :param color_by_network: color station markers based on network name
    :type figsize: tuple
    :param figsize: size of the figure
    :type dpi: int
    :param dpi: dots per inch for resolution
    :type show: bool
    :param show: show the plot once generated, defaults to False
    :type save: str
    :param save: absolute filepath and filename if figure should be saved
    :rtype f: matplotlib figure
    :return f: figure object
    :rtype m: Basemap
    :return m: basemap object
    """
    # Initiate matplotlib instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    m = initiate_basemap(map_corners=map_corners, scalebar=True)

    # TO DO: remove hard coding
    # Plot fault lines, hardcoded into structure
    if show_nz_faults:
        map_extras.plot_hikurangi_trench(m)
        map_extras.plot_active_faults(m)

    # If given, plot all background stations for this given event.
    if inv:
        plot_stations(m, inv=inv, event=None,
                      annotate_names=annotate_names,
                      color_by_network=color_by_network
                      )

    # TO DO: scale magnitudes
    # scaled_magnitudes = scale_magnitude(catalog)

    # If an event is given, try to plot the focal-mechanism beachball
    for event in catalog:
        event_beachball(m, event)

    # Finality
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show:
        if show == "hold":
            return f, m
        else:
            plt.show()
    plt.close()

    return f, m


def event_misfit_map(map_corners, ds, model, step=None,
                     annotate_station_info=False, contour_overlay=True,
                     show_nz_faults=False, figsize=(10, 9.4), dpi=100,
                     show=True, save=None):
    """
    To be used to plot misfit information from a pyasdf Dataset

    :type map_corners: dict of floats
    :param map_corners: [lat_bot,lat_top,lon_left,lon_right]
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing things to plot
    :type model: str
    :param model: model number, e.g. 'm00'
    :type step: str
    :param step: step number, e.g. 's00', used for title plotting only
    :type map_corners: dict of floats
    :param map_corners: {lat_min,lat_max,lon_min,lon_max}
    :type annotate_station_info: bool or str
    :param annotate_station_info: annotate station names and info next to marker
        can also be 'simple'
    :type contour_overlay: bool
    :param contour_overlay: interpolate z_values and create contour
    :type show_nz_faults: bool
    :param show_nz_faults: call hardcoded fault plotting scripts (TO DO change)
    :type figsize: tuple
    :param figsize: size of the figure
    :type dpi: int
    :param dpi: dots per inch for resolution
    :type show: bool
    :param show: show the plot once generated, defaults to False
    :type save: str
    :param save: absolute filepath and filename if figure should be saved
    :rtype f: matplotlib figure
    :return f: figure object
    :rtype m: Basemap
    :return m: basemap object
    """
    # Initiate matplotlib instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    m = initiate_basemap(map_corners=map_corners, scalebar=True)

    # If catalog in dataset, try plot the focal-mechanism beachball
    event = None
    if hasattr(ds, 'events'):
        event = ds.events[0]
        event.origins[0].time.precision = 0
        event_id = event.resource_id.id.split('/')[-1]
        event_beachball(m, event, zorder=102)
        plt.annotate(
            s=("{id}\n"
               "{date}\n"
               "{type}={mag:.2f}\n"
               "Depth(km)={depth:.2f}\n").format(
                id=event_id, date=event.preferred_origin().time,
                depth=event.preferred_origin().depth * 1E-3,
                type=event.preferred_magnitude().magnitude_type,
                mag=event.preferred_magnitude().mag,
            ),
            xy=(m.xmax - (m.xmax - m.xmin) * 0.3,
                m.ymin + (m.ymax - m.ymin) * 0.01),
            multialignment='right', fontsize=10
        )

    # Get some event-wide values
    if hasattr(ds, 'auxiliary_data'):
        if hasattr(ds.auxiliary_data, 'MisfitWindows'):
            window_dict = extractions.count_misfit_windows(
                ds, model, count_by_stations=True
            )

    # If waveforms in dataset, plot the stations
    x_values, y_values, z_values = [],[],[]
    if hasattr(ds, 'waveforms'):
        for sta in ds.waveforms.list():
            sta_anno = ""
            if 'StationXML' in ds.waveforms[sta].list():
                sta_anno += "{}\n".format(sta)
                # Get distance and backazimuth between event and station
                if event is not None:
                    gcdist, baz = gcd_and_baz(event,
                                              ds.waveforms[sta].StationXML[0][0]
                                              )
                    sta_anno += "{dist:.1f}km\n{baz:.1f}deg".format(
                        dist=gcdist, baz=baz
                    )
                # If auxiliary data is given, add some extra information to anno
                if hasattr(ds, 'auxiliary_data'):
                    # Determine the total number of windows for the given model
                    if 'MisfitWindows' in ds.auxiliary_data.list():
                        try:
                            num_windows = window_dict[sta]
                        except KeyError:
                            num_windows = 0
                        sta_anno += "\n{}".format(num_windows)

                    # Determine the total misfit for a given model
                    if 'AdjointSources' in ds.auxiliary_data:
                        total_misfit = extractions.sum_misfits(
                            ds, model, station=sta
                        )
                        if total_misfit:
                            sta_anno += "/{:.2E}".format(total_misfit)
                        else:
                            total_misfit = 0

                # Determine the station coordinates on the map object
                sta_x, sta_y = m(
                    ds.waveforms[sta].StationXML[0][0].longitude,
                    ds.waveforms[sta].StationXML[0][0].latitude
                )
                x_values.append(sta_x)
                y_values.append(sta_y)
                z_values.append(total_misfit)

                # Plot waveforms with data differently than those without
                if 'observed' in ds.waveforms[sta].get_waveform_tags():
                    markersize = 75
                    linewidth = 1.75
                    alpha = 1.0
                else:
                    markersize = 50
                    linewidth = 1.0
                    alpha = 0.5

            # Plot station locations
            m.scatter(sta_x, sta_y, marker='v', color='None', alpha=alpha,
                      edgecolor='k', linestyle='-', s=markersize,
                      linewidth=linewidth, zorder=100)

            # Annotate station information
            if annotate_station_info:
                # Just annotate the station names
                if annotate_station_info == "simple":
                    plt.annotate(s=sta_anno.split('\n')[0],
                                 xy=(sta_x - (m.xmax - m.xmin) * 0.01,
                                     sta_y + (m.ymax - m.ymin) * 0.01),
                                 fontsize=10, zorder=103,
                                 )
                # Annotate as much information as possible, use background box
                else:
                    plt.annotate(s=sta_anno,
                                 xy=(sta_x - (m.xmax - m.xmin) * 0.01,
                                     sta_y - (m.ymax - m.ymin) * 0.0025),
                                 bbox=dict(boxstyle="round", fc="white",
                                           ec="k", lw=1.25, alpha=0.75),
                                 multialignment='left', fontsize=6,
                                 zorder=103
                                 )

    # Plot a contour of a given Z value based on receiver locations
    # Fill can be changed to make the contour just lines without a colorbar
    if contour_overlay:
        interpolate_and_contour(m=m, x=x_values, y=y_values, z=z_values,
                                len_xi=100, len_yi=150, colormap='viridis',
                                interpolation_method='cubic', fill=True,
                                )

    # Plot fault lines, hardcoded into structure
    if show_nz_faults:
        map_extras.plot_hikurangi_trench(m)
        map_extras.plot_active_faults(m)

    # Set a title
    if step:
        title = "{m}.{s} misfit map, {r} stations".format(
            m=model, s=step, r=len(x_values)
        )
        if "Statistics" in ds.auxiliary_data.list():
            if model in ds.auxiliary_data.Statistics.list():
                if step in ds.auxiliary_data.Statistics[model].list():
                    pars = ds.auxiliary_data.Statistics[model][step].parameters
                    title += (
                        "\navg misfit: {a:.2f}, num windows: {w}\n"
                        "max misfit comp: {m1:.2f}/{m2}\n " 
                        "min misfit comp: {m3:.2f}/{m4}"
                       .format(
                            a=pars["average_misfit"], m1=pars["max_misfit"],
                            m2=pars["max_misfit_component"],
                            m3=pars["min_misfit"],
                            m4=pars["min_misfit_component"],
                            w=pars["number_misfit_windows"])
                    )

    else:
        title = "{m} misfit map, {r} stations".format(m=model,
                                                            r=len(x_values))
    plt.title(title)

    # Finality
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show:
        if show == "hold":
            return f, m
        else:
            plt.show()
    plt.close()

    return f, m


def manager_map(map_corners, inv=None, event=None, stations=None,
                annotate_names=False, show_nz_faults=False,
                color_by_network=False, figsize=(10, 9.4), dpi=100,
                show=True, save=None):
    """
    Initiate and populate a basemap object for New Zealands north island.
    Functionality to manually ignore stations based on user quality control
    Takes station coordinates and coloring from npz files
    Choice to annotate two stations which correspond to waveforms
    Calls beachball and trench tracer to populate basemap

    :type map_corners: dict of floats
    :param map_corners: {lat_min,lat_max,lon_min,lon_max}
    :type event: obspy.core.event.Event
    :param event: event object
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing relevant network and stations
    :type stations: obspy.core.inventory.Inventory or numpy.ndarray
    :param stations: background stations to plot on the map that will not
        interact with the source. if given as an inventory object, plotting
        will be more complex.
    :type annotate_names: bool
    :param annotate_names: annotate station names next to markers
    :type show_nz_faults: bool
    :param show_nz_faults: call hardcoded fault plotting scripts (TO DO change)
    :type color_by_network: bool
    :param color_by_network: color station markers based on network name
    :type figsize: tuple
    :param figsize: size of the figure
    :type dpi: int
    :param dpi: dots per inch for resolution
    :type show: bool
    :param show: show the plot once generated, defaults to False
    :type save: str
    :param save: absolute filepath and filename if figure should be saved
    :rtype f: matplotlib figure
    :return f: figure object
    :rtype m: Basemap
    :return m: basemap object
    """
    # Initiate matplotlib instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    m = initiate_basemap(map_corners=map_corners, scalebar=True)

    # TO DO: remove hard coding
    # Plot fault lines, hardcoded into structure
    if show_nz_faults:
        map_extras.plot_hikurangi_trench(m)
        map_extras.plot_geonet_active_faults(m)

    # If given, plot all background stations for this given event.
    if stations is not None:
        if isinstance(stations, np.ndarray):
            plot_stations_simple(m, lats=stations[:, 0], lons=stations[:, 1])
        else:
            plot_stations(m, inv=stations, event=event,
                          annotate_names=annotate_names,
                          color_by_network=color_by_network
                          )

    # If an event is given, try to plot the focal-mechanism beachball
    if event is not None:
        event_beachball(m, event)

    # If an inventory object is given, plot source receiver line, and info
    if inv is not None:
        connect_source_receiver(m, event=event, sta=inv[0][0])
        annotate_srcrcv_information(m, event=event, inv=inv)

    # Finality
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show:
        if show == "hold":
            return f, m
        else:
            plt.show()
    plt.close()

    return f, m

