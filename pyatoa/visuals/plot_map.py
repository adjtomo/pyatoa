#!/usr/bin/env python3
"""
misfit visualization tool to be called through adjointBuilder
produces a basemap with beachball and all available stations as well as the
relevant station highlighted. important information annotated (e.g.
misift information, distance, BAz etc.)
"""
import os
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from obspy import read_inventory
from obspy.imaging.beachball import beach

from pyatoa.utils.gathering.grab_auxiliaries import _hardcode_paths, \
    grab_geonet_moment_tensor
from pyatoa.utils.operations.source_receiver import gcd_and_baz, parse_inventory
from pyatoa.utils.operations.calculations import myround


mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def place_scalebar(m, map_corners):
    """
    put the scale bar in the corner at a reasonable distance from each edge
    :param m:
    :return:
    map_corners [lat_bot,lat_top,lon_left,lon_right]
    """
    # scale set for bottom right-hand corner
    latscale = map_corners[0] + (map_corners[1] - map_corners[0]) * 0.94  # up
    lonscale = map_corners[2] + (map_corners[3] - map_corners[2]) * 0.875 # left
    m.drawmapscale(lonscale, latscale, lonscale, latscale, 100,
                   yoffset=0.01 * (m.ymax-m.ymin))


def build_colormap(array):
    """build a custom range colormap
    """
    vmax = myround(np.nanmax(array), base=1, choice='up')
    vmin = myround(np.nanmin(array), base=1, choice='down')
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.plasma
    colormap = cm.ScalarMappable(norm=norm, cmap=cmap)

    return colormap


def plot_hikurangi_trench(m, path_):
    """trace the hikurangi trench on a basemap object 'm'
    """
    trenchcoords = np.load(path)
    lats = trenchcoords['LAT']
    lons = trenchcoords['LON']
    x, y = m(lons, lats)

    # interpolate points to make a smoother curve
    xprime = np.flip(x, axis=0)
    yprime = np.flip(y, axis=0)
    xprimenew = np.linspace(x.min(), x.max(), 100)
    yprimenew = np.interp(xprimenew, xprime, yprime)

    m.plot(xprimenew, yprimenew, '--', linewidth=1.25, color='k', zorder=2)


def plot_active_faults(m, path_):
    """plot onshore and offshore fault coordinate files
    """
    active_faults = np.load(path)
    lats = active_faults['LAT']
    lons = active_faults['LON']
    faults = active_faults['FAULT']

    for i in range(faults.min(), faults.max()+1, 1):
        indices = np.where(faults == i)
        x, y = m(lons[indices], lats[indices])
        m.plot(x, y, '--', linewidth=0.5, color='k', zorder=2, alpha=0.25)


def event_beachball(m, moment_tensor):
    """
    plot event beachball on basemap 'm' object for a given geonet moment tensor
    list, read in the from the GeoNet moment tensor csv file
    """
    eventx, eventy = m(moment_tensor['Longitude'], moment_tensor['Latitude'])
    focal_mechanism = [moment_tensor['strike2'], moment_tensor['dip2'],
                       moment_tensor['rake2']]
    b = beach(focal_mechanism, xy=(eventx, eventy), width=2.5E4, linewidth=1,
              facecolor='r')
    b.set_zorder(1000)
    ax = plt.gca()
    ax.add_collection(b)


def plot_stations(m, inv, event=None, **kwargs):
    """
    fill map with stations based on station availability and network
    """
    annotate_names = kwargs.get("annotate_names", False)
    color_by_network = kwargs.get("color_by_network", False)

    network_codes, station_codes, latitudes, longitudes = parse_inventory(
        inv, event)

    x, y = m(longitudes, latitudes)

    if color_by_network:
        color_dict = {"NZ": "black", "XX": "red", "X1": "green",
                      "X2": "blue", "YH": "yellow"}
        for net_, sta_, x_, y_ in zip(network_codes, station_codes, x, y):
            m.scatter(x_, y_, marker='v', color=color_dict[net_], edgecolor='k',
                      linestyle='-', s=60, zorder=5
                      )
    else:
        m.scatter(x, y, marker='v', color='w', edgecolor='k', linestyle='-',
                  s=60, zorder=999
                  )
    if annotate_names:
        for n_, x_, y_ in zip(station_codes, x, y):
            plt.annotate(n_, xy=(x_, y_), xytext=(x_, y_), zorder=6, fontsize=7,
                         bbox=dict(
                             facecolor='w', edgecolor='k', boxstyle='round')
                         )


def annotate_srcrcv_information(m, event, inv, config):
    """annotate event information into hard coded map area
    """
    event.origins[0].time.precision = 0
    gcdist, baz = gcd_and_baz(event, inv)

    plt.annotate(s=("{id} / {sta}\n"
                    "{date}\n"
                    "{type}={mag:.2f}\n"
                    "Depth(km)={depth:.2f}\n"
                    "Dist(km)={dist:.2f}\n"
                    "BAz(deg)={baz:.2f}").format(
        id=config.event_id, sta=inv[0][0].code,
        date=event.preferred_origin().time,
        depth=event.preferred_origin().depth*1E-3,
        type=event.preferred_magnitude().magnitude_type,
        mag=event.preferred_magnitude().mag, dist=gcdist, baz=baz
    ),
        xy=(m.xmin + (m.xmax-m.xmin) * 0.675,
            m.ymin + (m.ymax-m.ymin) * 0.035),
        multialignment='right', fontsize=10)


def connect_source_receiver(m, event, inv):
    """
    draw a dashed line connecting the station and receiver, highlight station
    :param m:
    :param event:
    :param inv:
    :return:
    """
    event_x, event_y = m(event.preferred_origin().longitude,
                         event.preferred_origin().latitude)
    station_x, station_y = m(inv[0][0].longitude, inv[0][0].latitude)
    m.plot([event_x, station_x], [event_y, station_y], '--', linewidth=1.1,
           c='k', zorder=998)
    m.scatter(station_x, station_y, marker='v', color='r', edgecolor='k',
              linestyle='-', s=75, zorder=1000)


def initiate_basemap(map_corners):
    """
    set up local map
    """
    continent_color = 'w'
    lake_color = 'w'

    # initiate map
    m = Basemap(projection='stere', resolution='h', rsphere=6371200,
                lat_0=np.mean(map_corners[:2]), lon_0=np.mean(map_corners[2:]),
                llcrnrlat=map_corners[0], llcrnrlon=map_corners[2],
                urcrnrlat=map_corners[1], urcrnrlon=map_corners[3],
                )
    m.drawcoastlines(linewidth=1.5)
    m.fillcontinents(color=continent_color, lake_color=lake_color)
    m.drawparallels(np.arange(int(map_corners[0]), int(map_corners[1]), 1),
                    labels=[1, 0, 0, 0], linewidth=0.0)
    m.drawmeridians(np.arange(int(map_corners[2]), int(map_corners[3])+1, 1),
                    labels=[0, 0, 0, 1], linewidth=0.0)
    place_scalebar(m, map_corners)

    return m


def generate_map(config, event, inv,
                 map_corners=[-42.5007,-36.9488,172.9998,179.5077],
                 show_faults=False, **kwargs):
    """
    TODO: change map corners to reflect the new mesh created in August

    initiate and populate a basemap object for New Zealands north island.
    Functionality to manually ignore stations based on user quality control
    Takes station coordinates and coloring from npz files
    Choice to annotate two stations which correspond to waveforms
    Calls beachball and trench tracer to populate basemap
    :type corners: list of floats
    :param corners: values for map corners to set bounds
     e.g. [lat_bot,lat_top,lon_left,lon_right]
     TODO: make map_corners a dictionary to keep it less ambiguous?
    """
    figsize = kwargs.get("figsize", (10, 9.4))
    dpi = kwargs.get("dpi", 100)
    show = kwargs.get("show", True)
    save = kwargs.get("save", False)

    f = plt.figure(figsize=figsize, dpi=dpi)
    m = initiate_basemap(map_corners=map_corners)

    # next section contains hardcoded paths
    background_inv = read_inventory(_hardcode_paths()['stations'])
    moment_tensor = grab_geonet_moment_tensor(config.event_id)
    if show_faults:
        warnings.warn("Plotting active faults takes some time, "
                      "please be patient", UserWarning)
        plot_hikurangi_trench(m, "./fault_coordinates/hikurangi_trench.npz")
        for path_ in ["./fault_coordinates/north_island_550_641_onshore.npz",
                      "./fault_coordinates/north_island_550_641_offshore.npz"]:
            plot_active_faults(m, path_)

    plot_stations(m, inv=background_inv, event=event, annotate_names=False,
                  color_by_network=False)
    event_beachball(m, moment_tensor)
    connect_source_receiver(m, event, inv)
    annotate_srcrcv_information(m, event, inv, config)

    if show:
        plt.show()
    if save:
        plt.savefig(save)

    return f, m


