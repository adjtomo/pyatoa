#!/usr/bin/env python3
"""
Mapping functionality with Basemap
"""
import matplotlib.pyplot as plt
from numpy import ndarray
from pyatoa.visuals.map_tools import (initiate_basemap, plot_stations,
                                      plot_stations_simple, event_beachball,
                                      connect_source_receiver,
                                      annotate_srcrcv_information)


def standalone_map(map_corners, inv=None, catalog=None, annotate_names=False,
                   color_by_network=False, show=True, save=None, **kwargs):
    """
    To be used in a standalone mapmaker. Plots a catalog and an inventory
    to show all events and all stations for e.g. a given tomographic inversion

    kwargs can be passed to the matplotlib.pyplot.save() function

    :type map_corners: dict of floats
    :param map_corners: {lat_bot,lat_top,lon_left,lon_right}
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
    :type color_by_network: bool
    :param color_by_network: color station markers based on network name
    :type show: bool
    :param show: show the plot once generated, defaults to False
    :type save: str
    :param save: absolute filepath and filename if figure should be saved
    :rtype f: matplotlib figure
    :return f: figure object
    :rtype m: Basemap
    :return m: basemap object
    """
    figsize = kwargs.get("figsize", (8,10))
    dpi = kwargs.get("dpi", 100)

    # Initiate matplotlib instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    m = initiate_basemap(map_corners=map_corners, scalebar=True, **kwargs)

    # If given, plot all background stations for this given event.
    if inv:
        plot_stations(m, inv=inv, event=None, annotate_names=annotate_names,
                      color_by_network=color_by_network)

    # If an event is given, try to plot the focal-mechanism beachball
    if catalog:
        for event in catalog:
            event_beachball(m, event)

    # Finality
    f.tight_layout()
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
                annotate_names=False, color_by_network=False, show=True,
                save=None, **kwargs):
    """
    Initiate and populate a basemap object.

    Functionality to manually ignore stations based on user quality control
    Choice to annotate station and receiver which correspond to waveforms
    Calls beachball and fault tracer (optional) to populate basemap

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
    :type color_by_network: bool
    :param color_by_network: color station markers based on network name
    :type show: bool
    :param show: show the plot once generated, defaults to False
    :type save: str
    :param save: absolute filepath and filename if figure should be saved
    :rtype f: matplotlib figure
    :return f: figure object
    :rtype m: Basemap
    :return m: basemap object
    """
    figsize = kwargs.get("figsize", (6, 8))
    dpi = kwargs.get("dpi", 100)

    # Initiate matplotlib instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    m = initiate_basemap(map_corners=map_corners, scalebar=True)

    # If given, plot all background stations for this given event.
    if stations is not None:
        if isinstance(stations, ndarray):
            plot_stations_simple(m, lats=stations[:, 0], lons=stations[:, 1])
        else:
            plot_stations(m, inv=stations, event=event,
                          annotate_names=annotate_names,
                          color_by_network=color_by_network
                          )

    # Plot the station
    if inv is not None:
        plot_stations(m, inv)

    # If an event is given, try to plot the focal-mechanism beachball
    if event is not None:
        event_beachball(m, event)

    # If an inventory object is given, plot source receiver line, and info
    if (inv and event) is not None:
        connect_source_receiver(m, event=event, sta=inv[0][0])
        annotate_srcrcv_information(m, event=event, inv=inv)

    # Finality
    f.tight_layout()
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show:
        if show == "hold":
            return (f, m)
        else:
            plt.show()
    else:
        plt.close()

