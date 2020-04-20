#!/usr/bin/env python3
"""
Mapping functionality with Basemap
"""
import matplotlib.pyplot as plt
from numpy import ndarray
from pyatoa.utils.srcrcv import gcd_and_baz
from pyatoa.utils.asdf.extractions import count_misfit_windows, sum_misfits
from pyatoa.visuals.map_tools import (initiate_basemap, plot_stations,
                                      plot_stations_simple, event_beachball,
                                      interpolate_and_contour,
                                      connect_source_receiver,
                                      annotate_srcrcv_information)


def default_kwargs(**kwargs):
    """
    all the maps use the same kwargs so it's best to set it in one place

    :type kwargs: dict
    :param kwargs: kwargs given by the workflow or user
    :return:
    """
    figsize = kwargs.get("figsize", (6, 8))
    dpi = kwargs.get("dpi", 100)

    return figsize, dpi


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


def event_misfit_map(map_corners, ds, model, step=None, normalize=None,
                     annotate_station_info=False, contour_overlay=True,
                     filled_contours=True, show=True, save=None, **kwargs):
    """
    To be used to plot misfit information from a Pyasdf Dataset

    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing things to plot
    :type model: str
    :param model: model number, e.g. 'm00'
    :type step: str
    :param step: step number, e.g. 's00', used for title plotting only
    :type normalize: float
    :param normalize: a value to normalize all misfits by, useful if comparing
        misfits across various events
    :type map_corners: dict of floats
    :param map_corners: {lat_min,lat_max,lon_min,lon_max}
    :type annotate_station_info: bool or str
    :param annotate_station_info: annotate station names and info next to marker
        can also be 'simple'
    :type contour_overlay: bool
    :param contour_overlay: interpolate z_values and create contour
    :type filled_contours: bool
    :param filled_contours: if True, use countourf, else use contour
    :type show: bool
    :param show: show the plot once generated, defaults to False
    :type save: str
    :param save: absolute filepath and filename if figure should be saved
    :rtype f: matplotlib figure
    :return f: figure object
    :rtype m: Basemap
    :return m: basemap object
    """
    figsize, dpi = default_kwargs(**kwargs)

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
            s=(f"{event_id}\n"
               f"{event.preferred_origin().time}\n"
               f"{event.preferred_magnitude().magnitude_type}="
               f"{event.preferred_magnitude().mag:.2f}\n"
               f"Depth(km)={event.preferred_origin().depth * 1E-3:.2f}\n"),
            xy=(m.xmax - (m.xmax - m.xmin) * 0.3,
                m.ymin + (m.ymax - m.ymin) * 0.01),
            multialignment='right', fontsize=10
        )

    # Get some event-wide values
    if hasattr(ds, 'auxiliary_data'):
        if hasattr(ds.auxiliary_data, 'MisfitWindows'):
            window_dict = count_misfit_windows(ds, model,
                                               count_by_stations=True)

    # If waveforms in dataset, plot the stations
    x_values, y_values, z_values = [], [], []
    if hasattr(ds, 'waveforms'):
        for sta in ds.waveforms.list():
            # Get information to annotate next to station
            sta_anno = ""
            if 'StationXML' in ds.waveforms[sta].list():
                sta_anno += "{}\n".format(sta)
                # Get distance and backazimuth between event and station
                if event is not None:
                    gcdist, baz = gcd_and_baz(event,
                                              ds.waveforms[sta].StationXML[0][0]
                                              )
                    sta_anno += f"{gcdist:.1f}km\n{baz:.1f}deg"
                # If auxiliary data is given, add some extra information to anno
                if hasattr(ds, 'auxiliary_data'):
                    # Determine the total number of windows for the given model
                    if 'MisfitWindows' in ds.auxiliary_data.list():
                        try:
                            num_windows = window_dict[sta]
                        except KeyError:
                            num_windows = 0
                        sta_anno += f"\n{num_windows}"

                    # Determine the total misfit for a given model
                    if 'AdjointSources' in ds.auxiliary_data:
                        total_misfit = sum_misfits(ds, model, station=sta)
                        if total_misfit:
                            # If misfit should be normalized to e.g. average
                            # misfit for this event/model/step combination
                            if normalize:
                                total_misfit /= normalize
                            else:
                                sta_anno += f"/{total_misfit:.2E}"
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
        if normalize:
            cbar_label = "misfit (normalized)"
        else:
            cbar_label = "misfit"
        interpolate_and_contour(m=m, x=x_values, y=y_values, z=z_values,
                                len_xi=100, len_yi=150, colormap='viridis',
                                interpolation_method='cubic',
                                fill=filled_contours, cbar_label=cbar_label
                                )

    # Set a title
    if step:
        title = f"{model}.{step} misfit map, {len(x_values)} stations"
        if "Statistics" in ds.auxiliary_data.list():
            if model in ds.auxiliary_data.Statistics.list():
                if step in ds.auxiliary_data.Statistics[model].list():
                    pars = ds.auxiliary_data.Statistics[model][step].parameters
                    title += (
                        f"\navg misfit: {pars['average_misfit']:.2f},"
                        f"num windows: {pars['number_misfit_windows']}\n"
                        f"max misfit comp: {pars['max_misfit']:.2f}/"
                        f"{pars['max_misfit_component']}\n "
                        f"min misfit comp: {pars['min_misfit']:.2f}/"
                        f"{pars['min_misfit_component']}"
                    )
    else:
        title = f"{model} misfit map, {len(x_values)} stations"
    plt.title(title)

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
    figsize, dpi = default_kwargs(**kwargs)

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
            return f, m
        else:
            f.tight_layout()
            plt.show()
    plt.close()

    return f, m


def ray_density(inv, catalog, map_corners, show=True, save=None, **kwargs):
    """
    Take all events in a catalog and all stations in an inventory and connects
    them by straight lines to show ray coverage. This takes way too long for
    a large number of source-receiver paths so something new needs to be written

    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing relevant network and stations
    :type catalog: obspy.core.event.catalog.Catalog
    :param catalog: a catalog of events to plot
    :type map_corners: dict of floats
    :param map_corners: [lat_bot,lat_top,lon_left,lon_right]
    :type show: bool
    :param show: show the plot once generated, defaults to False
    :type save: str
    :param save: absolute filepath and filename if figure should be saved
    """
    figsize, dpi = default_kwargs(**kwargs)

    f = plt.figure(figsize=figsize, dpi=dpi)
    m = initiate_basemap(map_corners=map_corners, scalebar=True,
                         coastline_zorder=101
                         )

    for event in catalog:
        event_beachball(m, event, type="focal_mechanism", width=2.6E4)
        for net in inv:
            for sta in net:
                connect_source_receiver(m, event, sta, linestyle="-",
                                        markercolor="w",
                                        linewidth=0.5, linecolor="k",
                                        zorder=100
                                        )
    f.tight_layout()
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show:
        if show == "hold":
            return f, m
        else:
            plt.show()
    plt.close()

