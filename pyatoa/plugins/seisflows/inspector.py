"""
A class to analyze the outputs of a Seisflows inversion by
looking at misfit information en masse and producing text files
and images related to analysis of data
"""
import os
import json
import pyasdf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from glob import glob
from obspy import Catalog, read_events
from obspy.geodetics import gps2dist_azimuth

from pyatoa.utils.tools.srcrcv import eventid, lonlat_utm
from pyatoa.utils.asdf.extractions import count_misfit_windows


class Inspector:
    """
    This plugin object will collect information from a Pyatoa run folder and
    allow the User to easily understand statistical information or generate
    statistical plots to help understand a seismic inversion
    """
    def __init__(self, path_to_datasets=None, tag=None, misfits=True, coords=True,
                 cat=True):
        """
        Inspector only requires the path to the datasets, it will then read in
        all the datasets and store the data internally. This is a long process
        but should only need to be done once.

        Allows parameters to determine what quantities are queried from dataset


        :type misfits: bool
        :param misfits: collect misfit information
        :type coords: bool
        :param coords: collect coordinate information
        :type cat: bool
        :param cat: collect events into a Catalog object
        :type path_to_datasets: str
        :param path_to_datasets: path to the ASDFDataSets that were outputted
            by Pyaflowa in the Seisflows workflow
        """
        # If a tag is given, load rather than reading from datasets
        if tag is not None:
            self.load(tag)
        else:
            # If no tag given, create dictionaries based on datasets
            self.coords = {}
            self.misfits = {}
            self.cat = Catalog()

            for dsfid in glob(os.path.join(path_to_datasets, "*.h5")):
                with pyasdf.ASDFDataSet(dsfid) as ds:
                    if cat:
                        self.cat += ds.events
                    if coords:
                        self.get_coords(ds)
                    if misfits:
                        self.get_misfits(ds)

    @property
    def event_ids(self):
        """Return a list of all event ids"""
        return [eventid(_) for _ in self.cat]

    @property
    def stations(self):
        """Return a list of all stations"""
        stas = []
        for event in self.coords:
            for sta in self.coords[event]:
                if sta not in ["lat", "lon"]:
                    stas.append(sta)
        return list(set(stas))

    @property
    def models(self):
        """Return a list of all models"""
        return list(self.sort_misfit_by_model().keys())

    def misfit_values(self, model):
        """
        Return a list of misfit values for a given model

        :type model: str
        :param model: model to query e.g. 'm00'
        :rtype list:
        :return: list of misfit values for a given model
        """
        misfit = []
        for event in self.misfits:
            for model_ in self.misfits[event]:
                if model_ != model:
                    continue
                for sta in self.misfits[event][model]:
                    misfit.append(self.misfits[event][model][sta]["msft"])
        return misfit

    def get_coords(self, ds):
        """
        Get source receiver coordinates, distances and BAz from a datset

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for distances
        """
        # Initialize the event as a dictionary
        eid = eventid(ds.events[0])
        self.coords[eid] = {}
        self.coords[eid]["lat"] = ds.events[0].preferred_origin().latitude
        self.coords[eid]["lon"] = ds.events[0].preferred_origin().longitude

        # Loop through all the stations in the dataset
        for sta, coords in ds.get_all_coordinates().items():
            self.coords[eid][sta] = {}

            gcd, _, baz = gps2dist_azimuth(lat1=self.coords[eid]["lat"],
                                           lon1=self.coords[eid]["lon"],
                                           lat2=coords["latitude"],
                                           lon2=coords["longitude"]
                                           )

            # Append information to specific dictionary entry
            self.coords[eid][sta]["lat"] = coords["latitude"]
            self.coords[eid][sta]["lon"] = coords["longitude"]
            self.coords[eid][sta]["elv_m"] = coords["elevation_in_m"]
            self.coords[eid][sta]["dist_km"] = gcd * 1E-3
            self.coords[eid][sta]["baz"] = baz

    def get_misfits(self, ds):
        """
        Get Misfit information from a dataset

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for misfit
        """
        eid = eventid(ds.events[0])

        self.misfits[eid] = {}
        for model in ds.auxiliary_data.AdjointSources.list():
            self.misfits[eid][model] = {}
            num_win = count_misfit_windows(ds, model, count_by_stations=True)

            # For each station, determine the number of windows and total misfit
            for station in ds.auxiliary_data.AdjointSources[model]:
                sta_id = station.parameters["station_id"]
                misfit = station.parameters["misfit_value"]

                # One time initiatation of a new dictionary object
                if sta_id not in self.misfits[eid][model]:
                    self.misfits[eid][model][sta_id] = {"msft": 0,
                                                        "nwin": num_win[sta_id]
                                                        }

                # Append the total number of windows, and the total misfit
                self.misfits[eid][model][sta_id]["msft"] += misfit

            # Scale the misfit of each station by the number of windows
            for sta_id in self.misfits[eid][model].keys():
                self.misfits[eid][model][sta_id]["msft"] /= \
                                    2 * self.misfits[eid][model][sta_id]["nwin"]

    def save(self, tag):
        """
        Save the downloaded attributes into JSON files for re-loading
        """
        if self.coords:
            with open(f"{tag}_coords.json", "w") as f:
                json.dump(self.coords, f, indent=4, sort_keys=True)
        if self.misfits:
            with open(f"{tag}_misfits.json", "w") as f:
                json.dump(self.misfits, f, indent=4, sort_keys=True)
        if self.cat:
            self.cat.write(f"{tag}_catalog.xml", format="quakeml")

    def load(self, tag):
        """
        Load previously saved attributes to avoid re-processing data
        """
        try:
            with open(f"{tag}_coords.json", "r") as f:
                self.coords = json.load(f)
        except FileNotFoundError:
            pass
        try:
            with open(f"{tag}_misfits.json", "r") as f:
                self.misfits = json.load(f)
        except FileNotFoundError:
            pass
        try:
            self.cat = read_events(f"{tag}_catalog.xml")
        except FileNotFoundError:
            pass

    def sort_misfit_by_station(self):
        """
        Sort the misfits collected by get_misfits() by model rather than
        by event. Returns a dictionary of misfit sorted by model

        :rtype dict:
        :return: misfits sorted by model
        """
        misfits = {}
        for event in self.misfits:
            for model in self.misfits[event]:
                if model not in misfits:
                    misfits[model] = {}
                for sta in self.misfits[event][model]:
                    if sta not in misfits[model]:
                        misfits[model][sta] = {"msft": 0, "nwin": 0,
                                               "nevents": 0}

                    # Append misfit info from each station-event in same model
                    misfits[model][sta]["msft"] += (
                        self.misfits[event][model][sta]["msft"]
                    )
                    misfits[model][sta]["nwin"] += (
                        self.misfits[event][model][sta]["nwin"]
                    )
                    misfits[model][sta]["nevents"] += 1

            # Scale the total misfit per station by the number of events
            for sta in misfits[model]:
                misfits[model][sta]["msft"] /= misfits[model][sta]["nevents"]

        return misfits

    def sort_misfit_by_model(self):
        """
        Rearrage misfits by model rather than event

        :rtype dict:
        :return: misfits sorted by model
        """
        misfits = {}
        for event in self.misfits:
            for model in self.misfits[event]:
                if model not in misfits:
                    misfits[model] = {}
                if event not in misfits:
                    misfits[model][event] = {}
                misfits[model][event] = self.misfits[event][model]

        return misfits

    def plot_misfit_by_distance(self, model, show=True, save=False,
                                event_id=None, sta_code=None):
        """
        Make a plot of misfit versus source-receiver distance

        :type model: str
        :param model: model to query, e.g. "m00"
        :type event_id: str
        :param event_id: only plot for a given event id
        :type sta_code: str
        :param sta_code: only plot for a given station
        :type show: bool
        :param show: Show the plot
        :type save: str
        :param save: fid to save the given figure
        :type event_id: str
        :param event_id: only plot for a given event
        """
        assert self.coords, "No distance information"
        assert(model in self.models), f"Model must be in {self.models}"
        if event_id:
            assert(event_id in self.event_ids), \
                f"event_id must be in {self.event_ids}"
        if sta_code:
            assert(sta_code in self.stations), \
                f"sta_code must be in {self.stations}"

        # Collect misfit and distance information per event
        f, ax = plt.subplots()
        stations, distances, misfits = [], [], []
        for i, event in enumerate(self.misfits):
            if event_id and event != event_id:
                continue
            for sta in self.misfits[event][model]:
                if sta_code and sta != sta_code:
                    continue
                misfit = self.misfits[event][model][sta]["msft"]
                stations.append(f"{sta}\n{event}\n{misfit:.2f}")
                distances.append(self.coords[event][sta]["dist_km"])
                misfits.append(misfit)

        # Quick plot so each event gets a different color
        sc = plt.scatter(distances, misfits, s=7, label=event, marker="o",
                         edgecolors="k", color="w")

        plt.title(f"Misfit vs. Event-Station distance\n"
                  f"model: {model} / event: {event_id} / station: {sta_code}")
        plt.xlabel("Distances (km)")
        plt.ylabel("Misfit")
        plt.ylim([0, np.ceil(max(misfits))])
        plt.grid(which="both", linestyle="--", alpha=0.5, linewidth=.5)

        hover = hover_on_plot(f, ax, sc, stations)

        if save:
            plt.savefig(save)
        if show:
            f.canvas.mpl_connect("motion_notify_event", hover)
            plt.show()

        return f, ax

    def plot_misfit_by_path(self, model, event_id=None, sta_code=None,
                            threshold=None, hover_on_lines=False,
                            colormap=plt.cm.viridis, show=True, save=None):
        """
        Plot misfit by source-receiver path to try to highlight portions of
        the model that may give rise to larger misfit

        :type model: str
        :param model: model to choose for misfit
        :type event_id: str
        :param event_id: only plot for a given event id
        :type sta_code: str
        :param sta_code: only plot for a given station
        :type threshold: float
        :param threshold: normalized misfit value below which, paths will not be
            plotted. Good for looking at only high misfit values. Values must
            be between 0 and 1
        :type hover_on_lines: bool
        :param hover_on_lines: for interactive plots, show misfit values when
            hovering over the source-receiver raypath lines. This can get a bit
            messy with a lot of lines
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        misfits = self.sort_misfit_by_model()
        assert(model in misfits)
        assert(model in self.models), f"Model must be in {self.models}"
        if event_id:
            assert(event_id in self.event_ids), \
                f"event_id must be in {self.event_ids}"
        if sta_code:
            assert(sta_code in self.stations), \
                f"sta_code must be in {self.stations}"
        if threshold:
            assert(0 <= threshold <= 1), "Threshold must be between 0 and 1"

        f, ax = plt.subplots()
        cmap, norm, cbar = colormap_colorbar(colormap, 0, 1.1, .1)

        # Empty lists to be filled when looping through data
        stations_plotted = []
        event_x, event_y, event_s = [], [], []
        station_x, station_y, station_s = [], [], []

        for event in misfits[model]:
            if event_id and event != event_id:
                continue
            ev_x, ev_y = lonlat_utm(lon_or_x=self.coords[event]["lon"],
                                    lat_or_y=self.coords[event]["lat"],
                                    utm_zone=-60, inverse=False
                                    )
            event_s.append(event)
            event_x.append(ev_x)
            event_y.append(ev_y)

            # Plot each station and a connecting line
            for sta in misfits[model][event]:
                if sta_code and sta != sta_code:
                    continue
                # Convert station coordinates and append to list
                sta_x, sta_y = lonlat_utm(
                    lon_or_x=self.coords[event][sta]["lon"],
                    lat_or_y=self.coords[event][sta]["lat"],
                    utm_zone=-60, inverse=False
                )
                # Ensure we only plot each station once
                if (sta_x, sta_y) not in stations_plotted:
                    stations_plotted.append((sta_x, sta_y))
                    station_x.append(sta_x)
                    station_y.append(sta_y)
                    station_s.append(sta)

                # Normalize misfit by the largest value, plot srcrcv as line
                misfit = (misfits[model][event][sta]["msft"] /
                          max(self.misfit_values(model))
                          )

                # Ignore misfit values below a certain threshold
                if threshold and misfit < threshold:
                    continue
                line, = plt.plot([ev_x, sta_x], [ev_y, sta_y],
                                 c=cmap(norm(misfit)), alpha=misfit,
                                 zorder=10 + misfit)
                if hover_on_lines:
                    hover_on_plot(f, ax, line, [f"{misfit:.2f}"],
                                  dissapear=True)

        # Plot sources and receivers as scatterplots
        sc_events = plt.scatter(event_x, event_y, marker="o", c="w",
                                edgecolors="r", s=10, zorder=100)
        sc_stations = plt.scatter(station_x, station_y, marker="v",
                                  edgecolors="g", c="w", s=10, zorder=100)

        # Make the plot a bit prettier
        plt.grid(which="both", linestyle=":")
        plt.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
        plt.xlabel("Easting (m)")
        plt.ylabel("Northing (m)")
        plt.title("Source-Receiver misfit\n"
                  f"model: {model} / event: {event_id} / station: {sta_code}")

        hover_on_plot(f, ax, sc_events, event_s)
        hover_on_plot(f, ax, sc_stations, station_s)
        plt.show()
        # Make source and receiver markers interactive
        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, sc_events, event_s)
            hover_on_plot(f, ax, sc_stations, station_s)
            plt.show()

        return f, ax


def colormap_colorbar(cmap, vmin, vmax, dv):
    """
    Create a custom colormap and colorbar

    :type cmap: matplotlib.colors.ListedColormap
    :param cmap: colormap to use, called like plt.cm.viridis
    :type vmin: float
    :param vmin: min value for colormap
    :type vmax: float
    :param vmax: max value for colormap
    :type dv: float
    :param dv: colormap boundary separations
    :rtype:
    :return:
    """
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    sm.set_clim(vmin, vmax)
    cbar = plt.colorbar(sm, boundaries=np.arange(vmin, vmax, dv), shrink=0.9)

    return cmap, norm, cbar


def hover_on_plot(f, ax, obj, values, dissapear=False):
    """
    Allow for hover on a plot for custom annotated information

    From Stackoverflow:
        https://stackoverflow.com/questions/7908636/possible-to-make-labels-
        appear-when-hovering-over-a-point-in-matplotlib

    :type f: matplotlib.figure.Figure
    :param f: figure object for hover
    :type ax: matplotlib.axes._subplot.AxesSubplot
    :param ax: axis object for hover
    :type obj: matplotlib.collections.PathCollection or
                matplotlib.lines.Line2D
    :param obj: scatter plot, returned from plt.scatter() or plt.plot()
    :type values: list of str
    :param values: list of annotations
    :type dissapear: bool
    :param dissapear: annotations dissapear when mouse moves off
    :rtype hover: function
    :return hover: the hover function to be passed to matplotlib
    """
    # Make some objects to be used for hover-over capabilities
    anno = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                       textcoords="offset points",
                       bbox=dict(boxstyle="round", fc="w"),
                       arrowprops=dict(arrowstyle="->")
                       )
    anno.set_visible(False)

    def update_anno(ind):
        """Functionality for getting info when hovering over a point
        during an interacting mpl session
        """
        # Choice between a 2D line and a scatter plot
        if isinstance(obj, mpl.lines.Line2D):
            x, y = obj.get_data()
            anno.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
        elif isinstance(obj, mpl.collections.PathCollection):
            pos = obj.get_offsets()[ind["ind"][0]]
            anno.xy = pos

        text = "{}".format("\n".join([values[n] for n in ind["ind"]]))
        anno.set_text(text)
        anno.get_bbox_patch().set_facecolor("w")
        # anno.get_bbox_patch().set_alpha(0)

    def hover(event):
        """Functionality for getting info when hovering over a point
        during an interacting mpl session
        """
        vis = anno.get_visible()
        if event.inaxes == ax:
            cont, ind = obj.contains(event)
            if cont:
                update_anno(ind)
                anno.set_visible(True)
                f.canvas.draw_idle()
            # This code snippet will make the annotation dissapear when
            # the mouse moves away
            else:
                if vis and dissapear:
                    anno.set_visible(False)
                    f.canvas.draw_idle()

    f.canvas.mpl_connect("motion_notify_event", hover)
    return hover

