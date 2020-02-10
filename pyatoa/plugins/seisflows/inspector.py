"""
A class to analyze the outputs of a Seisflows inversion by
looking at misfit information en masse and producing text files
and images related to analysis of data
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from obspy import Catalog
from pyasdf import ASDFDataSet as asdf
from obspy.geodetics import gps2dist_azimuth
from pyatoa.utils.tools.srcrcv import eventid
from pyatoa.utils.asdf.extractions import count_misfit_windows


class Inspector:
    """
    This plugin object will collect information from a Pyatoa run folder and
    allow the User to easily understand statistical information or generate
    statistical plots to help understand a seismic inversion
    """
    def __init__(self, path_to_datasets, misfits=True, coords=True,
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
        self.coords = {}
        self.misfits = {}
        self.cat = Catalog()

        for dsfid in glob(os.path.join(path_to_datasets, "*.h5")):
            with asdf(dsfid) as ds:
                if cat:
                    self.cat += ds.events
                if coords:
                    self.get_coords(ds)
                if misfits:
                    self.get_misfits(ds)

    @property
    def event_ids(self):
        return [eventid(_) for _ in self.cat]

    @property
    def stations(self):
        stas = []
        for event in self.coords:
            for sta in self.coords[event]:
                if sta not in ["lat", "lon"]:
                    stas.append(sta)
        return list(set(stas))

    @property
    def models(self):
        return list(self.sort_misfit_by_model().keys())

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
            pos = sc.get_offsets()[ind["ind"][0]]
            anno.xy = pos
            text = "{}".format("\n".join([stations[n] for n in ind["ind"]]))
            anno.set_text(text)
            anno.get_bbox_patch().set_alpha(1)
            anno.get_bbox_patch().set_facecolor("w")

        def hover(event):
            """Functionality for getting info when hovering over a point
            during an interacting mpl session
            """
            vis = anno.get_visible()
            if event.inaxes == ax:
                cont, ind = sc.contains(event)
                if cont:
                    update_anno(ind)
                    anno.set_visible(True)
                    f.canvas.draw_idle()
                # else:
                #     if vis:
                #         anno.set_visible(False)
                #         f.canvas.draw_idle()

        plt.title(f"Misfit vs. Event-Station distance\n"
                  f"model: {model} event: {event_id} station: {sta_code}")
        plt.xlabel("Distances (km)")
        plt.ylabel("Misfit")
        plt.ylim([0, np.ceil(max(misfits))])
        plt.grid(which="both", linestyle="--", alpha=0.5, linewidth=.5)

        f.canvas.mpl_connect("motion_notify_event", hover)
        plt.show()
        if save:
            plt.savefig(save)
        if show:
            f.canvas.mpl_connect("motion_notify_event", hover)
            plt.show()
