"""
A class to analyze the outputs of a Seisflows inversion by
looking at misfit information en masse and producing text files
and images related to analysis of data
"""
import os
import json
import pyasdf
from glob import glob
from obspy.geodetics import gps2dist_azimuth

from pyatoa.plugins.seisflows import visuals
from pyatoa.utils.tools.srcrcv import eventid
from pyatoa.utils.asdf.extractions import count_misfit_windows


class Inspector:
    """
    This plugin object will collect information from a Pyatoa run folder and
    allow the User to easily understand statistical information or generate
    statistical plots to help understand a seismic inversion
    """
    def __init__(self, path_to_datasets=None, tag=None, misfits=True,
                 srcrcv=True, windows=True):
        """
        Inspector only requires the path to the datasets, it will then read in
        all the datasets and store the data internally. This is a long process
        but should only need to be done once.

        Allows parameters to determine what quantities are queried from dataset


        :type misfits: bool
        :param misfits: collect misfit information
        :type srcrcv: bool
        :param srcrcv: collect coordinate information
        :type cat: bool
        :param cat: collect events into a Catalog object
        :type path_to_datasets: str
        :param path_to_datasets: path to the ASDFDataSets that were outputted
            by Pyaflowa in the Seisflows workflow
        """
        # If no tag given, create dictionaries based on datasets
        self.srcrcv = {}
        self.misfits = {}
        self.windows = {}

        # If a tag is given, load rather than reading from datasets
        if tag is not None:
            self.load(tag)
        else:
            for dsfid in glob(os.path.join(path_to_datasets, "*.h5")):
                with pyasdf.ASDFDataSet(dsfid) as ds:
                    if windows:
                        self.get_windows(ds)
                    if srcrcv:
                        self.get_srcrcv(ds)
                    if misfits:
                        self.get_misfits(ds)

    @property
    def event_ids(self):
        """Return a list of all event ids"""
        return list(self.srcrcv.keys())

    @property
    def stations(self):
        """Return a list of all stations"""
        stas = []
        for event in self.srcrcv:
            for sta in self.srcrcv[event]:
                # Station names are uppercase, attributes are lowercase
                if sta.isupper():
                    stas.append(sta)
        return list(set(stas))

    @property
    def models(self):
        """Return a list of all models"""
        return list(self.sort_misfits_by_model().keys())

    @property
    def mags(self):
        """Return a dictionary of event magnitudes"""
        return self.event_info("mag")

    @property
    def times(self):
        """Return a dictionary of event origin times"""
        return self.event_info("time")

    @property
    def depths(self):
        """Return a dictionary of event depths"""
        return self.event_info("depth_m")

    def event_info(self, choice):
        """
        Return event information in a dictionary object

        :type choice: str
        :param choice: choice of key to query dictionary
        """
        info = {}
        for event in self.srcrcv.keys():
            info[event] = self.srcrcv[event][choice]
        return info

    def window_values(self, model, choice):
        """
        Return a list of all time shift values for a given model

        :type model: str
        :param model: model to query e.g. 'm00'
        :type choice: str
        :param choice: key choice for window query
        :rtype list:
        :return: list of time shift values for a given model
        """
        choices = ["cc_shift_sec", "dlna", "max_cc", "length_s", "weight"]
        assert(choice in choices), f"choice must be in {choices}"

        ret = []
        windows = self.sort_windows_by_model()
        for event in windows[model]:
            for sta in windows[model][event]:
                for cha in windows[model][event][sta]:
                    ret += windows[model][event][sta][cha][choice]
        return ret

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

    def get_srcrcv(self, ds):
        """
        Get source receiver info including coordinates, distances and BAz
        from a given dataset.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for distances
        """
        # Initialize the event as a dictionary
        eid = eventid(ds.events[0])
        self.srcrcv[eid] = {}
        self.srcrcv[eid]["lat"] = ds.events[0].preferred_origin().latitude
        self.srcrcv[eid]["lon"] = ds.events[0].preferred_origin().longitude
        self.srcrcv[eid]["depth_m"] = ds.events[0].preferred_origin().depth
        self.srcrcv[eid]["time"] = str(ds.events[0].preferred_origin().time)
        self.srcrcv[eid]["mag"] = ds.events[0].preferred_magnitude().mag

        # Loop through all the stations in the dataset
        for sta, srcrcv in ds.get_all_coordinates().items():
            self.srcrcv[eid][sta] = {}

            gcd, _, baz = gps2dist_azimuth(lat1=self.srcrcv[eid]["lat"],
                                           lon1=self.srcrcv[eid]["lon"],
                                           lat2=srcrcv["latitude"],
                                           lon2=srcrcv["longitude"]
                                           )

            # Append information to specific dictionary entry
            self.srcrcv[eid][sta]["lat"] = srcrcv["latitude"]
            self.srcrcv[eid][sta]["lon"] = srcrcv["longitude"]
            self.srcrcv[eid][sta]["elv_m"] = srcrcv["elevation_in_m"]
            self.srcrcv[eid][sta]["dist_km"] = gcd * 1E-3
            self.srcrcv[eid][sta]["baz"] = baz

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
                
    def get_windows(self, ds):
        """
        Get Window information from auxiliary_data.MisfitWindows
        
        :return: 
        """
        eid = eventid(ds.events[0])
    
        self.windows[eid] = {}
        for model in ds.auxiliary_data.MisfitWindows.list():
            self.windows[eid][model] = {}

            # For each station, determine the number of windows and total misfit
            for window in ds.auxiliary_data.MisfitWindows[model]:
                cha_id = window.parameters["channel_id"]
                net, sta, loc, cha = cha_id.split(".")
                sta_id = f"{net}.{sta}"

                dlna = window.parameters["dlnA"]
                weight = window.parameters["window_weight"]
                max_cc = window.parameters["max_cc_value"]
                length_s = (window.parameters["relative_endtime"] -
                          window.parameters["relative_starttime"]
                          )
                cc_shift_sec = window.parameters["cc_shift_in_seconds"]

                # One time initiatations of a new dictionary object
                win = self.windows[eid][model]
                if sta_id not in win:
                    win[sta_id] = {}
                if cha not in self.windows[eid][model][sta_id]:
                    win[sta_id][cha] = {"cc_shift_sec": [], "dlna": [],
                                        "weight": [], "max_cc": [],
                                        "length_s": []
                                        }

                # Append values from the parameters into dictionary object
                win[sta_id][cha]["dlna"].append(dlna)
                win[sta_id][cha]["weight"].append(weight)
                win[sta_id][cha]["max_cc"].append(max_cc)
                win[sta_id][cha]["length_s"].append(length_s)
                win[sta_id][cha]["cc_shift_sec"].append(cc_shift_sec)
                
    def save(self, tag):
        """
        Save the downloaded attributes into JSON files for re-loading

        :type tag: str
        :param tag: unique naming tag for saving json files
        """
        def write(self, suffix):  # NOQA
            """
            Convenience function to save internal attributes
            """
            obj = getattr(self, suffix)
            if obj:
                with open(f"{tag}_{suffix}.json", "w") as f:
                    print(f"writing {suffix}")
                    json.dump(obj, f, indent=4, sort_keys=True)

        for s in ["srcrcv", "misfits", "windows"]:
            write(self, s)

    def load(self, tag):
        """
        Load previously saved attributes to avoid re-processing data

        :type tag: str
        :param tag: tag to look for json files
        """
        def read(self, suffix):  # NOQA
            """
            Convenience function to read in saved files

            :type suffix: str
            :param suffix: suffix of file name
            """
            print(f"reading {suffix} file", end="... ")
            try:
                with open(f"{tag}_{suffix}.json", "r") as f:
                    setattr(self, suffix, json.load(f))
                    print("found")
            except FileNotFoundError:
                print("not found")
                pass

        for s in ["srcrcv", "misfits", "windows"]:
            read(self, s)

    def sort_misfits_by_station(self):
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

    def sort_misfits_by_model(self):
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
                if event not in misfits[model]:
                    misfits[model][event] = {}
                misfits[model][event] = self.misfits[event][model]

        return misfits

    def sort_windows_by_model(self):
        """
        Rearrage windows by model rather than event

        :rtype dict:
        :return: windows sorted by model
        """
        windows = {}
        for event in self.windows:
            for model in self.windows[event]:
                if model not in windows:
                    windows[model] = {}
                if event not in windows[model]:
                    windows[model][event] = {}
                windows[model][event] = self.windows[event][model]

        return windows

    def plot(self, choice, **kwargs):
        """
        Convenience plot function that calls to functions contained in the
        external module `pyatoa.plugins.seisflows.visuals`

        :type choice: str
        :param choice: choice of plot to create, see choices
        """
        # Set the plotting functions into an easily accesible dictionary
        choices = {
            "windows_by_distance": visuals.windows_by_distance,
            "misfit_by_distance": visuals.misfit_by_distance,
            "misfit_by_path": visuals.misfit_by_path,
            "window_by_path": visuals.window_by_path,
            "event_depths": visuals.event_depths
            }

        assert(choice in choices)
        choices[choice](self, **kwargs)


