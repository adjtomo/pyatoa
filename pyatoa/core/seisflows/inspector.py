#!/usr/bin/env python3
"""
A class to analyze the outputs of a Seisflows inversion by
looking at misfit information en masse and producing dictionary objects that
can quickly be queried by built-in functions.to look at stats and figures to
understand the progress of an inversion

To Do:
    -Incorporate Pandas into this class to make it easier to work with
     such large scale data. Dictionaries and lists will be cumbersome for large
     scale inversions

The Inspector carries around information about:

srcrcv: Sources and receivers, locations, depths, distances, backazimuths
    origintimes for sources
windows: All information regarding individual misfit windows, e.g. time shift,
    amplitude anomaly, window start and end, correlations
misfit: En masse misfit information from Pyadjoint, including number of windows
    and total misfit
"""
import os
import json
import pyasdf
import traceback
import numpy as np
import pandas as pd
from glob import glob
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth

from pyatoa.utils.form import event_name
from pyatoa.utils.calculate import abs_max
from pyatoa.utils.srcrcv import lonlat_utm
from pyatoa.utils.asdf.extractions import count_misfit_windows
from pyatoa.visuals.gadget import Gadget


class Inspector(Gadget):
    """
    This plugin object will collect information from a Pyatoa run folder and
    allow the User to easily understand statistical information or generate
    statistical plots to help understand a seismic inversion
    
    Inherits plotting capabilities from the Artist class to reduce clutter.
    """
    def __init__(self, tag=None, path=None, misfits=True, srcrcv=True,
                 windows=True, utm=-60):
        """
        Inspector only requires the path to the datasets, it will then read in
        all the datasets and store the data internally. This is a long process
        but should only need to be done once.

        Allows parameters to determine what quantities are queried from dataset
        Inherits plotting functionality from the Visuals class

        :type misfits: bool
        :param misfits: collect misfit information
        :type srcrcv: bool
        :param srcrcv: collect coordinate information
        :type path: str
        :param path: path to the ASDFDataSets that were outputted
            by Pyaflowa in the Seisflows workflow
        """
        # If no tag given, create dictionaries based on datasets
        self.srcrcv = {}
        self.misfits = {}
        self.windows = {}
        self._utm = utm

        # Internal attributes for properties, to be filled by getters
        self._stations = None
        self._events = None
        self._str = None
        self._models = None
        self._steps = None
        self._iterations = 0

        # If a tag is given, load rather than reading from datasets
        if tag is not None:
            try:
                self.read(tag)
            except FileNotFoundError:
                print("file not found")
        elif path is not None:
            dsfids = glob(os.path.join(path, "*.h5"))
            for i, dsfid in enumerate(dsfids):
                print(f"{os.path.basename(dsfid):<25} {i:0>2}/{len(dsfids)}", 
                      end="...")
                try:
                    self.append(dsfid, windows, srcrcv, misfits)
                    print("done")
                except KeyError as e:
                    print(f"error: {e}")
                    traceback.print_exc()
                    continue
        self._get_info()
        self._get_str()
    
    def _get_str(self):
        """
        Get the string representation once and save as internal attribute
        """
        # Get a list of internal public methods
        method_list = [func for func in dir(self)
                       if callable(getattr(self, func))
                       and not func.startswith("_")]
        # Get a list of internal public variables
        variable_list = [var for var in vars(self).keys()
                         if not var.startswith("_")]
        str_out = "INSPECTOR\n"
        str_out += ("Attributes:\n"
                    f"\tevents: {len(self._events)}\n"
                    f"\tstations: {len(self._stations)}\n"
                    f"\tmodels: {len(self._models)}\n")
        str_out += "Variables:\n"
        for var in variable_list:
            str_out += f"\t{var}\n"
        str_out += "Methods:\n"
        for meth in method_list:
            if meth in self.__class__.__dict__.keys():
                str_out += f"\t{meth}\n"
        str_out += "Plotting (Gadget):\n"
        for meth in method_list:
            if meth not in self.__class__.__dict__.keys():
                str_out += f"\t{meth}\n"

        self._str = str_out

    def __str__(self, flush=False):
        """
        Return a list of all variables and functions available for quick ref
        """
        return self._str

    @property
    def events(self):
        """Return a list of all event ids"""
        return self._events

    @property
    def stations(self):
        """Return a list of all stations"""
        return self._stations

    @property
    def models(self):
        """Return a list of all models"""
        return self._models

    @property
    def steps(self):
        """Returns a dictionary of models with values listing steps for each"""
        return self._steps

    @property
    def iterations(self):
        """Returns the number of iterations, or the sum of all step counts"""
        return self._iterations

    @property
    def mags(self):
        """Return a dictionary of event magnitudes"""
        return self._event_info("mag")

    @property
    def times(self):
        """Return a dictionary of event origin times"""
        return self._event_info("time")

    @property
    def depths(self):
        """Return a dictionary of event depths in units of meters"""
        return self._event_info("depth_m")

    def _get_info(self):
        """
        One-time retrieve lists of station names and event ids, based on the 
        fact that stations are separated by a '.' and events are not.
        Also grabs model and step information based on organization of misfits
        """
        events, stations = [], []
        for key in self.srcrcv.keys():
            if "." in key:
                stations.append(key)
            else:
                events.append(key)
        self._stations = stations
        self._events = events

        # Only need to check one event because they share same model/step info
        _event = list(self.misfits.keys())[0]
        models = list(self.misfits[_event].keys())
        steps = {m: [s for s in self.misfits[_event][m].keys()] for m in models}
        iterations = 0
        for m in models:
            iterations += len(steps[m])

        self._models = models
        self._steps = steps
        self._iterations = iterations

    def _event_info(self, choice):
        """
        Return chosen event information in a dictionary object, called by
        class properties.

        :type choice: str
        :param choice: choice of key to query dictionary
            choices are: depth_m, lat, lon, mag, time, utm_x, utm_y
        """
        choices = ["depth_m", "lat", "lon", "mag", "time", "utm_x", "utm_y"]
        assert (choice in choices), f"Choice must be in {choices}"
        info = {}
        for event in self.srcrcv.keys():
            try:
                info[event] = self.srcrcv[event][choice]
            except KeyError:
                continue
        return info

    def to_dataframe(self, choice):
        """
        Convert internal nested dictionaries into a multiindex Pandas dataframe.
        Currently not used but may be useful in the future if nested
        dictionaries become too cumbersome.

        :param choice:
        :return:
        """
        mydict = self.sort_by_model(choice=choice)

        return pd.DataFrame.from_dict(
            {(model, step, event, station): mydict[model][step][event][station]
             for model in mydict.keys() for step in mydict[model].keys()
             for event in mydict[model][step].keys()
             for station in mydict[model][step][event].keys()},
            orient='index')

    def event_stats(self, event, model, step):
        """
        Get misfit, number of stations, number of windows and total misfit
        on an event wide basis
        :return:
        """
        assert (model is not None and step is not None), \
            "model and step must be specified for event-wise printing"
        if model:
            assert (model in self.models), "model not available"
        if step:
            assert (step in self.steps[model]), "step not available"

        for e in self.misfits:
            event_msft, nwin = 0, 0
            event = self.misfits[e][model][step]
            nsta = len(event.keys())
            for sta in event:
                event_msft += event[sta]["msft"]
                nwin += event[sta]["nwin"]
            event_msft /= (2 * nwin)
            print(event_str.format(e=e, st=nsta, w=nwin,
                                   m=event_msft))

    def pprint(self, choice, model=None, step=None, event=None, station=None):
        """
        Pretty print misfit information for assessment of inversion behavior

        By model header:
            model: model number e.g. 'm00'
            step: step count e.g. 's00'
            nevt: number of events for a given model and step
            nsta: number of stations for a given model and step
            nwndw: number of misfit windows for a given model and step
            dt_mu: mean time shift in units of seconds
            dt_std: 1-sigma for time shifts in units of seconds
            da_mu: mean amplitude anomaly dlnA
            da_std: 1-sigma for amplitude anomaly dlna
            sum(dt): sum of absolute values of time shifts in units of seconds
            misfit: total misfit scaled by number of windows and events

        :param model:
        :param event:
        :param station:
        :return:
        """
        # Print statements for each model, which states the number of events,
        # stations and windows contained for each model
        if choice == "model":
            misfits = self.sort_by_model(choice="misfits")
            windows = self.sort_by_model(choice="windows")
            # Preformatted strings
            model_str = (
                "{mdl:>5}{stp:>5}{evt:>5}{sta:>5}{wdw:>6}{dt_mu:6.2f}"
                "{dt_std:7.2f}{da_mu:7.2f}{da_std:7.2f}{abs_dt:8.1f}  "
                "{mft:6.2E}"
            )
            header = (
                "model step nevt nsta nwndw dt_mu dt_std  da_mu da_std  sum(dt)" 
                "   misfit"
            )
            print(header)
            for m in misfits:
                if model and m != model:
                    continue
                print("=" * len(header))
                for i, s in enumerate(misfits[m]):
                    if step and s != step:
                        continue
                    nevents = len(misfits[m][s].keys())
                    # All information should be collected per step count
                    stations, tshifts, dlnas = [], [], []
                    total_windows, total_misfit = 0, 0
                    for e in misfits[m][s]:
                        event_misfit, event_nwin = 0, 0
                        if event and e != event:
                            continue
                        stations += list(misfits[m][s][e].keys())
                        for sta in misfits[m][s][e]:
                            if station and sta != station:
                                continue
                            event_nwin += misfits[m][s][e][sta]["nwin"]
                            event_misfit += misfits[m][s][e][sta]["msft"]
                            # Get list of time shift and dlna vals from windows
                            for c in windows[m][s][e][sta]:
                                window = windows[m][s][e][sta][c]
                                tshifts += window["cc_shift_sec"]
                                dlnas += window["dlna"]
                        # Event misfit: Tape et al. (2010) Eq. 6
                        if event_nwin:
                            total_misfit += event_misfit / (2 * event_nwin)
                        total_windows += event_nwin
                        assert (total_windows == len(tshifts)), \
                            "number of windows and measurements mismatch"

                    # Total misfit: Tape et al. (2010) Eq. 7
                    total_misfit /= nevents
                    # Only print the model once in the table
                    if i == 0:
                        mdl = m
                    else:
                        mdl = ""
                    print(model_str.format(
                        mdl=mdl, stp=s, evt=nevents, sta=len(set(stations)),
                        wdw=total_windows, dt_mu=np.mean(tshifts),
                        dt_std=np.std(tshifts), da_mu=np.mean(dlnas),
                        da_std=np.std(dlnas), abs_dt=sum(np.abs(tshifts)),
                        mft=total_misfit,)
                    )
            print("\n")

        elif choice == "event":
            assert(model is not None and step is not None), \
                "model and step must be specified for event-wise printing"
            if model:
                assert(model in self.models), "model not available"
            if step:
                assert(step in self.steps[model]), "step not available"
            header = "   event_id nsta nwin    misfit"

            event_str = "{e:>11}{st:>5}{w:>5}{m:10.2E}"
            print(header + "\n" + "=" * len(header))
            for e in self.misfits:
                event_msft, nwin = 0, 0
                event = self.misfits[e][model][step]
                nsta = len(event.keys())
                for sta in event:
                    event_msft += event[sta]["msft"]
                    nwin += event[sta]["nwin"]
                event_msft /= (2 * nwin)
                print(event_str.format(e=e, st=nsta, w=nwin,
                                       m=event_msft))

        elif choice == "station":
            if model:
                assert(model in self.models), "model not available"
            if step:
                assert(step in self.steps[model]), "step not available"
            header = " station nevents nwin    misfit"
            station_str = "{s:>8}{e:>8}{w:>5}{m:10.2E}"

            misfits = self.sort_misfits_by_station()
            for m in misfits:
                if model and m != model:
                    continue
                for s in misfits[m]:
                    if step and s != step:
                        continue
                    print(f"\nModel: {m} / Step: {s}")
                    print(header + "\n" + "=" * len(header))
                    for st in misfits[m][s]:
                        sta = misfits[m][s][st]
                        print(
                            station_str.format(s=st, e=sta["nevents"],
                                               w=sta["nwin"],
                                               m=sta["msft"] / (2 * sta["nwin"])
                                               )
                              )
        else:
            print("'choice' must be model, event, or station")

    def append(self, dsfid, windows=True, srcrcv=True, misfits=True):
        """
        Append a new pyasdf.ASDFDataSet file to the current set of internal
        statistics.

        :type dsfid: str
        :param dsfid: fid of the dataset
        :type windows: bool
        :param windows: get window info
        :type srcrcv: bool
        :param srcrcv: get srcrcv info
        :type misfits: bool
        :param misfits: get misfit info
        """
        try:
            with pyasdf.ASDFDataSet(dsfid) as ds:
                if windows:
                    self._get_windows_from_dataset(ds)
                if srcrcv:
                    self._get_srcrcv_from_dataset(ds)
                if misfits:
                    self._get_misfits_from_dataset(ds)
                return
        except OSError:
            print(f"error: already open")
            return

    def save(self, tag, path="./"):
        """
        Save the downloaded attributes into JSON files for easier re-loading.

        :type tag: str
        :param tag: unique naming tag for saving json files
        :type path: str
        :param path: optional path to save to, defaults to cwd
        """
        variables = [_ for _ in vars(self).keys() if '_' not in _]
        # Save all components into a single dictionary
        save_dict = {}
        for v in variables:
            if hasattr(self, v):
                save_dict[v] = getattr(self, v)

        # Save all outputs
        with open(os.path.join(path, f"{tag}.json"), "w") as f:
            print("writing file")
            json.dump(save_dict, f, indent=4, sort_keys=True)

    def write(self, tag, path="./"):
        """
        Same as save(), but I kept writing .write() so I figured i'd have it

        :type tag: str
        :param tag: unique naming tag for saving json files
        :type path: str
        :param path: optional path to file, defaults to cwd
        """
        self.save(tag, path)

    def read(self, tag, path="./"):
        """
        Load previously saved attributes to avoid re-processing data.

        :type tag: str
        :param tag: tag to look for json files
        :type path: str
        :param path: optional path to file, defaults to cwd
        """
        variables = [_ for _ in vars(self).keys() if '_' not in _]
        if tag.endswith(".json"):
            tag = tag.split(".")[0]

        print(f"reading file", end="... ")
        with open(os.path.join(path, f"{tag}.json"), "r") as f:
            loaded_variables = json.load(f)
            print("found")
            for v in variables:
                setattr(self, v, loaded_variables[v])

    def _get_srcrcv_from_dataset(self, ds):
        """
        Get source receiver info including coordinates, distances and BAz
        from a given dataset. Appends to internal variable srcrcv.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for distances
        """
        # Initialize the event as a dictionary
        eid = event_name(ds=ds)

        # Get UTM projection of event coordinates
        ev_x, ev_y = lonlat_utm(
            lon_or_x=ds.events[0].preferred_origin().longitude,
            lat_or_y=ds.events[0].preferred_origin().latitude,
            utm_zone=self._utm, inverse=False
        )

        self.srcrcv[eid] = {"lat": ds.events[0].preferred_origin().latitude,
                            "lon": ds.events[0].preferred_origin().longitude,
                            "depth_m": ds.events[0].preferred_origin().depth,
                            "time": str(ds.events[0].preferred_origin().time),
                            "mag": ds.events[0].preferred_magnitude().mag,
                            "utm_x": ev_x,
                            "utm_y": ev_y
                            }

        # Loop through all the stations in the dataset
        for sta, sta_info in ds.get_all_coordinates().items():
            # Append station location information one-time to dictionary
            if sta not in self.srcrcv:
                sta_x, sta_y = lonlat_utm(lon_or_x=sta_info["longitude"],
                                          lat_or_y=sta_info["latitude"],
                                          utm_zone=self._utm, inverse=False
                                          )
                self.srcrcv[sta] = {"lat": sta_info["latitude"],
                                    "lon": sta_info["longitude"],
                                    "elv_m": sta_info["elevation_in_m"],
                                    "utm_x": sta_x,
                                    "utm_y": sta_y
                                    }

            # Append src-rcv distance and backazimuth to specific event
            gcd, _, baz = gps2dist_azimuth(lat1=self.srcrcv[eid]["lat"],
                                           lon1=self.srcrcv[eid]["lon"],
                                           lat2=self.srcrcv[sta]["lat"],
                                           lon2=self.srcrcv[sta]["lon"]
                                           )
            self.srcrcv[eid][sta] = {"dist_km": gcd * 1E-3, "baz": baz}

    def _get_misfits_from_dataset(self, ds):
        """
        Get Misfit information from a dataset.
        Appends to internal variable misfit.

        Note:
            Misfit is collected unscaled, and will need to be scaled by
            the number of windows and number of events to get values matching
            those of the inversion.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for misfit
        """
        eid = event_name(ds=ds)

        self.misfits[eid] = {}
        for model in ds.auxiliary_data.AdjointSources.list():
            self.misfits[eid][model] = {}
            for step in ds.auxiliary_data.AdjointSources[model].list():
                self.misfits[eid][model][step] = {}
                num_win = count_misfit_windows(ds, model, step,
                                               count_by_stations=True)

                # For each station, determine the number of windows and misfit
                for station in \
                        ds.auxiliary_data.AdjointSources[model][step]:
                    sta_id = station.parameters["station_id"]
                    misfit = station.parameters["misfit_value"]

                    # One time initiatation of a new dictionary object which
                    # contains the number of windows for a given station
                    if sta_id not in self.misfits[eid][model][step]:
                        self.misfits[eid][model][step][sta_id] = {
                            "msft": 0, "nwin": num_win[sta_id]
                        }

                    # Append the total number of windows, and the total misfit
                    self.misfits[eid][model][step][sta_id]["msft"] += misfit

    def _get_windows_from_dataset(self, ds):
        """
        Get Window information from auxiliary_data.MisfitWindows
        Appends to internal variable windows.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for misfit
        """
        eid = event_name(ds=ds)

        self.windows[eid] = {}
        for model in ds.auxiliary_data.MisfitWindows.list():
            self.windows[eid][model] = {}
            for step in ds.auxiliary_data.MisfitWindows[model].list():
                self.windows[eid][model][step] = {}
                # For each station, determine number of windows and total misfit
                for window in \
                        ds.auxiliary_data.MisfitWindows[model][step]:
                    cha_id = window.parameters["channel_id"]
                    net, sta, loc, cha = cha_id.split(".")
                    sta_id = f"{net}.{sta}"

                    dlna = window.parameters["dlnA"]
                    weight = window.parameters["window_weight"]
                    max_cc = window.parameters["max_cc_value"]
                    length_s = (window.parameters["relative_endtime"] -
                                window.parameters["relative_starttime"]
                                )
                    rel_start = window.parameters["relative_starttime"]
                    rel_end = window.parameters["relative_endtime"]
                    cc_shift_sec = window.parameters["cc_shift_in_seconds"]

                    # One time initiatations of a new dictionary object
                    win = self.windows[eid][model][step]
                    if sta_id not in win:
                        win[sta_id] = {}
                    if cha not in self.windows[eid][model][step][sta_id]:
                        win[sta_id][cha] = {"cc_shift_sec": [], "dlna": [],
                                            "weight": [], "max_cc": [],
                                            "length_s": [], "rel_start": [],
                                            "rel_end": []
                                            }

                    # Append values from the parameters into dictionary object
                    win[sta_id][cha]["dlna"].append(dlna)
                    win[sta_id][cha]["weight"].append(weight)
                    win[sta_id][cha]["max_cc"].append(max_cc)
                    win[sta_id][cha]["length_s"].append(length_s)
                    win[sta_id][cha]["rel_end"].append(rel_end)
                    win[sta_id][cha]["rel_start"].append(rel_start)
                    win[sta_id][cha]["cc_shift_sec"].append(cc_shift_sec)

    def window_values(self, model, step, choice="cc_shift_sec",
                      values_only=False):
        """
        Sorts through misfit windows and returns a sorted list of `choice`.
        Useful for looking for maximum time shift or misfit and the
        corresponding event, component, etc.

        If values_only == True:
            returns a single list of values for the given choice
        elif values_only == False:
            returns a list of lists, where each list contains corresponding
            event, station, component and window information in the following:
            [event_id, station, component, window_number, value]

        Choices are: "cc_shift_sec", "dlna", "max_cc", "length_s", "weight"
                     "rel_end", "rel_start"

        :type model: str
        :param model: model to query, e.g. 'm00'
        :type step: str
        :param step: step count to query, e.g. 's00'
        :type choice: str
        :param choice: choice of measurement to return
        :type values_only: bool
        :param values_only: only return the values, not the corresponding info
        :rtype a_out: list
        :return a_out: values of windows for given choice
        """
        choices = ["cc_shift_sec", "dlna", "max_cc", "length_s", "weight",
                   "rel_start", "rel_end"]
        assert(choice in choices), f"choice must be in {choices}"

        windows = self.sort_by_model("windows")[model][step]
        a_out = []
        for event in windows:
            for sta in windows[event]:
                for comp in windows[event][sta]:
                    # Each window gets its own line in the array
                    for w, val in enumerate(windows[event][sta][comp][choice]):
                        if values_only:
                            a_out.append(float(val))
                        else:
                            a_out.append([event, sta, comp, int(w), float(val)])

        # Sort by the values
        if values_only:
            a_out.sort(reverse=True)
        else:
            a_out.sort(key=lambda x: x[-1], reverse=True)

        return a_out

    def misfit_values(self, model, step):
        """
        Return a single list of misfit values for a given model and step.

        :type model: str
        :param model: model to query e.g. 'm00'
        :type step: str
        :param step: step count to query, e.g. 's00'
        :rtype list:
        :return: list of misfit values for a given model
        """
        misfit = []
        for event in self.misfits:
            for model_ in self.misfits[event]:
                if model_ != model:
                    continue
                for step_ in self.misfits[event][model][step]:
                    if step_ != step:
                        continue
                    for sta in self.misfits[event][model][step]:
                        misfit.append(
                            self.misfits[event][model][step][sta]["msft"])
        return misfit

    def sort_by_model(self, choice):
        """
        Rearrage chosen dictionary so that the first layer corresponds to model
        rather than the default sorting of by event.

        Nesting layers goes from:
            event > model > step
        to
            model > step > event

        :rtype dict:
        :return: misfits sorted by model
        """
        assert(choice in ["misfits", "windows"]), \
            "choice must be in 'misfits', 'windows'"
        if choice == "misfits":
            d = self.misfits
        elif choice == "windows":
            d = self.windows

        # Nasty one liner to rearrage a double-nested dictionary from
        # event/model/step - > model/step/event
        return {m: {s: {e: d[e][m][s] for e in d if m in d[e]}
                    for s in self.steps[m]}
                for m in self.models
                }

    def sort_misfits_by_station(self):
        """
        Rearrage misfits so that the first dictionary layer is sorted by station
        rather than the default sorting of by event. Misfit is returned scaled
        by the number of windows collected for the given station.

        :rtype dict:
        :return: misfits sorted by station
        """
        misfits = {}
        for event in self.misfits:
            for model in self.misfits[event]:
                if model not in misfits:
                    misfits[model] = {}
                for step in self.misfits[event][model]:
                    if step not in misfits[model]:
                        misfits[model][step] = {}
                    for sta in self.misfits[event][model][step]:
                        if sta not in misfits[model][step]:
                            misfits[model][step][sta] = {"msft": 0, "nwin": 0,
                                                         "nevents": 0}

                        # Append misfit info from each station-event 
                        misfits[model][step][sta]["msft"] += (
                            self.misfits[event][model][step][sta]["msft"]
                        )
                        misfits[model][step][sta]["nwin"] += (
                            self.misfits[event][model][step][sta]["nwin"]
                        )
                        misfits[model][step][sta]["nevents"] += 1

                # Scale the total misfit per station by number of windows
                for sta in misfits[model][step]:
                    misfits[model][step][sta]["msft"] /= \
                        2 * misfits[model][step][sta]["nwin"]

        return misfits

    def measurement_length(self, choice):
        """
        Find the cumulative length of misfit windows for a given model/step,
        or the number of misfit windows for a given model/step.

        :rtype: dict
        :return: cumulative window length in seconds for each model/step
        """
        assert(choice in ["cum_win_len", "num_windows"]), \
            "choice must be 'cum_win_len or 'num_windows"
        self._get_str()
        windows = self.sort_by_model("windows")

        cumulative_window_length = {m: {s: 0 for s in windows[m]}
                                    for m in windows
                                    }
        number_of_windows = {m: {s: 0 for s in windows[m]}
                             for m in windows
                             }
        for model in windows:
            for step in windows[model]:
                for event in windows[model][step]:
                    for sta in windows[model][step][event]:
                        for cha in windows[model][step][event][sta]:
                            w_ = windows[model][step][event][sta][cha]
                            # Each of these entries is a list of the same len
                            number_of_windows[model][step] += len(w_["dlna"])
                            for length in w_["length_s"]:
                                cumulative_window_length[model][step] += length
        if choice == "cum_win_len":
            return cumulative_window_length
        elif choice == "num_windows":
            return number_of_windows

    def sum_misfits(self):
        """
        Sum the total misfit for a given model based on the individual
        misfits for each misfit window, and the number of sources used.

        :rtype: dict
        :return: total misfit for each model in the class
        """
        misfits = self.sort_by_model("misfits")

        # One liner for nested dictionaries for each step to have a misfit val
        cmsft = {key: {step: 0 for step in misfits[key]} for key in misfits}
        for model in misfits:
            for step in misfits[model]:
                total_misfit = 0
                for e, event in enumerate(misfits[model][step]):
                    ev_msft, ev_nwin = 0, 0
                    for sta in misfits[model][step][event]:
                        # need to scale misfit by nwin (Tape 2010 eq. 6)
                        ev_msft += misfits[model][step][event][sta]["msft"]
                        ev_nwin += misfits[model][step][event][sta]["nwin"]
                    total_misfit += ev_msft / (2 * ev_nwin)
                # Divide by the number of sources (Tape 2010 eq. 7)
                cmsft[model][step] = total_misfit / (e + 1)

        return cmsft

    def exclude_events(self, coords=None, depth_min=None, depth_max=None,
                       mag_min=None, mag_max=None, starttime=None,
                       endtime=None):
        """
        Go through misfits and windows and remove events that fall outside
        a certain bounding box. This is useful for looking at certain regions
        of the map, or certain depths, magnitudes, origintimes.

        NOTE:
            Makes edits to arrays in place so data should be saved beforehand.

        :type coords: list of floats
        :param coords: [lat_min, lat_max, lon_min, lon_max]
        :type depth_min: float
        :param depth_min: minimum depth of event in km
        :type depth_max: float
        :param depth_max: maximum depth of event in km
        :type mag_min: float
        :param mag_min: minimum magnitude
        :type mag_max: float
        :param mag_max: maximum magnitude
        :type starttime: obspy.UTCDateTime()
        :param starttime: minimum origintime of event
        :type endtime: obspy.UTCDateTime()
        :param endtime: maximum origintime of event
        """
        # Non-mutable conversion for coords
        if not coords:
            coords = []

        # Determine list of events that fall outside box
        events_outside = []
        for event in self.srcrcv.keys():
            # Skip over receivers
            if "." in event:
                continue
            # Skip events that fall outside the bouding box
            elif coords and \
                    (self.srcrcv[event]["lat"] < coords[0]) or \
                    (self.srcrcv[event]["lat"] > coords[1]) or \
                    (self.srcrcv[event]["lon"] < coords[2]) or \
                    (self.srcrcv[event]["lon"] > coords[3]):
                events_outside.append(event)
            # Skip events that are outside the given depth bounds
            elif depth_min or depth_max:
                source_depth = self.srcrcv[event]["depth_m"] * 1E-3
                if depth_min and source_depth < depth_min:
                    events_outside.append(event)
                elif depth_max and source_depth > depth_max:
                    events_outside.append(event)
            # Skip events that are outside the given magnitude bounds
            elif mag_min or mag_max:
                if mag_min and self.srcrcv[event]["mag"] < mag_min:
                    events_outside.append(event)
                elif mag_max and self.srcrcv[event]["mag"] > mag_max:
                    events_outside.append(event)
            # Skip events that are outside the given times
            elif starttime or endtime:
                source_time = UTCDateTime(self.srcrcv[event]["time"])
                if starttime and source_time < starttime:
                    events_outside.append(event)
                elif endtime and source_time < endtime:
                    events_outside.append(event)

        # Go through the misfits and remove in place
        print(f"excluding {len(events_outside)} events")
        for event in events_outside:
            del self.misfits[event]
            del self.windows[event]
            del self.srcrcv[event]

        # Regather attribute information
        self._get_info()
        self._get_str()





