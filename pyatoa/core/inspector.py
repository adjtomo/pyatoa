#!/usr/bin/env python3
"""
A class to aggregate time windows, source-receiver information and misfit
using Pandas.
"""
import os
import pyasdf
import traceback
import numpy as np
import pandas as pd
from glob import glob
from fnmatch import filter as fnf
from obspy.geodetics import gps2dist_azimuth
from pyatoa.utils.form import format_event_name
from pyatoa.visuals.inspector_plotter import InspectorPlotter


class Inspector(InspectorPlotter):
    """
    This plugin object will collect information from a Pyatoa run folder and
    allow the User to easily understand statistical information or generate
    statistical plots to help understand a seismic inversion.

    Inherits plotting capabilities from InspectorPlotter class to reduce clutter
    """

    def __init__(self, tag="default", verbose=True):
        """
        Inspector will automatically search for relevant file names using the
        tag attribute. If nothing is found, internal dataframes will be empty.

        :type tag: str
        :param tag: tag of a previously saved workflow to be used for reading
            in existing data from disk
        :type verbose: bool
        :param verbose: detail the files that are being read and their status
        """
        self.windows = pd.DataFrame()
        self.sources = pd.DataFrame()
        self.receivers = pd.DataFrame()
        self.tag = tag
        self.verbose = verbose

        # Placeholder attributes for getters
        self._models = None
        self._srcrcv = None

        # Try to load an already created Inspector
        try:
            self.read(tag=self.tag)
        except FileNotFoundError:
            pass

    def _get_str(self):
        """
        Get the string representation once and save as internal attribute
        """
        # Get a list of internal public methods
        try:
            str_out = (f"{len(self.events):<4} event(s)\n"
                       f"{len(self.stations):<4} station(s)\n"
                       f"{len(self.iterations):<4} iteration(s)\n"
                       f"{self.evaluations:<4} evaluation(s)")

        except KeyError:
            str_out = (f"{0:<4} event(s)\n"
                       f"{0:<4} station(s)\n"
                       f"{0:<4} iteration(s)\n"
                       f"{0:<4} evaluation(s)\n")
        return str_out

    def __str__(self):
        """
        Return a list of all variables and functions available for quick ref
        """
        return self._get_str()

    def __repr__(self):
        return self._get_str()

    def _try_print(self, a):
        """Try-except catch for property print statements"""
        try:
            return self.windows.loc[:, a].unique()
        except KeyError:
            try:
                return self.sources.loc[:, a]
            except KeyError:
                return []

    @property
    def events(self):
        """Return an array of all event ids"""
        return self._try_print("event")

    @property
    def stations(self):
        """Return an array of all stations"""
        return self._try_print("station")

    @property
    def networks(self):
        """Return an array of all stations"""
        return self._try_print("network")

    @property
    def netsta(self):
        """Return a Dataframe containing unique network-station idents"""
        try:
            return pd.concat([self.windows.loc[:, "network"],
                              self.windows.loc[:, "station"]],
                             axis=1).drop_duplicates().reset_index(drop=True)
        except KeyError:
            return []

    @property
    def srcrcv(self):
        """Return a dataframe with source-receiver information, dists and baz"""
        if self._srcrcv is None:
            self.get_srcrcv()
        return self._srcrcv

    @property
    def iterations(self):
        """Return an array of all iteration"""
        return self._try_print("iteration")

    @property
    def steps(self):
        """Returns a pandas. Series of iteration with values listing steps"""
        try:
            return self.windows.groupby("iteration").apply(
                lambda x: x["step"].unique()
            )
        except KeyError:
            return []
    
    @property
    def models(self):
        """Return a dict of model numbers related to a unique iteration/step"""
        if self._models is None:
            self.get_models()
        return self._models

    @property
    def initial_model(self):
        """Return tuple of the iteration and step count corresponding M00"""
        if self._models is None:
            self.get_models()
        return self._models.iteration.iloc[0], self._models.step_count.iloc[0]

    @property
    def final_model(self):
        """Return tuple of iteration and step count for final accepted model"""
        if self._models is None:
            self.get_models()
        return self.models.iteration.iloc[-1], self._models.step_count.iloc[-1]

    @property
    def good_models(self):
        """Return models that are only status 0 or 1 (initial or success)"""
        if self._models is None:
            self.get_models()
        return self.models[self.models.status.isin([0, 1])].reset_index(
                                                                    drop=True)

    @property
    def evaluations(self):
        """Returns the number of iterations, or the sum of all step counts"""
        try:
            return sum(self.steps.apply(len).values)
        except AttributeError:
            return 0

    @property
    def mags(self):
        """Return a dictionary of event magnitudes"""
        return self._try_print("magnitude")

    @property
    def times(self):
        """Return a dictionary of event origin times"""
        return self._try_print("time")

    @property
    def depths(self):
        """Return a dictionary of event depths in units of meters"""
        return self._try_print("depth_km")

    def _get_srcrcv_from_dataset(self, ds):
        """
        Get source and receiver information from dataset, this includes
        latitude and longitude values for both, and event information including
        magnitude, origin time, id, etc.

        Returns Dataframes for sources and receivers iff they are not already
        contained in the class dataframes, to avoid duplicates.

        Returns empty DataFrames if no unique info was found.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for distances
        :rtype source: pandas.core.frame.DataFrame
        :return source: single row Dataframe containing event info from dataset
        :rtype receivers: multiindexed dataframe containing unique station info
        """
        # Create a dataframe with source information, ignore duplicates
        event_id = format_event_name(ds)
        if event_id not in self.sources.index:
            src = {
                "event_id": format_event_name(ds),
                "time": str(ds.events[0].preferred_origin().time),
                "magnitude": ds.events[0].preferred_magnitude().mag,
                "depth_km": ds.events[0].preferred_origin().depth * 1E-3,
                "latitude": ds.events[0].preferred_origin().latitude,
                "longitude": ds.events[0].preferred_origin().longitude,
                }
            source = pd.DataFrame([list(src.values())],
                                  columns=list(src.keys())
                                  )
            source.set_index("event_id", inplace=True)

            self.sources = pd.concat([self.sources, source])

        # Loop through all the stations in the dataset to create a dataframe
        networks, stations, latitudes, longitudes = [], [], [], []
        for sta, sta_info in ds.get_all_coordinates().items():
            # Append station information one time globally by checking name
            net, sta = sta.split(".")
            if not (net, sta) in self.receivers.index:
                networks.append(net)
                stations.append(sta)
                latitudes.append(sta_info["latitude"])
                longitudes.append(sta_info["longitude"])

        # Create a list of tuples for multiindexing
        if networks:
            tuples = list(zip(*[networks, stations]))
            idx = pd.MultiIndex.from_tuples(tuples,
                                            names=["network", "station"])
            receivers = pd.DataFrame([latitudes, longitudes],
                                     index=["latitude", "longitude"],
                                     columns=idx
                                     )
            self.receivers = pd.concat([self.receivers, receivers.T])

    def _get_windows_from_dataset(self, ds):
        """
        Get window and misfit information from dataset auxiliary data
        Model and Step information should match between the two
        auxiliary data objects MisfitWindows and AdjointSources

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for misfit:
        :rtype: pandas.DataFrame
        :return: a dataframe object containing information per misfit window
        """
        eid = format_event_name(ds)

        # Initialize an empty dictionary that will be used to initalize
        # a Pandas DataFrame
        window = {"event": [], "iteration": [], "step": [], "network": [],
                  "station": [], "channel": [], "component": [], "misfit": [],
                  "length_s": [],
                  }
        # These are direct parameter names of the MisfitWindow aux data objects
        winfo = {"dlnA": [], "window_weight": [], "max_cc_value": [],
                 "relative_endtime": [], "relative_starttime": [],
                 "cc_shift_in_seconds": [], "absolute_starttime": [],
                 "absolute_endtime": [],
                 }

        misfit_windows = ds.auxiliary_data.MisfitWindows
        adjoint_sources = ds.auxiliary_data.AdjointSources

        for iter_ in misfit_windows.list():
            for step in misfit_windows[iter_].list():
                # If any entries exist for a given event/model/step
                # ignore appending them to the internal structure as they've
                # already been collected
                if not self.windows.empty and \
                        not self.isolate(iter_, step, eid).empty:
                    continue

                for win in misfit_windows[iter_][step]:
                    # pick apart information from this window
                    cha_id = win.parameters["channel_id"]
                    net, sta, loc, cha = cha_id.split(".")
                    component = cha[-1]

                    try:

                        # Workaround for potential mismatch between channel
                        # names of windows and adjsrcs, search for w/ wildcard
                        adj_tag = fnf(adjoint_sources[iter_][step].list(),
                                      f"{net}_{sta}_*{component}"
                                      )[0]

                        # This misfit value will be the same for mult windows
                        window["misfit"].append(adjoint_sources[iter_][step][
                            adj_tag].parameters["misfit"])
                    except IndexError:
                        if self.verbose:
                            print(f"No matching adjoint source for {cha_id}")
                        window["misfit"].append(np.nan)

                    # winfo keys match the keys of the Pyflex Window objects
                    for par in winfo:
                        winfo[par].append(win.parameters[par])

                    # get identifying information for this window
                    window["event"].append(eid)
                    window["network"].append(net)
                    window["station"].append(sta)
                    window["channel"].append(cha)
                    window["component"].append(component)
                    window["iteration"].append(iter_)
                    window["step"].append(step)

                    # useful to get window length information
                    window["length_s"].append(
                        win.parameters["relative_endtime"] -
                        win.parameters["relative_starttime"]
                    )

        # Only add to internal structure if something was collected
        if window["event"]:
            window.update(winfo)
            self.windows = pd.concat([self.windows, pd.DataFrame(window)],
                                     ignore_index=True)

    def discover(self, path="./"):
        """
        Allow the Inspector to scour through a path and find relevant files,
        appending them to the internal structure as necessary.

        :type path: str
        :param path: path to the pyasdf.asdf_data_set.ASDFDataSets that were
            outputted by the Seisflows workflow
        """
        dsfids = glob(os.path.join(path, "*.h5"))
        for i, dsfid in enumerate(dsfids):
            if self.verbose:
                print(
                    f"{os.path.basename(dsfid):<25} {i:0>3}/{len(dsfids):0>3}",
                    end="...")
            try:
                self.append(dsfid)
                if self.verbose:
                    print("done")
            except KeyError as e:
                if self.verbose:
                    print(f"error: {e}")
                    traceback.print_exc()
                continue

        return self

    def append(self, dsfid, srcrcv=True, windows=True):
        """
        Simple function to parse information from a
        pyasdf.asdf_data_setASDFDataSet file and append it to the currect
        collection of information.

        :type dsfid: str
        :param dsfid: fid of the dataset
        :type srcrcv: bool
        :param srcrcv: gather source-receiver information
        :type windows: bool
        :param windows: gather window information
        """
        try:
            with pyasdf.ASDFDataSet(dsfid) as ds:
                if srcrcv:
                    self._get_srcrcv_from_dataset(ds)
                if windows:
                    try:
                        self._get_windows_from_dataset(ds)
                    except AttributeError as e:
                        if self.verbose:
                            print("error reading dataset: "
                                  "missing auxiliary data")
                return
        except OSError:
            if self.verbose:
                print(f"error reading dataset: already open")
            return

    def extend(self, windows):
        """
        Extend the current Inspector data frames with the windows from another
        Inspector. This is useful for when an inversion has been run in legs, so
        two individual inspectors constitute a single inversion.

        .. note::
            The current inspector is considered leg A, and the argument
            'windows' is considered leg B. Leg B will have its iteration numbers
            changed to reflect this

        .. warning::
            This will only work if all the events and stations are the same.
            That is, only two identical inversion scenarios can be used.

        :type windows: pandas.core.data_frame.DataFrame or list of DataFrames
        :param windows: Windows from a separate inspector object that will be
            used to extend the current Inspector. Can also be provided as a list
            of DataFrames to extend multiple times.
        """
        def convert(val):
            """Convenience function to convert between int and str repr"""
            if isinstance(val, str):
                return int(val[1:])
            elif isinstance(val, int):
                return f"i{val:0>2}"

        # To allow for list arguments
        if not isinstance(windows, list):
            windows = [windows]

        for win in windows:
            # Ensure that inplace changes won't affect original data
            windows_ext = win.copy()

            # Determine the new B iteration values based on the
            # final iteration of leg A
            final_iter_a = self.iterations[-1]
            for iter_ in windows_ext.iteration.unique():
                shifted_iter = convert(convert(iter_) + convert(final_iter_a))
                windows_ext.iteration.replace(iter_, shifted_iter, inplace=True)

            self.windows = pd.concat([self.windows, windows_ext])

        # Redo get models since iterations have changed
        if self._models is not None:
            self.get_models()

        return self

    def save(self, path="./", fmt="csv", tag=None):
        """
        Save the downloaded attributes into JSON files for easier re-loading.

        .. note::
            fmt == 'hdf' requires 'pytables' to be installed in the environment

        :type tag: str
        :param tag: tag to use to save files, defaults to the class tag
            but allows for the option of overwriting that
        :type path: str
        :param path: optional path to save to, defaults to cwd
        :type fmt: str
        :param fmt: format of the files to write, default csv
        """
        if tag is None:
            tag = self.tag
        if fmt == "hdf":
            try:
                import pytables
            except ImportError:
                fmt = "csv"
                print("format 'hdf' requires pytables, defaulting to 'csv'")

        if fmt == "csv":
            if not self.sources.empty:
                self.sources.to_csv(os.path.join(path, f"{tag}_src.csv"))
            if not self.receivers.empty:
                self.receivers.to_csv(os.path.join(path, f"{tag}_rcv.csv"))
            if not self.windows.empty:
                self.windows.to_csv(os.path.join(path, f"{tag}.csv"),
                                    index=False)
        elif fmt == "hdf":
            with pd.HDFStore(os.path.join(path, f"{tag}.hdf")) as s:
                s["sources"] = self.sources
                s["receivers"] = self.receivers
                s["windows"] = self.windows
        else:
            raise NotImplementedError

    def write(self, **kwargs):
        """Same as Inspector.save(), but I kept writing .write()"""
        self.save(**kwargs)

    def read(self, path="./", fmt=None, tag=None):
        """
        Load previously saved attributes to avoid re-processing data.

        :type tag: str
        :param tag: tag to use to look for files, defaults to the class tag
            but allows for the option of overwriting that
        :type path: str
        :param path: optional path to file, defaults to cwd
        :type fmt: str
        :param fmt: format of the files to read, default csv
        """
        if tag is None:
            tag = self.tag

        # Dynamically determine file format
        if not fmt:
            tag = tag.split(".")[0]  # remove extension if there is one
            if os.path.exists(os.path.join(path, f"{tag}.csv")):
                fmt = "csv"
            elif os.path.exists(os.path.join(path, f"{tag}.hdf")):
                fmt = "hdf"
            else:
                raise FileNotFoundError

        if fmt == "csv":
            self.sources = pd.read_csv(os.path.join(path, f"{tag}_src.csv"))
            self.sources.set_index("event_id", inplace=True)

            self.receivers = pd.read_csv(os.path.join(path, f"{tag}_rcv.csv"))
            self.receivers.set_index(["network", "station"], inplace=True)

            self.windows = pd.read_csv(os.path.join(path, f"{tag}.csv"))
        elif fmt == "hdf":
            with pd.HDFStore(os.path.join(path, f"{tag}.hdf")) as s:
                self.sources = s["sources"]
                self.receivers = s["receivers"]
                self.windows = s["windows"]
        else:
            raise NotImplementedError

    def reset(self):
        """
        Simple function to wipe out all the internal attributes, not super
        useful but may come in handy somewhere
        """
        self.windows = pd.DataFrame()
        self.sources = pd.DataFrame()
        self.receivers = pd.DataFrame()

    def isolate(self, iteration=None, step_count=None,  event=None,
                network=None, station=None, channel=None, comp=None, keys=None,
                exclude=None, unique_key=None):
        """
        Returns a new dataframe that is grouped by a given index if variable is
        None, defaults to returning all available values

        :type event: str
        :param event: event id e.g. '2018p130600' (optional
        :type iteration: str
        :param iteration: iteration e.g. 'i00' (optional)
        :type step_count: str
        :param step_count: step count e.g. 's00' (optional)
        :type station: str
        :param station: station name e.g. 'BKZ' (optional)
        :type network: str
        :param network: network name e.g. 'NZ' (optional)
        :type channel: str
        :param channel: channel name e.g. 'HHE' (optional)
        :type comp: str
        :param comp: component name e.g. 'Z' (optional)
        :type unique_key: str
        :param unique_key: isolates model, event and station information, 
            alongside a single info key, such as dlnA.
            Useful for looking at one variable without have to write out long 
            lists to 'exclude' or 'keys'
        :type keys: list
        :param keys: list of keys to retain in returned dataset, 'exclude'
            will override this variable, best to use them separately
        :type exclude: list
        :param exclude: list of keys to remove from returned dataset
        :rtype: pandas.DataFrame
        :return: DataFrame with selected rows based on selected column values
        """
        df = self.windows
        df = df.loc[(df["event"] == (event or df["event"].to_numpy())) &
                    (df["iteration"] == (
                                    iteration or df["iteration"].to_numpy())) &
                    (df["step"] == (step_count or df["step"].to_numpy())) &
                    (df["station"] == (station or df["station"].to_numpy())) &
                    (df["network"] == (network or df["network"].to_numpy())) &
                    (df["channel"] == (channel or df["channel"].to_numpy())) &
                    (df["component"] == (comp or df["component"].to_numpy())) 
                    ]
        if unique_key is not None:
            # return the unique key alongside identifying information
            unique_keys = ["event", "iteration", "step", "network", "station", 
                           "channel", "comp", unique_key]
            df = df.loc[:, df.columns.intersection(unique_keys)]
        if exclude is not None:
            # delete excluded keys from key list one by one
            df_keys = df.keys().to_numpy()
            for e in exclude:
                df_keys = df_keys[df_keys != e]
            if keys is not None:
                keys = np.append(df_keys, keys)
            else:
                keys = df_keys
        if keys is not None:
            # 'exclude' may produce repeat keys so run unique beforehand
            df = df.loc[:, df.columns.intersection(np.unique(keys))]
        return df

    def nwin(self, level="step"):
        """
        Find the cumulative length of misfit windows for a given iter/step,
        or the number of misfit windows for a given iter/step.

        .. note::
            Neat trick to select just by station:
            insp.windows(level='station').query("station == 'BFZ'")

        :type level: str
        :param level: Default is 'step'
            'step': to get the total window length and number of windows for the
                    given step count.
            'station': to get this on a per-station basis,
                    useful for identifying sta quality.
        :rtype: pandas.DataFrame
        :return: a DataFrame with indices corresponding to iter, step,
            columns listing the number of windows (n_win) and the cumulative
            length of windows in seconds (length_s)
        """
        group_list = ["iteration", "step", "length_s"]
        if level in ["station", "event"]:
            group_list.insert(2, level)
        elif level == "step":
            pass
        else:
            raise TypeError(
                "nwin() argument 'level' must be 'station', 'event', 'step'")

        windows = self.windows.loc[:, tuple(group_list)]
        windows.sort_values(group_list, inplace=True)

        group = windows.groupby(group_list[:-1]).length_s
        return pd.concat([group.apply(len).rename("nwin"), group.sum()],
                         axis=1)

    def misfit(self, level="step"):
        """
        Sum the total misfit for a given iteration based on the individual
        misfits for each misfit window, and the number of sources used.

        .. note::
            To get per-station misfit on a per-step basis
                df = insp.misfits(level="station").query("station == 'TOZ'")
                df.groupby(['iteration', 'step']).sum()

        :type level: str
        :param level:  Default is 'step'
            'station': unscaled misfit on a per-station basis
            'step': to get total misfit for a given step count.
            'event': to get this on a per-event misfit.
        :rtype: dict
        :return: total misfit for each iteration in the class
        """
        # Various levels to sort the misfit by
        group_list = ["iteration", "step", "event", "station", "component", 
                      "misfit"]
        misfits = self.windows.loc[:, tuple(group_list)]

        # Count the number of windows on a per station basis
        nwin = misfits.groupby(
                group_list[:-1]).misfit.apply(len).rename("nwin")

        # Misfit is unique per component, not window, drop repeat components
        misfits.drop_duplicates(subset=group_list[:-1], keep="first", 
                                inplace=True)

        # Group misfit and window on a per station basis, collect together
        nwin = nwin.groupby(group_list[:-2]).sum()
        misfits = misfits.groupby(
                group_list[:-2]).misfit.sum().rename("unscaled_misfit")
        df = pd.concat([misfits, nwin], axis=1)

        # No formal definition of station misfit so we just define it as the
        # misfit for a given station, divided by number of windows
        if level == "station":
            df["misfit"] = df.apply(
                lambda row: row.unscaled_misfit / row.nwin, axis=1
            )
        # Event misfit function defined by Tape et al. (2010) Eq. 6
        elif level in ["event", "step"]:
            # Group misfits to the event level and sum together windows, misfit
            df = df.groupby(group_list[:3]).sum() 
            df["misfit"] = df.apply(
                lambda row: row.unscaled_misfit / (2 * row.nwin), axis=1
            )
            if level == "step":
                # Sum the event misfits if step-wise misfit is requested
                misfits = df.loc[:, "misfit"]
                group = misfits.groupby(["iteration", "step"])
                df = pd.concat([group.apply(len).rename("n_event"),
                                group.sum().rename("summed_misfit")], axis=1)
                # Misfit function a la Tape et al. (2010) Eq. 7
                df["misfit"] = df.apply(
                    lambda row: row.summed_misfit / row.n_event, axis=1
                )
                df.drop(labels="summed_misfit", axis=1)
        else:
            raise NotImplementedError(
                "level must be 'station', 'event' or 'step'")

        return df

    def compare_misfit(self, iter_init=None, step_init=None, iter_final=None, 
                       step_final=None):
        """
        Calculate the difference in misfit for each source receiver pair,
        comparing one iteration/step count with another.

        .. note::
            Comparison is defined as DELTA = M_FINAL - M_INIT 
            Therefore a negative misfit value corresponds to a reduction/
            improvement in misfit from init to final. Additionally a positive 
            nwin value corresponds to more misfit windows chosen in the final 
            model.

        :type iter_init: str
        :param iter_init: initial iteration to use in comparison
        :type step_init: str
        :param step_init: initial step count to use in comparison
        :type iter_final: str
        :param iter_final: final iteration to use in comparison
        :type step_final: str
        :param step_final: final step count to use in comparison
        :rtype: pandas.core.data_frame.DataFrame
        :return: a sorted data frame containing the difference of misfit and
            number of windows between final and initial
        """
        # Assuming if first arg isnt given, default to first/last model
        if iter_init is None:
            iter_init, step_init = self.initial_model
            iter_final, step_final = self.final_model

        misfit = self.misfit(level="station")
        msft_1 = misfit.loc[iter_init].loc[step_init]
        msft_2 = misfit.loc[iter_final].loc[step_final]

        # Doesn't really make sense to compare unscaled misfit so drop column
        msft_1.drop(["unscaled_misfit"], axis=1, inplace=True)
        msft_2.drop(["unscaled_misfit"], axis=1, inplace=True)

        # Subtract 2 from 1, drop NaN values that occur when no matching info
        # and reverse sort so that improved misfit is shown at the top
        return msft_2.subtract(msft_1).dropna().sort_values(by="misfit",
                                                            ascending=True)

    def filter_sources(self, lat_min=None, lat_max=None, lon_min=None,
                       lon_max=None, depth_min=None, depth_max=None,
                       mag_min=None, mag_max=None, min_start=None,
                       max_start=None):
        """
        Go through misfits and windows and remove events that fall outside
        a certain bounding box. Return sources that fall within the box.
        Bounds are inclusive of given values.

        :type lat_min: float
        :param lat_min: minimum latitude in degrees
        :type lat_max: float
        :param lat_max: maximum latitude in degrees
        :type lon_min: float
        :param lon_min: minimum longitude in degrees
        :type lon_max: float
        :param lon_max: maximum longitude in degrees
        :type depth_min: float
        :param depth_min: minimum depth of event in km, depth is positive
        :type depth_max: float
        :param depth_max: maximum depth of event in km, depth is positive
        :type mag_min: float
        :param mag_min: minimum magnitude
        :type mag_max: float
        :param mag_max: maximum magnitude
        :type min_start: obspy.UTCDateTime()
        :param min_start: minimum origintime of event
        :type max_start: obspy.UTCDateTime()
        :param max_start: maximum origintime of event
        """
        sources = self.sources.copy()
        if lat_min:
            sources = sources.loc[sources["latitude"] >= lat_min]
        if lat_max:
            sources = sources.loc[sources["latitude"] <= lat_max]
        if lon_min:
            sources = sources.loc[sources["longitude"] >= lon_min]
        if lon_max:
            sources = sources.loc[sources["longitude"] <= lon_max]
        if depth_min:
            sources = sources.loc[sources["depth_km"] >= depth_min]
        if depth_max:
            sources = sources.loc[sources["depth_km"] <= depth_max]
        if mag_min:
            sources = sources.loc[sources["magnitude"] >= mag_min]
        if mag_max:
            sources = sources.loc[sources["magnitude"] <= mag_max]
        if min_start or max_start:
            # Convert strings to datetime objects for datetime manipulations
            sources["time"] = pd.to_datetime(sources["time"])
            if min_start:
                sources = sources.loc[
                    sources["time"] >= min_start].set_index("event_id")
            if max_start:
                sources = sources.loc[
                    sources["time"] <= max_start].set_index("event_id")

        return sources

    def get_models(self):
        """
        Return a sorted list of misfits which correspond to accepted models,
        label discards of the line search, and differentiate the final accepted
        line search evaluation from the previous iteration and the initial
        evaluation of the current iteration.

        .. note::
            State and status is given as:
            0 == INITIAL function evaluation for the model;
            1 == SUCCESS -ful function evaluation for the model;
            -1 == DISCARD trial step from line search.

        :rtype: pandas.core.data_frame.DataFrame
        :return: a dataframe containing model numbers, their corresponding
            iteration, step count and misfit value, and the status of the
            function evaluation.
        """
        misfit = self.misfit()
        models = {"model": [], "iteration": [], "step_count": [], "misfit": [],
                  "status": [], "state": []
                  }

        # Model lags iteration by 1
        for m, iter_ in enumerate(self.iterations):
            # First we collect misfit values for each step for reference
            misfits_ = [float(misfit.loc[iter_].loc[_].misfit) for _ in
                        self.steps[iter_]
                        ]

            # Then we loop through the steps and pick out the smallest misfit
            for s, step in enumerate(self.steps[iter_]):
                # Initial evaluation, accepted misfits
                if step == "s00":
                    model = m
                    status = 0
                # Line search, mix of discards and final misfit
                else:
                    model = m + 1
                    if misfits_[s] == min(misfits_):
                        status = 1
                    else:
                        status = -1

                models["model"].append(f"m{model:0>2}")
                models["misfit"].append(misfits_[s])
                models["iteration"].append(iter_)
                models["step_count"].append(step)
                models["state"].append(status)
                models["status"].append({0: "INITIAL", 
                                        1: "SUCCESS", 
                                        -1: "DISCARD"}[status]
                                        )

        self._models = pd.DataFrame(models)

    def get_srcrcv(self):
        """
        Retrieve information regarding source-receiver pairs including distance,
        backazimuth and theoretical traveltimes for a 1D Earth model.

        :rtype: pandas.core.frame.DataFrame
        :return: separate dataframe with distance and backazimuth columns, that
            may be used as a lookup table
        """
        if self.sources.empty or self.receivers.empty:
            return []

        srcrcv_dict = {"event": [], "network": [], "station": [],
                       "distance_km": [], "backazimuth": []
                       }

        for eid, elat, elon, edpth in zip(self.sources.index.to_numpy(),
                                          self.sources.latitude.to_numpy(),
                                          self.sources.longitude.to_numpy(),
                                          self.sources.depth_km.to_numpy()
                                          ):
            for rid, rlat, rlon in zip(self.receivers.index,
                                       self.receivers.latitude.to_numpy(),
                                       self.receivers.longitude.to_numpy()
                                       ):
                gcd, _, baz = gps2dist_azimuth(lat1=elat, lon1=elon,
                                               lat2=rlat, lon2=rlon,
                                               )
                net, sta = rid
                srcrcv_dict["event"].append(eid)
                srcrcv_dict["network"].append(net)
                srcrcv_dict["station"].append(sta)
                srcrcv_dict["distance_km"].append(gcd * 1E-3)
                srcrcv_dict["backazimuth"].append(baz)

        self._srcrcv = pd.DataFrame(srcrcv_dict)







