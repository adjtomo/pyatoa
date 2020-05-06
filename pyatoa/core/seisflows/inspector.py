#!/usr/bin/env python3
"""
A class to analyze misfit windows using Pandas.
"""
import os
import pyasdf
import traceback
import pandas as pd
from glob import glob

from pyatoa.utils.form import event_name
from pyatoa.visuals.gadget import Gadget


class Inspector(Gadget):
    """
    This plugin object will collect information from a Pyatoa run folder and
    allow the User to easily understand statistical information or generate
    statistical plots to help understand a seismic inversion

    Inherits plotting capabilities from the Gadget class to reduce clutter.
    """

    def __init__(self, tag=None, path=None):
        """
        Inspector only requires the path to the datasets, it will then read in
        all the datasets and store the data internally. This is a long process
        but should only need to be performed once.

        :type tag: str
        :param tag: tag of a previously saved workflow to be used for reading
            in existing data from disk
        :type path: str
        :param path: path to the ASDFDataSets that were outputted
            by Pyaflowa in the Seisflows workflow
        """
        # If no tag given, create dictionaries based on datasets
        self.windows = pd.DataFrame()
        self.sources = pd.DataFrame()
        self.receivers = pd.DataFrame()

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
                    self.append(dsfid)
                    print("done")
                except KeyError as e:
                    print(f"error: {e}")
                    traceback.print_exc()
                    continue

    def _get_str(self):
        """
        Get the string representation once and save as internal attribute
        """
        # Get a list of internal public methods
        try:
            str_out = (f"{len(self.events):<4} event(s)\n"
                       f"{len(self.stations):<4} station(s)\n"
                       f"{len(self.models):<4} model(s)\n"
                       f"{self.iterations:<4} iteration(s)")

        except KeyError:
            str_out = (f"{0:<4} event(s)\n"
                       f"{0:<4} station(s)\n"
                       f"{0:<4} model(s)\n"
                       f"{0:<4} iteration(s)\n")
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
    def models(self):
        """Return an array of all models"""
        return self._try_print("model")

    @property
    def steps(self):
        """Returns a pandas.Series of models with values listing steps"""
        try:
            return self.windows.groupby("model").apply(
                lambda x: x["step"].unique()
            )
        except KeyError:
            return []

    @property
    def iterations(self):
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
        :rtype source: pandas.DataFrame
        :return source: single row Dataframe containing event info from dataset
        :rtype receivers: multiindexed dataframe containing unique station info
        """
        # Create a dataframe with source information
        event_id = event_name(ds)
        if event_id not in self.sources.index:
            src = {
                "event_id": event_name(ds),
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
            # Append station information one time globally
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

    def _get_windows_from_dataset(self, ds, adj_src_fmt="{net}_{sta}_BX{cmp}"):
        """
        Get window and misfit information from dataset auxiliary data
        Model and Step information should match between the two
        auxiliary data objects MisfitWindows and AdjointSources

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for misfit:
        :type adj_src_fmt: str
        :param adj_src_fmt: adjoint sources will be named different to the
            channel naming they are derived from. This string will be formatted
            in order to access the corresponding AdjointSource auxiliary data.
        :rtype: pandas.DataFrame
        :return: a dataframe object containing information per misfit window
        """
        eid = event_name(ds=ds)

        # Initialize an empty dictionary that will be used to initalize
        # a Pandas DataFrame
        window = {"event": [], "model": [], "step": [], "network": [],
                  "station": [], "channel": [], "component": [], "misfit": [],
                  "length_s": [],
                  }
        winfo = {"dlnA": [], "window_weight": [], "max_cc_value": [],
                 "relative_endtime": [], "relative_starttime": [],
                 "cc_shift_in_seconds": []
                 }

        misfit_windows = ds.auxiliary_data.MisfitWindows
        adjoint_sources = ds.auxiliary_data.AdjointSources

        for model in misfit_windows.list():
            for step in misfit_windows[model].list():
                for win in misfit_windows[model][step]:
                    # pick apart information from this window
                    cha_id = win.parameters["channel_id"]
                    net, sta, loc, cha = cha_id.split(".")
                    component = cha[-1]

                    # get information from corresponding adjoint source
                    # This will be the same for multiple windows
                    adj_src = adjoint_sources[model][step][adj_src_fmt.format(
                        net=net, sta=sta, cmp=component
                    )]
                    window["misfit"].append(adj_src.parameters["misfit_value"])

                    # winfo keys match the keys of the Pyflex Window objects
                    for par in winfo:
                        winfo[par].append(win.parameters[par])

                    # get identifying information for this window
                    window["event"].append(eid)
                    window["network"].append(net)
                    window["station"].append(sta)
                    window["channel"].append(cha)
                    window["component"].append(component)
                    window["model"].append(model)
                    window["step"].append(step)

                    # useful to get window length information
                    window["length_s"].append(
                        win.parameters["relative_endtime"] -
                        win.parameters["relative_starttime"]
                    )
        window.update(winfo)

        self.windows = pd.concat([self.windows, pd.DataFrame(window)],
                                 ignore_index=True)

    def append(self, dsfid, srcrcv=True, windows=True):
        """
        Append a new pyasdf.ASDFDataSet file to the current set of internal
        statistics.

        :type dsfid: str
        :param dsfid: fid of the dataset
        :type srcrcv: bool
        :param srcrcv: gather source-receiver information
        :type windows: bool
        :param windows: gather window information
        """
        try:
            with pyasdf.ASDFDataSet(dsfid) as ds:
                if windows:
                    self._get_windows_from_dataset(ds)
                if srcrcv:
                    self._get_srcrcv_from_dataset(ds)
                return
        except OSError:
            print(f"error: already open")
            return

    def save(self, tag, path="./", fmt="csv"):
        """
        Save the downloaded attributes into JSON files for easier re-loading.

        fmt == 'hdf' requires 'pytables'

        :type tag: str
        :param tag: unique naming tag for saving json files
        :type path: str
        :param path: optional path to save to, defaults to cwd
        :type fmt: str
        :param fmt: format of the files to write, default csv
        """
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

    def write(self, tag, **kwargs):
        """Same as Inspector.save(), but I kept writing .write()"""
        self.save(tag, **kwargs)

    def read(self, tag, path="./", fmt=None):
        """
        Load previously saved attributes to avoid re-processing data.

        :type tag: str
        :param tag: tag to look for json files
        :type path: str
        :param path: optional path to file, defaults to cwd
        :type fmt: str
        :param fmt: format of the files to read, default csv
        """
        # Dynamically determine file format
        if not fmt:
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

    def isolate(self, event=None, model=None, step=None, station=None,
                network=None):
        """
        Returns a new dataframe that is grouped by a given index
        if variable is None, defaults to returning all available values

        :type event: str
        :param event: event id e.g. '2018p130600' (optional
        :type model: str
        :param model: model number e.g. 'm00' (optional)
        :type step: str
        :param step: step count e.g. 's00' (optional)
        :type station: str
        :param station: station name e.g. 'BKZ' (optional)
        :type network: str
        :param network: network name e.g. 'NZ' (optional)
        :rtype: pandas.DataFrame
        :return: DataFrame with selected rows based on selected column values
        """
        df = self.windows
        return df.loc[(df["event"] == (event or df["event"].values)) &
                      (df["model"] == (model or df["model"].values)) &
                      (df["step"] == (step or df["step"].values)) &
                      (df["station"] == (station or df["station"].values)) &
                      (df["network"] == (network or df["network"].values))
                      ]

    def window_lengths(self, level="step"):
        """
        Find the cumulative length of misfit windows for a given model/step,
        or the number of misfit windows for a given model/step.

        :type level: str
        :param level: Default is 'step'
            'step': to get the total window length and number of windows for the
                    given step count.
            'station': to get this on a per-station basis,
                    useful for identifying sta quality.
        :rtype: pandas.DataFrame
        :return: a DataFrame with indices corresponding to model, step,
            columns listing the number of windows (n_win) and the cumulative
            length of windows in seconds (length_s)
        """
        group_list = ["model", "step", "length_s"]
        if level == "station":
            group_list.insert(2, "station")

        windows = self.windows.loc[:, tuple(group_list)]
        windows.sort_values(group_list, inplace=True)

        group = windows.groupby(group_list[:-1]).length_s
        return pd.concat([group.apply(len).rename("n_win"), group.sum()],
                         axis=1)

    def misfits(self, level="step"):
        """
        Sum the total misfit for a given model based on the individual
        misfits for each misfit window, and the number of sources used.

        :type level: str
        :param level:  Default is 'step'
            'step': to get total misfit for a given step count.
            'event': to get this on a per-event misfit.
        :rtype: dict
        :return: total misfit for each model in the class
        """
        group_list = ["model", "step", "event", "misfit"]

        misfits = self.windows.loc[:, tuple(group_list)]
        group = misfits.groupby(group_list[:-1]).misfit
        df = pd.concat([group.sum().rename("unscaled_misfit"),
                        group.apply(len).rename("n_win")], axis=1)

        # Event misfit function a la Tape et al. (2010) Eq. 6
        df["misfit"] = df.apply(
            lambda row: row.unscaled_misfit / (2 * row.n_win), axis=1
        )
        if level == "event":
            pass
        elif level == "step":
            # Sum the event misfits if step-wise misfit is requested
            misfits = df.loc[:, "misfit"]
            group = misfits.groupby(["model", "step"])
            df = pd.concat([group.apply(len).rename("n_event"),
                            group.sum().rename("summed_misfit")], axis=1)
            # Misfit function a la Tape et al. (2010) Eq. 7
            df["misfit"] = df.apply(
                lambda row: row.summed_misfit / row.n_event, axis=1
            )
            df.drop(labels="summed_misfit", axis=1)
        else:
            raise NotImplementedError("level must be 'event' or 'step'")

        return df

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

    def convert_coordinates(self, utm):
        """
        Convert the coordinates of the sources and receivers from the default
        latitude longitude values to a UTM projection of choice

        :param utm:
        :return:
        """



