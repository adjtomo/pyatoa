#!/usr/bin/env python3
"""
A class to aggregate time windows, source-receiver information and misfit
using Pandas.
"""
import os
import pyasdf
import traceback
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from glob import glob
from copy import deepcopy
from fnmatch import filter as fnf
from obspy.geodetics import gps2dist_azimuth
from pyatoa import logger
from pyatoa.utils.form import format_event_name
from pyatoa.visuals.insp_plot import InspectorPlotter


class Inspector(InspectorPlotter):
    """
    This plugin object will collect information from a Pyatoa run folder and
    allow the User to easily understand statistical information or generate
    statistical plots to help understand a seismic inversion.

    Inherits plotting capabilities from InspectorPlotter class to reduce clutter
    """

    def __init__(self, tag="inspector", verbose=False):
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

        # Placeholder attributes for getters
        self._models = None
        self._srcrcv = None
        self._step_misfit = None
        self._event_misfit = None
        self._station_misfit = None

        if verbose:
            logger.setLevel("DEBUG")
        else:
            logger.setLevel("CRITICAL")

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
                       f"{self.evaluations:<4} evaluation(s)"
                       )

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
    def keys(self):
        """Shorthand to access the keys of the Windows dataframe"""
        return self.windows.keys()

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
    def pairs(self):
        """Determine the number of unique source-receiver pairs"""
        cats = ["iteration", "step", "event", "station"]
        df = self.windows.groupby(cats).count()
        # Pick an arbitrary label as all the counts will be the same
        df = df.groupby(cats[:2]).count()[["network"]]

        return df.rename({"network": "count"}, axis=1)

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
        try:
            return self.steps.index[0], self.steps[0][0]
        except TypeError:
            logger.warning("Inspector has no 'steps' data, returning None")
            return None, None

    @property
    def final_model(self):
        """Return tuple of iteration and step count for final accepted model"""
        try:
            return self.steps.index[-1], self.steps[-1][-1]
        except TypeError:
            logger.warning("Inspector has no 'steps' data, returning None")
            return None, None

    @property
    def good_models(self):
        """Return models that are only status 0 or 1 (initial or success)"""
        if self._models is None:
            self.get_models()
        return self.models[self.models.state.isin([0, 1])]

    @property
    def restarts(self):
        """
        Try to guess the indices of restarts for convergence plot based on 
        misfit increase in adjacent good models as well as discontinous misfit 
        values for the final line search model and subsequent initial model.
        Not guaranteed to catch everything so may require manual review using 
        the convergence() function
        """
        if self._models is None:
            self.get_models()

        # Find out where the misfit values increase instead of decrease
        misfit = self.good_models.misfit.round(decimals=3)
        misfit_increase = np.where(misfit.diff() > 0)[0]
        mi_idx = misfit.iloc[misfit_increase].index.values

        # Find out where the same model shows a discontinuous misfit
        dm = self.good_models[["model", "misfit"]].groupby(
                                                 "model").diff().misfit.round(3)
        dm_idx = dm[abs(dm) > 0].index.values

        restart_indices = np.concatenate((mi_idx, dm_idx))

        return self.models.iloc[np.unique(restart_indices)]

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

    def generate_report(self, path_report=None, iteration=None,
                        step_count=None, geographic=True, outliers=True,
                        scatter=True, summary=True, nstd=2,
                        dpi=200, **kwargs):
        """
        An aggregate function that generates a "report" by creating a number
        of figures that summarize the misfit of the inversion. Makes it easier
        for the User as they don't have to remember each of the functions in
        the Inspector's wheelhouse, all relevant figures will be generated
        automatically.

        :type path_report: str
        :param path_report: The path where the report will be saved. 
            Defaults to "./report".
        :type iteration: int
        :param iteration: The iteration number.
        :type step_count: int
        :param step_count: The step count.
        :type geographic: bool
        :param geographic: If True, includes geographic data in the report. 
            Defaults to True.
        :type outliers: bool
        :param outliers: If True, includes outliers in the report. 
            Defaults to True.
        :type scatter: bool
        :param scatter: If True, includes a scatter plot in the report. 
            Defaults to True.
        :type summary: bool
        :param summary: If True, includes a summary in the report. 
            Defaults to True.
        :type nstd: int
        :param nstd: The number of standard deviations for outlier detection.  
            Defaults to 2.
        :type dpi: int
        :param dpi: The resolution in dots per inch for the figures in the 
            report. Defaults to 200.
        :type kwargs: dict
        :param kwargs: Additional keyword arguments.
        """
        # By default we generate a report for the final model 
        iteration, step_count = self.validate_evaluation(iteration, step_count,
                                                         choice="final")
        if path_report is None:
            path_report = f"./report_{iteration}{step_count}"
        if not os.path.exists(path_report):
            os.makedirs(path_report)

        # Generate some geographic information
        if geographic:
            geographic_plot_functions = ["map", "travel_times", "raypaths",
                                         "event_depths"]
            for plot_function in geographic_plot_functions:
                save = os.path.join(path_report, f"{plot_function}.png")
                if os.path.exists(save):
                    continue
                getattr(self, plot_function)(iteration=iteration,
                                             step_count=step_count,
                                             show=False, dpi=dpi)
                plt.close()

        # Plot misfit spider plots of event misfit for events that are outside
        # N standard deviations of the mean w.r.t misfit value
        if outliers:
            upper_outliers, lower_outliers, mean ,std = \
                self.event_outliers(iteration, step_count, nstd=nstd)
            
            if upper_outliers.empty and lower_outliers.empty:
                logger.warning("No outliers found, skipping outlier plots, " 
                               "reduce `nstd` to reevaluate for outliers")
            else:
                for outliers, tag in zip([upper_outliers, lower_outliers],
                                        ["upper_outlier", "lower_outlier"]):
                    for event_name in outliers.index.to_list():
                        self.event_station_misfit_map(
                            event=event_name, iteration=iteration,
                            step_count=step_count, 
                            save=os.path.join(path_report,
                                            f"{tag}_{event_name}.png"),
                            show=False, dpi=dpi
                        )
                        plt.close()

        # Create a few scatterplots comparing some parameters
        if scatter:
            for xy in [
                ("distance_km", "cc_shift_in_seconds"),
                ("backazimuth", "cc_shift_in_seconds"),
                ("length_s", "cc_shift_in_seconds"),
            ]:
                x, y = xy
                self.scatter(x=x, y=y, show=False, dpi=dpi,
                             save=os.path.join(path_report, f"{x}_v_{y}.png")
                             )

        # Plot summary figures that show the status of inversion holistically
        if summary:
            self.convergence(normalize=True, show=False, dpi=dpi,
                             save=os.path.join(path_report, "convergence.png")
                             )

            summary_functions = ["event_station_hist2d", "event_comparison",
                                 "window_stack", "histogram_summary"]
            for plot_function in summary_functions:
                save = os.path.join(path_report, f"{plot_function}.png")
                if os.path.exists(save):
                    continue
                # We want the histogram summary to be comparative
                if plot_function == "histogram_summary":
                    getattr(self, plot_function)(
                        iteration="i01", step_count="s00", 
                        iteration_comp=iteration, step_count_comp=step_count,
                        save=save, show=False, dpi=dpi
                        )
                else:
                    getattr(self, plot_function)(
                        iteration=iteration, step_count=step_count, save=save, 
                        show=False, dpi=dpi
                        )
                plt.close()

        plt.close("all")

        self.generate_report_text(path_report, nstd)
            
    def generate_report_text(self, path_report="./", nstd=1):
        """
        Generate a text report highlighting good/bad performing events and 
        stations that will provide the User a quickly accessible summary of 
        their inversion and may motivate looking at some waveforms
        """       
        line_break = "\n" + "=" * 80 +"\n"     
        iter_end, step_end = self.validate_evaluation(
            iteration=None, step_count=None, choice="final"
            )
        
        # Get event mean and std for current evaluation
        _, _, mean ,std = \
                self.event_outliers(iter_end, step_end, nstd=nstd)
        
        # Get an ascended list of misfit/windows per event
        _windows = self.windows  
        self.windows = self.isolate(iteration=iter_end, step_count=step_end)

        # Event specific window information
        win_per_event = self.nwin(level="event")
        avg_win_per_event = win_per_event.nwin.mean()

        _tenwin = win_per_event.iloc[0].nwin
        top_event_by_win = win_per_event.iloc[0].name[-1]

        _benwin = win_per_event.iloc[-1].nwin
        bot_event_by_win = win_per_event.iloc[-1].name[-1]

        win_per_event_str = win_per_event.to_string()
        
        # Station specific window information
        win_per_station = self.nwin(level="station")
        avg_win_per_station = win_per_station.nwin.mean()

        _tsnwin = win_per_station.iloc[0].nwin
        top_sta_by_win = win_per_station.iloc[0].name[-1]

        _bsnwin = win_per_station.iloc[-1].nwin
        bot_sta_by_win = win_per_station.iloc[-1].name[-1]

        win_per_station_str = win_per_station.to_string()

        # Event specific misfit information
        misfit_per_event = self.misfit(level="event", reset=True)

        _temsft = misfit_per_event.iloc[0].misfit
        highest_misfit_event = misfit_per_event.iloc[0].name[-1]
        
        _bemsft = misfit_per_event.iloc[-1].misfit
        lowest_misfit_event = misfit_per_event.iloc[-1].name[-1]

        misfit_per_event_str = misfit_per_event.sort_values(
            "misfit", ascending=False).to_string()

        # Station specific misfit information
        misfit_per_sta = self.misfit(level="station", reset=True)

        _tsmsft = misfit_per_sta.iloc[0].misfit
        highest_misfit_sta = misfit_per_sta.iloc[0].name[-1]
        
        _bsmsft = misfit_per_sta.iloc[-1].misfit
        lowest_misfit_sta = misfit_per_sta.iloc[-1].name[-1]

        misfit_per_sta_str = misfit_per_sta.sort_values(
            "misfit", ascending=False).to_string()
                
        # Get windows per component for the current evaluation
        _windows_eval = self.windows
        win_str = ""
        for component in self.windows.component.unique():
            self.windows = self.isolate(component=component)
            # Sort of a hacky way of getting what we know is a single value
            nwin_per_comp = self.nwin().nwin.to_list()[0]
            win_str += f"- {component}: {nwin_per_comp}\n"
            self.windows = _windows_eval

        # Restore the original windows, just incase but likely not needed?
        self.windows = _windows  

        # Compile all the above information into a nice text output to be writ
        srcrcv_summary = (
            f"Avg windows per event:  {avg_win_per_event:.2f}\n"
            f"- Evt w/ max win:       {top_event_by_win} ({_tenwin:.0f})\n"
            f"- Evt w/ min win:       {bot_event_by_win} ({_benwin:.0f})\n"
            f"- Evt w/ max msft:      {highest_misfit_event} ({_temsft:.2f})\n"
            f"- Evt w/ min msft:      {lowest_misfit_event} ({_bemsft:.2f})\n"
            "\n"
            f"Avg windows per sta:    {avg_win_per_station:.2f}\n"
            f"- Sta w/ max win:       {top_sta_by_win} ({_tsnwin:.0f})\n"
            f"- Sta w/ min win:       {bot_sta_by_win} ({_bsnwin:.0f})\n"
            f"- Sta w/ max msft:      {highest_misfit_sta} ({_tsmsft:.2f})\n"
            f"- Sta w/ min msft:      {lowest_misfit_sta}  ({_bsmsft:.2f})\n"
            "\n"
            f"Windows per component:\n"
            f"{win_str}"
            )

        # Header contains general information for understanding inversion
        report = [
            f"{'INSPECTOR REPORT':^80}",
            f"{'SUMMARY':^80}",
            f"{self._get_str()}", 
            f"{f'SRCRCV SUMMARY [{iter_end}{step_end}]':^80}",
            f"{srcrcv_summary}",
            f"{'TOTAL WINDOWS':^80}",
            f"{self.nwin().to_string()}",  
            f"{'TOTAL MISFIT':^80}",
            f"{self.misfit().to_string()}", 
            f"{'WINDOWS PER EVENT':^80}",
            f"{win_per_event_str}",  
            f"{f'MISFIT PER EVENT (MEAN={mean:.2f}, {nstd}STD={std:.2f})':^80}",
            f"{misfit_per_event_str}", 
            f"{'WINDOWS PER STATION':^80}",
            f"{win_per_station_str}",
            f"{f'MISFIT PER STATION':^80}",
            f"{misfit_per_sta_str}", 
        ]
    
        with open(os.path.join(path_report, "inspector_report.txt"), "w") as f:
            f.writelines(f"{line_break}".join(report))

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
        event_id = format_event_name(ds.events[0])
        # Some events, like FORCESOLUTIONS, do not contain information on magni.
        try:
            magnitude = ds.events[0].preferred_magnitude().mag
        except AttributeError:
            magnitude = None

        if event_id not in self.sources.index:
            src = {
                "event_id": format_event_name(ds.events[0]),
                "time": str(ds.events[0].preferred_origin().time),
                "magnitude": magnitude,
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
        networks, stations, locations = [], [], []
        latitudes, longitudes = [], []
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

        TODO: break this into _get_windows_from_dataset and 
              _get_adjsrcs_from_dataset?

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to query for misfit:
        :rtype: pandas.DataFrame
        :return: a dataframe object containing information per misfit window
        """
        eid = format_event_name(ds.events[0])

        # Initialize an empty dictionary that will be used to initalize
        # a Pandas DataFrame
        window = {"event": [], "iteration": [], "step": [], "network": [],
                  "station": [], "location": [], "channel": [], "component": [],
                  "misfit": [], "length_s": [],
                  }
        # These are direct parameter names of the MisfitWindow aux data objects
        winfo = {"dlnA": [], "window_weight": [], "max_cc_value": [],
                 "relative_endtime": [], "relative_starttime": [],
                 "cc_shift_in_seconds": [], "absolute_starttime": [],
                 "absolute_endtime": [],
                 }

        misfit_windows = ds.auxiliary_data.MisfitWindows
        adjoint_sources = ds.auxiliary_data.AdjointSources

        # Initiation loop to get iteration and step count, allows for the case
        # where no step count is given (e.g., iteration == 'default')
        iters, steps = [], []
        for iter_ in misfit_windows.list():
            for step in misfit_windows[iter_].list():
                # Ensure that step counts are formatted like: 's00'
                # if not then we DONT have step counts in the dataset
                if not step.startswith("s") and not len(step) == 3:
                    step = ""

                iters.append(iter_)
                steps.append(step)

        # Pulling out important information from the windows and adj src.
        for iter_, step in zip(iters, steps):
            # If any entries exist for a given event/model/step
            # ignore appending them to the internal structure as they've
            # already been collected
            if not self.windows.empty and \
                    not self.isolate(iter_, step, eid).empty:
                continue

            # Explicitely allow for case with no step count in dataset
            misfit_window_eval = misfit_windows[iter_]
            adjoint_source_eval = adjoint_sources[iter_]
            if step:
                misfit_window_eval = misfit_window_eval[step]
                adjoint_source_eval = adjoint_source_eval[step]

            for win in misfit_window_eval:
                # pick apart information from this window
                cha_id = win.parameters["channel_id"]
                net, sta, loc, cha = cha_id.split(".")
                component = cha[-1]

                try:
                    # Workaround for potential mismatch between channel
                    # names of windows and adjsrcs, search for w/ wildcard
                    adj_tag = fnf(adjoint_source_eval.list(),
                                  f"{net}_{sta}_*{component}"
                                  )[0]

                    # This misfit value will be the same for mult windows
                    window["misfit"].append(adjoint_source_eval[
                                                adj_tag].parameters["misfit"])
                except IndexError:
                    logger.warning(f"No matching adjoint source for {cha_id}")
                    window["misfit"].append(np.nan)

                # winfo keys match the keys of the Pyflex Window objects
                for par in winfo:
                    winfo[par].append(win.parameters[par])

                # get identifying information for this window
                window["event"].append(eid)
                window["network"].append(net)
                window["station"].append(sta)
                window["location"].append(loc)
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

    def validate_evaluation(self, iteration, step_count, choice="final"):
        """
        Provide acceptable values for 'iteration' and 'step_count' to underlying
        functions that require it.
        Whenever a user does not choose an iteration or step count, e.g., in
        plotting functions, this function defines default values based on the
        initial model (if neither given), or the last step count for a given
        iteration (if only iteration is given). Only step count is not allowed.
        If both iteration and step count are provided, just check that these
        are acceptable values

        :type iteration: str
        :param iteration: chosen iteration, formatted as e.g., 'i01'
        :type step_count: str
        :param step_count: chosen step count, formatted as e.g., 's00'
        :type choice: str
        :param choice: 'initial' or 'final' to set the default behavior of
            NoneType iteration and step_count returning either the initial
            model or the final model evaluation
        :rtype: tuple of str
        :return: (iteration, step_count) default values for the iteration
            and step_count
        """
        # Default iteration and step count if None are given
        if iteration is None and step_count is None:
            if choice == "initial":
                iteration, step_count = self.initial_model
            elif choice == "final":
                iteration, step_count = self.final_model
            logger.debug(f"No iteration or step count given, defaulting to "
                         f"{choice} model: {iteration}{step_count}")
        elif iteration and (step_count is None):
            step_count = self.steps[iteration][-1]

            logger.debug(f"No step count given, defaulting to final step count "
                         f"within given iteration: {iteration}{step_count}")
        elif (iteration is None) and (step_count is not None):
            raise ValueError("'step_count' cannot be provided by itself, you "
                             "must also set the variable: 'iteration'")
        else:
            assert (iteration in self.iterations and
                    step_count in self.steps[iteration]), \
                f"{iteration}{step_count} does not exist in Inspector"
        return iteration, step_count

    def discover(self, path="./", ignore_symlinks=True):
        """
        Allow the Inspector to scour through a path and find relevant files,
        appending them to the internal structure as necessary.

        :type path: str
        :param path: path to the pyasdf.asdf_data_set.ASDFDataSets that were
            outputted by the Seisflows workflow
        :type ignore_symlinks: bool
        :param ignore_symlinks: skip over symlinked HDF5 files when discovering
        """
        dsfids = glob(os.path.join(path, "*.h5"))
        # remove symlinks from the list if requested
        if ignore_symlinks:
            dsfids = [_ for _ in dsfids if not os.path.islink(_)]
        for i, dsfid in enumerate(dsfids):

            try:
                self.append(dsfid)
                logger.info(f"{os.path.basename(dsfid):<25} "
                            f"{i + 1:0>3}/{len(dsfids):0>3}: done",
                            )
            except KeyError as e:
                logger.info(f"{os.path.basename(dsfid):<25} "
                            f"{i + 1:0>3}/{len(dsfids):0>3}: error {e}",
                            )
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
                        logger.warning("error reading dataset: missing "
                                       "auxiliary data")
                return
        except OSError:
            logger.warning(f"error reading dataset: already open")
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
            write_check = 0
            if not self.sources.empty:
                self.sources.to_csv(os.path.join(path, f"{tag}_src.csv"))
                write_check += 1
            if not self.receivers.empty:
                self.receivers.to_csv(os.path.join(path, f"{tag}_rcv.csv"))
                write_check += 1
            if not self.windows.empty:
                self.windows.to_csv(os.path.join(path, f"{tag}.csv"),
                                    index=False)
                write_check += 1
            if write_check == 0:
                logger.warning("Inspector empty, will not write to disk")
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

    def copy(self):
        """
        Return a deep copy of the Inspector
        """
        return deepcopy(self)

    def reset(self):
        """
        Simple function to wipe out all the internal attributes
        """
        self.windows = pd.DataFrame()
        self.sources = pd.DataFrame()
        self.receivers = pd.DataFrame()

    def isolate(self, iteration=None, step_count=None, event=None,
                network=None, station=None, channel=None, component=None,
                keys=None, exclude=None, unique_key=None):
        """
        Returns a new dataframe that is grouped by a given index if variable is
        None, defaults to returning all available values

        :type event: str
        :param event: event id e.g. '2018p130600' (optional)
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
        :type component: str
        :param component: component name e.g. 'Z' (optional)
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
                    (df["component"] == (
                            component or df["component"].to_numpy()))
                    ]
        if unique_key is not None:
            # return the unique key alongside identifying information
            unique_keys = ["event", "iteration", "step", "network", "station", 
                           "channel", "comp", unique_key]
            df = df.loc[:, df.columns.intersection(unique_keys)]
        if exclude is not None:
            if not isinstance(exclude, list):
                exclude = [exclude]
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
        :param level: Level to get number of windows by. Default is 'step'

            * step: to get the total window length and number of windows for the
              given step count.
            * station: to get this on a per-station basis,
              useful for identifying sta quality.
        :rtype: pandas.DataFrame
        :return: a DataFrame with indices corresponding to iter, step,
            columns listing the number of windows (nwin) and the cumulative
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
        df = pd.concat([group.apply(len).rename("nwin"), group.sum()],
                        axis=1)
        if level == "step":
            return df
        else:
            # Only sort by window number if level is 'station' or 'event'
            return df.sort_values("nwin", ascending=False)

    def misfit(self, level="step", reset=False):
        """
        Sum the total misfit for a given iteration based on the individual
        misfits for each misfit window, and the number of sources used.
        Calculated misfits are stored internally to avoid needing to recalculate
        each time this function is called

        .. note::
            To get per-station misfit on a per-step basis
                df = insp.misfits(level="station").query("station == 'TOZ'")
                df.groupby(['iteration', 'step']).sum()

        :type level: str
        :param level:  Default is 'step'
            'station': unscaled misfit on a per-station basis
            'step': to get total misfit for a given step count.
            'event': to get this on a per-event misfit.
        :type reset: bool
        :param reset: reset internally stored attribute and re-calculate misfit
        :rtype: dict
        :return: total misfit for each iteration in the class
        """
        # We will try to access internal attributes first to save time
        if not reset:
            if level == "step" and self._step_misfit is not None:
                return self._step_misfit
            elif level == "station" and self._station_misfit is not None:
                return self._station_misfit
            elif level == "event" and self._event_misfit is not None:
                return self._event_misfit

        # Various levels to sort the misfit by
        group_list = ["iteration", "step", "event", "network", "station",
                      "component", "misfit"]
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

        # Set internal attribute for easier access at next request
        if level == "step":
            self._step_misfit = df
        elif level == "station":
            self._station_misfit = df
        elif level == "event":
            self._event_misfit = df

        return df

    def stats(self, level="event", choice="mean", key=None, iteration=None,
              step_count=None):
        """
        Calculate the per-level statistical values for DataFrame

        :type level: str
        :param level: get statistical values per 'event' or 'station'
        :type choice: str
        :param choice: Pandas function, 'mean', 'std', 'var', etc.
        :type key: windows column header, e.g. 'cc_shift_in_seconds'
        :type iteration: str
        :param iteration: filter for a given iteration
        :type step_count: str
        :param step_count: filter for a given step count
        :rtype: pandas.DataFrame
        :return: DataFrame containing the `choice` of stats for given options
        """
        group_list = ["iteration", "step", level]

        df = getattr(self.windows.groupby(group_list), choice)(
                numeric_only=True)
        if iteration is not None:
            df = df.loc[iteration]
            if step_count is not None:
                df = df.loc[step_count]
        if key is not None:
            df = df[key]

        return df

    def minmax(self, iteration=None, step_count=None, keys=None,
               quantities=None, pprint=True):
        """
        Calculate and print the min/max values for a whole slew of parameters
        for a given iteration and step count. Useful for understanding the
        worst/ best case scenarios and their relation to the average.

        :type iteration: str
        :param iteration: filter for a given iteration
        :type step_count: str
        :param step_count: filter for a given step count
        :type keys: list of str
        :param keys: keys to calculate minmax values for, must be a subset of
            Inspector.windows.keys()
        :type quantities: list of str
        :param quantities: quantities to get values for, e.g. min, max, median,
            must be an attribute of pandas.core.series.Series
        :type pprint: bool
        :param pprint: pretty print the resulting values
        :rtype: dict
        :return: dictionary containing the minmax stats
        """
        if iteration is None:
            iteration, step_count = self.final_model
        if keys is None:
            keys = ["misfit", "length_s", "dlnA", "max_cc_value",
                    "cc_shift_in_seconds"]
        if quantities is None:
            quantities = ["min", "max", "mean", "median", "std"]

        minmax_dict = {}
        df = self.windows[self.windows.iteration == iteration]
        df = df[df.step == step_count]

        minmax_dict["nwin"] = len(df)
        minmax_dict["len"] = df.length_s.sum()

        for key in keys:
            for quantity in quantities:
                minmax_dict[f"{key}_{quantity}"] = getattr(df[key], quantity)()

        if pprint:
            max_key_len = max([len(_) for _ in minmax_dict.keys()])
            for key, val in minmax_dict.items():
                print(f"{key + ':':<{max_key_len}} {val:.4f}")

        return minmax_dict

    def compare_events(self, iteration_a=None, step_count_a=None, 
                       iteration_b=None, step_count_b=None):
        """
        Compare the misfit and number of windows on an event by event basis
        between two evaluations. Provides absolute values as well as
        differences. Final dataframe is sorted by the difference in misfit,
        showing the most and least improved events.

        :type iteration_a: str
        :param iteration_a: initial iteration to use in comparison
        :type step_count_a: str
        :param step_count_a: initial step count to use in comparison
        :type iteration_b: str
        :param iteration_b: final iteration to use in comparison
        :type step_count_b: str
        :param step_count_b: final step count to use in comparison
        :rtype: pandas.core.data_frame.DataFrame
        :return: a sorted data frame containing the difference of misfit and
            number of windows between final and initial
        """
        # Assuming if first arg isnt given, default to first/last model
        if iteration_a is None:
            iteration_a, step_count_a = self.initial_model
        if iteration_b is None:
            iteration_b, step_count_b = self.final_model

        # If initial or final models not given, nothing to compare
        if None in [iteration_a, step_count_a, iteration_b, step_count_b]:
            logger.warning("Cannot locate model indices to compare model data")
            return None

        misfit = self.misfit(level="event")
        msft_a = misfit.loc[iteration_a, step_count_a]
        msft_b = misfit.loc[iteration_b, step_count_b]

        # Doesn't really make sense to compare unscaled misfit so drop column
        msft_a = msft_a.drop(["unscaled_misfit"], axis=1).copy()
        msft_b = msft_b.drop(["unscaled_misfit"], axis=1).copy()

        # For renaming and access to renamed columns
        initial = f"{iteration_a}{step_count_a}"
        final = f"{iteration_b}{step_count_b}"

        msft_a.rename({"nwin": f"nwin_{initial}",
                       "misfit": f"misfit_{initial}"},
                      axis="columns", inplace=True)
        msft_b.rename({"nwin": f"nwin_{final}", "misfit": f"misfit_{final}"},
                      axis="columns", inplace=True)

        df = pd.merge(msft_a, msft_b, left_index=True, right_index=True)
        df["diff_misfit"] = df[f"misfit_{final}"] - df[f"misfit_{initial}"]
        df["diff_nwin"] = df[f"nwin_{final}"] - df[f"nwin_{initial}"]

        return df.sort_values(by="diff_misfit")

    def compare_windows(self, iteration_a=None, step_count_a=None, 
                        iteration_b=None, step_count_b=None):
        """
        Compare individual, matching misfit windows between two evaluations.
        
        .. note::
            This will only work/make sense if the windows were fixed between 
            the two evaluations, such that they share the exact same window
            selections.

        :type iteration_a: str
        :param iteration_a: initial iteration to use in comparison
        :type step_count_a: str
        :param step_count_a: initial step count to use in comparison
        :type iteration_b: str
        :param iteration_b: final iteration to use in comparison
        :type step_count_b: str
        :param step_count_b: final step count to use in comparison
        :rtype: pandas.core.data_frame.DataFrame
        :return: a data frame containing differences of windowing paramenters
            between final and initial models
        """
        # These are the window values that will be different between two evals
        comp_values = ["misfit", "dlnA", "window_weight", "max_cc_value",
                       "cc_shift_in_seconds"]

        # Assuming if first arg isnt given, default to first/last model
        if iteration_a is None:
            iteration_a, step_count_a = self.initial_model
        if iteration_b is None:
            iteration_b, step_count_b = self.final_model

        # Use copies to ensure any inplace changes don't make it back to self
        windows_a = self.isolate(iteration_a, step_count_a).copy()
        windows_b = self.isolate(iteration_b, step_count_b).copy()
    
        assert(len(windows_a) == len(windows_b)), \
                ("the number of windows does not match between the two "
                 "evaluations, windows cannot be compared")

        # We are using references of the windows to make inplace changes which
        # throws chained assigment warnings. This is acceptable so ignore
        evals = []
        with pd.option_context("mode.chained_assignment", None):
            for _win in [windows_a, windows_b]:
                eval = f"{_win.iteration.iloc[0]}{_win.step.iloc[0]}"
                evals.append(eval)
                # Drop unncessary columns that are not useful in comparison
                _win.drop(["length_s", "relative_endtime", "absolute_starttime",
                           "absolute_endtime", "iteration", "step"],
                             axis=1, inplace=True)
                # Rename columns so they don't get merged into one another
                for column in comp_values:
                    _win.rename({column: f"{column}_{eval}"},
                                       axis="columns", inplace=True)
                # Set the index as a column so that the user can figure out
                # the windows index in the original dataframe
                _win[f"index_{eval}"] = _win.index

        # Merge the evaluations using shared attributes i.e. src rcv info
        df = pd.merge(windows_a, windows_b)
        # Take differences of all the comparison values, 'final - initial'
        initial, final = evals
        for val in comp_values:
            df[f"diff_{val}"] = df[f"{val}_{final}"] - df[f"{val}_{initial}"]

        return df
    
    def compare_misfits(self, iteration_a=None, step_count_a=None,
                        iteration_b=None, step_count_b=None):
        """
        Compare the misfit values between the final and initial model for each
        source receiver pair. Returns differences of unscaled misfit and 
        difference in number of windows between two evaluations. 

        .. note::

            See also compare_events() for a similar function but for events only

        :type iteration_a: str
        :param iteration_a: initial iteration to use in comparison
        :type step_count_a: str
        :param step_count_a: initial step count to use in comparison
        :type iteration_b: str
        :param iteration_b: final iteration to use in comparison
        :type step_count_b: str
        :param step_count_b: final step count to use in comparison
        :rtype: pandas.core.data_frame.DataFrame
        :return: a data frame containing differences of misfit and number of
            windows between final and initial models
        """
        iter_start, step_start = self.validate_evaluation(
            iteration=iteration_a, step_count=step_count_a, choice="initial"
            )
        iter_end, step_end = self.validate_evaluation(
            iteration=iteration_b, step_count=step_count_b, choice="final"
            )
        
        _windows = self.windows  
        self.windows = self.isolate(iteration=iter_start, step_count=step_start)
        misfit_start = self.misfit(level="station", reset=True)

        self.windows = _windows
        self.windows = self.isolate(iteration=iter_end, step_count=step_end)
        misfit_end = self.misfit(level="station", reset=True)

        # Merge the two making sure they match event and station
        sfx_start = f"{iter_start}{step_start}"
        sfx_end = f"{iter_end}{step_end}"
        df = pd.merge(misfit_start, misfit_end, 
                      on=["event", "network", "station"], 
                      suffixes=(f"_{sfx_start}", f"_{sfx_end}")
                      )

        # Calculate the difference between the two misfits
        df["diff_misfit"] = \
            df[f"unscaled_misfit_{sfx_end}"] - df[f"unscaled_misfit_{sfx_start}"]
        df["diff_nwin"] = df[f"nwin_{sfx_end}"] - df[f"nwin_{sfx_start}"]

        # Drop the original values
        df.drop([f"unscaled_misfit_{sfx_end}", f"unscaled_misfit_{sfx_start}", 
                 f"misfit_{sfx_end}", f"misfit_{sfx_start}", 
                 f"nwin_{sfx_end}", f"nwin_{sfx_start}"], 
                axis=1, inplace=True)

        # Reset windows so User can still use the Inspector as advertised
        self.windows = _windows

        # Isolate the largest misfit offenders
        return df.sort_values(by="diff_misfit")
    
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

    def event_outliers(self, iteration=None, step_count=None, choice="misfit",
                       nstd=1):
        """
        Returns outliers for a given misfit measure (misfit or window number)
        by calculating mean and standard deviation and finding events that
        fall outside some integer multiple of standard deviations from the mean.
        Used for plotting in `event_comparison` but also useful for quickly
        assessing which events have anomalously low or high misfit values

        :type iteration: str
        :param iteration: iteration to choose for misfit
        :type step_count: str
        :param step_count: step count to query, e.g. 's00'
        :type choice: str
        :param choice: choice of misfit value, either 'misfit' or 'nwin' or
            'unscaled_misfit'
        :type nstd: int
        :param nstd: number of standard deviations to set upper and lower
            thresholds. Defaults to 1
        :rtype: (Pandas.series, Pandas.series, float, float)
        :return: (events above upper threshold, events below lower threshold,
                  mean, standard deviation)
        """
        iteration, step_count = self.validate_evaluation(iteration, step_count)

        arr = self.misfit(
            level="event")[choice][iteration][step_count].to_numpy()
        index = self.misfit(level="event")[choice][iteration][step_count]
        mean = np.mean(arr)
        std = np.std(arr)

        upper_thresh = mean + (nstd * std)
        lower_thresh = mean - (nstd * std)

        idx_upper_outliers = np.where(arr >= upper_thresh)
        idx_lower_outliers = np.where(arr <= lower_thresh)

        upper_outliers = index.iloc[idx_upper_outliers]
        lower_outliers = index.iloc[idx_lower_outliers]

        return upper_outliers, lower_outliers, mean, std

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

    def get_unique_models(self, float_precision=3):
        """
        Find all accepted models (status 0 or 1) that have a unique misfit
        value. Because some forward evaluations are repeats of the previous
        line search evaluation, they will effectively be the same evaluation so
        they can be removed

        :type float_precision: int
        :param float_precision: identical misfit values will differ after some
            decimal place. this value determines which decimal place to
            truncate the values for comparison
        """
        models = self.good_models
        models.reset_index(drop=True, inplace=True)
        misfit = models.misfit.round(decimals=float_precision)
        identical_misfit = np.where(misfit.diff() == 0)[0]
        models.drop(axis=0, index=identical_misfit, inplace=True)
        models.reset_index(drop=True, inplace=True)

        return models


if __name__ == "__main__":


    """
    Here we define a simple command-line tool for using the Inspector. Useful 
    for Users who have already generated the inspector, and want to quickly 
    make figures to explore the misfit of their inversion. Must be run in the
    directory containing your '.csv' Inspector files. Kwargs can be passed as
    later arguments in the format 'key=val'
    
    .. rubric::
        
        $ python inspector.py <function_name> <kwarg_key=kwarg_val> ...
        e.g.,
        $ python inspector.py event_station_hist2d iteration=1 step_count=2
    """
    import sys
    insp = Inspector()
    kwargs = {}
    if len(sys.argv) > 2:
        for arg in sys.argv[2:]:
            key, val = arg.split("=")
            kwargs[key] = val
    getattr(insp, sys.argv[1])(**kwargs)


