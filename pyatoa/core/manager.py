#!/usr/bin/env python3
"""
A class to control workflow and temporarily store and manipulate data
"""
import os
import obspy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pyflex
import pyadjoint
import warnings

from copy import deepcopy
from obspy.signal.filter import envelope

from pyatoa import logger
from pyatoa.core.config import Config
from pyatoa.utils.asdf.add import (add_misfit_windows, add_adjoint_sources,
                                   add_waveforms, add_config)

from pyatoa.utils.asdf.load import load_windows, load_adjsrcs
from pyatoa.utils.form import channel_code

from pyatoa.utils.process import (apply_filter, trim_streams, zero_pad,
                                  match_npts, normalize, stf_convolve,
                                  is_preprocessed)
from pyatoa.utils.srcrcv import gcd_and_baz
from pyatoa.utils.window import reject_on_global_amplitude_ratio

from pyatoa.scripts.load_example_data import load_example_data

from pyatoa.visuals.wave_maker import WaveMaker
from pyatoa.visuals.map_maker import MapMaker


class ManagerError(Exception):
    """
    A class-wide custom exception raised when functions fail gracefully
    """
    pass


class ManagerStats(dict):
    """
    A simple dictionary that can get and set keys as attributes and has a 
    cleaner looking print statement, used for storing internal statistics
    in the Manager class
    """
    def __init__(self):
        self.dataset_id = None 
        self.event_id = None 
        self.inv_name = None 
        self.nwin = None
        self.len_obs = None
        self.len_syn = None
        self.misfit = None
        self.half_dur = None
        self.time_offset_sec = None
        self.standardized = False 
        self.obs_processed = False
        self.syn_processed = False

    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)

    def __str__(self):
        str_ = ""
        for key, value in self.items():
            str_ += f"{key:>15}: {value}\n"
        return str_[:-1]

    def reset(self):
        """Convenience function to reset stats to None"""
        self.__init__()


class Manager:
    """
    Pyatoas core workflow object.

    Manager is the central workflow control object. It calls on mid and
    low level classes to gather data, standardize and preprocess stream objects,
    generate misfit windows, and calculate adjoint sources. Has a variety of
    internal sanity checks to ensure that the workflow stays on the rails.
    """
    def __init__(self, config=None, ds=None, event=None, st_obs=None,
                 st_syn=None, inv=None, windows=None, staltas=None,
                 adjsrcs=None, gcd=None, baz=None):
        """
        Initiate the Manager class with or without pre-defined attributes.

        .. note::

            If `ds` is not given in data can only be provided through init
            or by passing them directly to the Manager. Data will also not be
            saved.

        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: ASDF data set from which to read and write data
        :type event: obspy.core.event.Event
        :param event: An event object containing relevant earthquake information
        :type st_obs: obspy.core.stream.Stream
        :param st_obs: Stream object containing waveforms of observations
        :type st_syn: obspy.core.stream.Stream
        :param st_syn: Stream object containing waveforms of observations
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory that should only contain the station of interest,
            it's relevant channels, and response information
        :type windows: dict of pyflex.Window objects
        :param windows: misfit windows calculated by Pyflex, stored in a
            dictionary based on component naming
        :type adjsrcs: dict of pyadjoint.AdjointSource objects
        :param adjsrcs: adjoint source waveforms stored in dictionaries
        :type gcd: float
        :param gcd: great circle distance between source and receiver in km
        :type baz: float
        :param baz: Backazimuth between source and receiver in units of degrees
        """
        self.ds = ds
        self.inv = inv

        # Instantiate a Config object
        if config is not None:
            self.config = config
        else:
            logger.info("`config` not provided, initiating empty")
            self.config = Config()

        # Ensure any user-provided event is an Event object
        if isinstance(event, obspy.core.event.catalog.Catalog):
            logger.info(f"event given as catalog, taking zeroth entry")
            event = event[0]
        self.event = event

        # Try to get origin time information from the event
        if self.event is not None:
            origintime = self.event.preferred_origin().time
        else:
            origintime = None

        # Copy Streams to avoid affecting original data
        if st_obs is not None:
            self.st_obs = st_obs.copy()
        else:
            self.st_obs = None
        if st_syn is not None:
            self.st_syn = st_syn.copy()
        else:
            self.st_syn = None

        # Data produced by the workflow
        self.gcd = gcd
        self.baz = baz
        self.windows = windows
        self.staltas = staltas or {}
        self.adjsrcs = adjsrcs
        self.rejwins = {}

        # Internal statistics to keep track of the workflow progress
        self.stats = ManagerStats()

        # Run internal checks on data
        self.check()

    def __str__(self):
        """
        Print statement shows available data detailing workflow
        """
        self.check()
        return ("Manager Data\n"
                f"    dataset   [ds]:        {self.stats.dataset_id}\n"
                f"    quakeml   [event]:     {self.stats.event_id}\n"
                f"    station   [inv]:       {self.stats.inv_name}\n"
                f"    observed  [st_obs]:    {self.stats.len_obs}\n"
                f"    synthetic [st_syn]:    {self.stats.len_syn}\n"
                "Stats & Status\n"  
                f"    half_dur:              {self.stats.half_dur}\n"
                f"    time_offset_sec:       {self.stats.time_offset_sec}\n"
                f"    standardized:          {self.stats.standardized}\n"
                f"    obs_processed:         {self.stats.obs_processed}\n"
                f"    syn_processed:         {self.stats.syn_processed}\n"
                f"    nwin   [windows]:      {self.stats.nwin}\n"
                f"    misfit [adjsrcs]:      {self.stats.misfit}\n"
                )

    def __repr__(self):
        return self.__str__()

    @property
    def st(self):
        """
        Simplified call to return all streams available, observed and synthetic
        """
        if self.st_syn and self.st_obs:
            return self.st_syn + self.st_obs
        elif self.st_syn:
            return self.st_syn
        elif self.st_obs:
            return self.st_obs
        else:
            return None

    def check(self):
        """
        (Re)check the stats of the workflow and data within the Manager.

        Rechecks conditions whenever called, incase something has gone awry
        mid-workflow. Stats should only be set by this function.
        """
        # Give dataset filename if available
        if self.stats.dataset_id is None and self.ds is not None:
            self.stats.dataset_id = os.path.basename(self.ds.filename)

        # Determine the resource identifier for the Event object
        if self.stats.event_id is None and self.event is not None:
            self.stats.event_id = self.event.resource_id.id

        # Get the network and station name from the Inventory object
        if self.stats.inv_name is None and self.inv is not None:
            self.stats.inv_name = ".".join([self.inv[0].code,
                                            self.inv[0][0].code])

        # Check if waveforms are Stream objects, and if preprocessed
        if self.st_obs is not None:
            self.stats.len_obs = len(self.st_obs)
            self.stats.obs_processed = is_preprocessed(self.st_obs)
            if self.stats.len_obs > len(self.config.component_list):
                logger.warning("More observed traces than listed components, "
                               "this may need to be reviewed manually")

        if self.st_syn is not None:
            self.stats.len_syn = len(self.st_syn)
            self.stats.syn_processed = is_preprocessed(self.st_syn)
            if self.stats.len_syn > len(self.config.component_list):
                logger.warning("More synthetic traces than listed components, "
                               "this may need to be reviewed manually")

        # Check that component list matches components in streams, else some
        # functions that rely on Stream.select() will fail to execute
        if self.st_obs is not None and self.stats.obs_processed:
            st_obs_comps = [tr.stats.component for tr in self.st_obs]
            if not set(st_obs_comps).issubset(set(self.config.component_list)):
                logger.warning(f"`st_obs` components {set(st_obs_comps)} != "
                               f"component list {self.config.component_list}")

        if self.st_syn is not None and self.stats.syn_processed:
            st_syn_comps = [tr.stats.component for tr in self.st_syn]
            if not set(st_syn_comps).issubset(set(self.config.component_list)):
                logger.warning(f"`st_syn` components {set(st_syn_comps)} != "
                               f"component list {self.config.component_list}")

        # Check standardization by comparing waveforms against the first
        if not self.stats.standardized and self.st_obs and self.st_syn:
            for tr in self.st[1:]:
                for atr in ["sampling_rate", "npts", "starttime"]:
                    if getattr(tr.stats, atr) != getattr(self.st[0].stats, atr):
                        break
                break
            else:
                self.stats.standardized = True

        # Check for half duration used for source-time-function with synthetics
        if not self.stats.half_dur and self.event is not None:
            try:
                mt = self.event.preferred_focal_mechanism().moment_tensor
                self.stats.half_dur = mt.source_time_function.duration / 2
            except AttributeError:
                pass

        # Count how many misfit windows are contained in the dataset
        if self.stats.nwin is None and self.windows is not None:
            self.stats.nwin = sum([len(_) for _ in self.windows.values()])

        # Determine the unscaled misfit
        if not self.stats.misfit and self.adjsrcs is not None:
            self.stats.misfit = sum([_.misfit for _ in self.adjsrcs.values()])

    def reset(self):
        """
        Restart workflow by deleting all collected data in the Manager, but
        retain dataset, event, config, so a new station can be
        processed with the same configuration as the previous workflow.
        """
        self.__init__(ds=self.ds, event=self.event, config=self.config)

    def write_to_dataset(self, ds=None, choice=None):
        """
        Write the data collected inside Manager to an ASDFDataSet

        :type ds: pyasdf.asdf_data_set.ASDFDataSet or None
        :param ds: write to a given ASDFDataSet. If None, will look for
            internal attribute `self.ds` to write to. Allows overwriting to
            new datasets
        :type choice: list or None
        :param choice: choose which internal attributes to write, by default
            writes all of the following:
            'event': Event atttribute as a QuakeML
            'inv':  Inventory attribute as a StationXML
            'st_obs': Observed waveform under tag `config.observed_tag`
            'st_syn': Synthetic waveform under tag `config.synthetic_tag`
            'windows': Misfit windows collected by Pyflex are stored under
              `auxiliary_data.MisfitWindow`
            'adjsrcs': Adjoint sources created by Pyadjoint are stored under
              `auxiliary_data.AdjointSources`
            'config': the Pyatoa Config object is stored under
                'auxiliary_data.Config' and can be used to re-load the Manager
                and re-do processing
        """
        # Allow using both default and input datasets for writing. Check if we
        # can actually write to the dataset
        if ds is None:
            ds = self.ds
        if ds is None:
            logger.warning("no dataset found, cannot write")
            return
        elif ds._ASDFDataSet__file.mode == "r":  # NOQA
            logger.warning("dataset opened in read-only mode, cannot write")
            return

        if choice is None:
            choice = ["event", "inv", "st_obs", "st_syn", "windows", "adjsrcs",
                      "config"]
            
        logger.info(f"saving Manager attributes to ASDFDataSet")

        # Events should only need to be added once to the ASDFDataSet
        if self.event and "event" in choice:
            try:
                ds.add_quakeml(self.event)
                logger.info("Event object added to ASDFDataSet")
            except ValueError:
                pass
        
        # StationXML files only need to be added once to the ASDFDataSet
        if self.inv and "inv" in choice:
            try:
                ds.add_stationxml(self.inv)
                logger.info("StationXML object added to ASDFDataSet")
            except TypeError:
                pass

        # Observed waveforms only need to be added once to the ASDFDataSet,
        # do not overwrite existing waveforms because observed shouldn't change
        if self.st_obs and "st_obs" in choice:
            add_waveforms(st=self.st_obs, ds=ds, tag=self.config.observed_tag,
                          overwrite=False)

        # Synthetic waveforms should be allowed to overwrite, e.g., in the case
        # that we are rerunning synthetics or restarting iterations
        if self.st_syn and "st_syn" in choice:
            add_waveforms(st=self.st_syn, ds=ds, tag=self.config.synthetic_tag,
                          overwrite=True)

        # Windows will overwrite if windows already exist for this evaluation
        if self.windows and "windows" in choice:
            add_misfit_windows(windows=self.windows, ds=ds, 
                               path=self.config.aux_path, overwrite=True)

        # AdjointSources will overwrite if they already exist for evaluation
        if self.adjsrcs and "adjsrcs" in choice:
            add_adjoint_sources(adjsrcs=self.adjsrcs, ds=ds,
                                path=self.config.aux_path,
                                time_offset=self.stats.time_offset_sec,
                                overwrite=True)

        if self.config and "config" in choice:
            add_config(config=self.config, ds=ds, path=self.config.aux_path)
 
    def write_adjsrcs(self, path="./", write_blanks=True):
        """
        Write internally stored adjoint source traces into SPECFEM defined
        two-column ascii files. Filenames are based on what is expected by
        Specfem, that is: 'NN.SSS.CCC.adj'

        ..note::
            By default writes adjoint sources for ALL components if one
            component has an adjoint source. If an adjoint sourced doesn't exist
            for a given component, it will be written with zeros. This is to
            satisfy SPECFEM3D requirements.

        :type path: str
        :param path: path to save the
        :type write_blanks: bool
        :param write_blanks: write zeroed out adjoint sources for components
            with no adjoint sources to meet the requirements of SPECFEM3D.
            defaults to True
        """
        assert(self.adjsrcs is not None), f"No adjoint sources to write"

        for adj in self.adjsrcs.values():
            fid = f"{adj.network}.{adj.station}.{adj.component}.adj"
            adj.write(filename=os.path.join(path, fid), format="SPECFEM",
                      time_offset=self.stats.time_offset_sec
                      )
        if write_blanks:
            # To see if any blank adjoint sources required, check the difference
            # between internal component list and components with adjsrcs
            # Assumed here that everything is in upper case
            blank_comps = list(
                set(self.config.component_list).difference(
                    set(self.adjsrcs.keys()))
            )
            if blank_comps:
                # Deep copy so that zeroing data doesn't affect original data
                blank_adj = deepcopy(adj)
                blank_adj.adjoint_source *= 0
                for comp in blank_comps:
                    new_adj_comp = f"{adj.component[:-1]}{comp}"
                    fid = f"{adj.network}.{adj.station}.{new_adj_comp}.adj"
                    blank_adj.write(filename=os.path.join(path, fid),
                                    format="SPECFEM",
                                    time_offset=self.stats.time_offset_sec
                                    )

    def load(self, code=None, path=None, ds=None, synthetic_tag=None,
             observed_tag=None, config=True, windows=False,
             adjsrcs=False):
        """
        Populate the manager using a previously populated ASDFDataSet.
        Useful for re-instantiating an existing workflow that has already 
        gathered data and saved it to an ASDFDataSet.

        .. note::
            mgmt.load() will return example data with no dataset

        .. warning::
            Loading any floating point values may result in rounding errors.
            Be careful to round off floating points to the correct place before
            using in future work.

        :type code: str
        :param code: SEED conv. code, e.g. NZ.BFZ.10.HHZ
        :type path: str
        :param path: if no Config object is given during init, the User
            can specify the config path here to load data from the dataset.
            This skips the need to initiate a separate Config object.
        :type ds: None or pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset can be given to load from, will not set the ds
        :type synthetic_tag: str
        :param synthetic_tag: waveform tag of the synthetic data in the dataset
            e.g. 'synthetic_m00s00'. If None given, will use `config` attribute.
        :type observed_tag: str
        :param observed_tag: waveform tag of the observed data in the dataset
            e.g. 'observed'. If None given, will use `config` attribute.
        :type config: bool
        :param config: load config from the dataset, defaults to True but
            can be set False if Config should be instantiated by the User
        :type windows: bool
        :param windows: load misfit windows from the dataset, defaults to False
        :type adjsrcs: bool
        :param adjsrcs: load adjoint sources from the dataset, defaults to False
        """
        if code is None:
            logger.info("loading example data to Manager")
            self.cfg, self.st_obs, self.st_syn, self.event, self.inv = \
                                                            load_example_data() 
        else: 
            # Allows a ds to be provided outside the attribute
            if self.ds and ds is None:
                ds = self.ds
            else:
                raise TypeError("load requires a Dataset")

            # If no Config object in Manager, try to load from dataset
            if config:
                if path is None:
                    raise TypeError("load requires valid 'path' argument")
                logger.info(f"loading config from dataset {path}")
                try:
                    self.config = Config(ds=ds, path=path)
                except AttributeError:
                    logger.warning(f"no Config object in dataset path: {path}")

            assert(self.config is not None), "Config object required for load"
            assert len(code.split('.')) == 2, "'code' must be in form 'NN.SSS'"
            if windows or adjsrcs:
                assert(path is not None), "path required to load auxiliary data"
                iter_, step = path.split("/")

            # Reset and populate using the dataset
            self.__init__(config=self.config, ds=ds, event=ds.events[0])
            net, sta = code.split('.')
            sta_tag = f"{net}.{sta}"
            if sta_tag in ds.waveforms.list():
                self.inv = ds.waveforms[sta_tag].StationXML
                self.st_syn = ds.waveforms[sta_tag][synthetic_tag or
                                                    self.config.synthetic_tag]
                self.st_obs = ds.waveforms[sta_tag][observed_tag or
                                                    self.config.observed_tag]
                if windows:
                    self.windows = load_windows(ds, net, sta, iter_, step, 
                                                return_previous=False)
                if adjsrcs:
                    self.adjsrcs = load_adjsrcs(ds, net, sta, iter_, step)
            else:
                logger.warning(f"no data for {sta_tag} found in dataset")

        self.check()
        return self

    def flow(self, standardize_to="syn", fix_windows=False, iteration=None,
             step_count=None, **kwargs):
        """
        A convenience function to run the full workflow with a single command.
        Does not include gathering. Takes kwargs related to all underlying
        functions.

        .. code:: python

            mgmt = Manager()
            mgmt.flow() == mgmt.standardize().preprocess().window().measure()

        :type standardize_to: str
        :param standardize_to: choice of 'obs' or 'syn' to use one of the time
            series to standardize (resample, trim etc.) the other.
        :type fix_windows: bool
        :param fix_windows: if True, will attempt to retrieve saved windows from
            an ASDFDataSet under the `iteration` and `step_count` tags to use
            during misfit quantification rather than measuring new windows
        :type iteration: int or str
        :param iteration: if 'fix_windows' is True, look for windows in this
            iteration. If None, will check the latest iteration/step_count
            in the given dataset
        :type step_count: int or str
        :param step_count: if 'fix_windows' is True, look for windows in this
            step_count. If None, will check the latest iteration/step_count
            in the given dataset
        :raises ManagerError: for any controlled exceptions
        """
        if fix_windows:
            assert(self.ds is not None), \
                f"`fix_windows` requires ASDFDataSet `ds`"
            assert(iteration is not None and step_count is not None), (
                f"`fix_windows` requires 'iteration' and 'step_count' to access"
                f"windows from ASDFDataSet"
            )

        self.standardize(standardize_to=standardize_to)
        self.preprocess(**kwargs)
        if fix_windows:
            self.retrieve_windows_from_dataset(iteration=iteration,
                                               step_count=step_count)
        else:
            self.window()
        self.measure()

    def flow_multiband(self, periods, standardize_to="syn", fix_windows=False,
                       iteration=None, step_count=None, plot=False, **kwargs):
        """
        Run the full workflow for a number of distinct period bands, returning
        a final set of adjoint sources generated as a summation of adjoint
        sources from each of these period bands. Allows for re-using windows
        collected from the first set of period bands to evaluate adjoint sources
        from the remaining period bands.

        .. note::

            Kwargs are passed through to Manager.preprocess() function only

        .. rubric:: Basic Usage

        Manager.flow_multiband(periods=[(1, 5), (10, 30), (40, 100)])

        :type periods: list of tuples
        :param periods: a list of tuples that define multiple period bands to
            generate windows and adjoint sources for. Overwrites the Config's
            internal `min_period` and `max_period` parameters. The final
            adjoint source will be a summation of all adjoint sources generated.
        :type standardize_to: str
        :param standardize_to: choice of 'obs' or 'syn' to use one of the time
            series to standardize (resample, trim etc.) the other.
        :type fix_windows: bool
        :param fix_windows: if True, will attempt to retrieve saved windows from
            an ASDFDataSet under the `iteration` and `step_count` tags to use
            during misfit quantification rather than measuring new windows
        :type iteration: int or str
        :param iteration: if 'fix_windows' is True, look for windows in this
            iteration. If None, will check the latest iteration/step_count
            in the given dataset
        :type step_count: int or str
        :param step_count: if 'fix_windows' is True, look for windows in this
            step_count. If None, will check the latest iteration/step_count
            in the given dataset
        :type plot: str
        :param plot: name of figure if given, will plot waveform and map for
            each period band and append period band to figure name `plot`
        :rtype: tuple of dict
        :return: (windows, adjoint_sources), returns all the collected
            measurements from each of the period bands
        :raises ManagerError: for any controlled exceptions
        """
        if fix_windows:
            assert(self.ds is not None), \
                f"`fix_windows` requires ASDFDataSet `ds`"
            assert(iteration is not None and step_count is not None), (
                f"`fix_windows` requires 'iteration' and 'step_count' to access"
                f"windows from ASDFDataSet"
            )

        # Copy these waveforms to overwrite for each new period band
        st_obs_raw = self.st_obs.copy()
        st_syn_raw = self.st_syn.copy()

        # Do preprocessing once since
        self.check()
        self.standardize(standardize_to=standardize_to)
        self.preprocess(**kwargs)

        tags = []
        multiband_adjsrcs, multiband_windows = {}, {}
        for period in periods:
            tag = f"{period[0]}-{period[1]}"  # e.g., '5-10s'
            tags.append(tag)
            logger.info(f"calculating adjoint source period band {tag}s")

            self.config.min_period, self.config.max_period = period

            # Standard flow()
            try:
                self.check()
                self.standardize(standardize_to=standardize_to)
                self.preprocess(**kwargs)
                if fix_windows:
                    self.retrieve_windows_from_dataset(iteration=iteration,
                                                       step_count=step_count)
                else:
                    self.window()
                self.measure()
                if plot:
                    save = f"{plot}_{tag}.png"
                    self.plot(choice="both", save=save)
            except ManagerError as e:
                logger.warning(f"period band {tag}s error {e}, "
                               f"cannot return adjoint source for {tag}")

            # Save results of the processing step
            multiband_adjsrcs[tag] = self.adjsrcs or None
            multiband_windows[tag] = self.windows or None

            # Reset for the next run. Don't do a full reset because that gets
            # rid of metadata too, which we need for response removal
            self.windows = None
            self.adjsrcs = None
            self.st_obs = st_obs_raw.copy()
            self.st_syn = st_syn_raw.copy()
            self.stats.reset()

        # Average all adjoint sources into a single object and collect windows
        self.windows, self.adjsrcs = \
            self._combine_mutliband_results(multiband_windows,
                                            multiband_adjsrcs)

    def _combine_mutliband_results(self, windows, adjsrcs):
        """
        Function flow_multiband() generates multiple sets of adjoint sources
        for a variety of period bands, however the User is only interested in a
        single adjoint source which is the average of all of these adjoint
        sources.

        This function will take the multiple sets of adjoint sources and sum
        them accordingly, returning a single set of AdjointSource objects which
        can be used the same as any `adjsrc` attribute returned from `measure`.

        :type adjsrcs: dict of dicts
        :param adjsrcs: a collection of dictionaries whose keys are the
            period band set in `flow_multiband(periods)` and whose values are
            dictionaries returned in `Manager.adjsrcs` from `Manager.measure()`
        :rtype: (dict of Windows, dict of AdjointSource)
        :return: a dictionary of Windows, and AdjointSource objects for each
            component in the componet list. Adjoint sources and misfits
            are the average of all input `adjsrcs` for the given `periods` range
        """
        adjsrcs_out = {}
        windows_out = {}
        for comp in self.config.component_list:
            n = 0
            windows_out[comp] = []
            for tag, adjsrc_dict in adjsrcs.items():
                # Sometimes adjoint sources are empty for a given period range
                # which likely means windows are also empty
                if not adjsrc_dict:
                    continue

                adjsrc = adjsrc_dict[comp]
                # Set a template adjoint source whose attrs. that will change
                if comp not in adjsrcs_out:
                    adjsrcs_out[comp] = deepcopy(adjsrc)
                    adjsrcs_out[comp].min_period = 1E6  # very large number
                    adjsrcs_out[comp].max_period = 0
                    adjsrcs_out[comp].misfit = 0
                    adjsrcs_out[comp].window_stats = []
                    adjsrcs_out[comp].windows = []
                    adjsrcs_out[comp].adjoint_source = adjsrc.adjoint_source * 0

                # Windows are easy, simply append Windows objects to a list
                try:
                    windows_out[comp] += windows[tag][comp]
                except TypeError:
                    pass

                # Set internal attributes as collections of other attributes
                adjsrcs_out[comp].min_period = min(adjsrc.min_period,
                                                   adjsrcs_out[comp].min_period)
                adjsrcs_out[comp].max_period = max(adjsrc.max_period,
                                                   adjsrcs_out[comp].max_period)
                adjsrcs_out[comp].misfit += adjsrc.misfit
                adjsrcs_out[comp].windows += adjsrc.windows
                adjsrcs_out[comp].window_stats += adjsrc.window_stats
                adjsrcs_out[comp].adjoint_source += adjsrc.adjoint_source
                n += 1
            # Normalize based on the number of input adjoint sources
            adjsrcs_out[comp].misfit /= n
            adjsrcs_out[comp].adjoint_source /= n

        return windows_out, adjsrcs_out

    def standardize(self, force=False, standardize_to="syn", normalize_to=None):
        """
        Standardize the observed and synthetic traces in place. 
        Ensures Streams have the same starttime, endtime, sampling rate, npts.

        :type force: bool
        :param force: allow the User to force the function to run even if checks
            say that the two Streams are already standardized
        :type standardize_to: str
        :param standardize_to: allows User to set which Stream conforms to which
            by default the Observed traces should conform to the Synthetic ones
            because exports to Specfem should be controlled by the Synthetic
            sampling rate, npts, etc. Choices are 'obs' and 'syn'.
        :type normalize_to: str
        :param normalize_to: allow for normalizing the amplitudes of the two
            traces. Choices are:
            'obs': normalize synthetic waveforms to the max amplitude of obs
            'syn': normalize observed waveform to the max amplitude of syn
            'one': normalize both waveforms so that their max amplitude is 1
        """
        self.check()
        if not self.stats.len_obs or not self.stats.len_syn:
            raise ManagerError("cannot standardize, not enough waveform data")
        elif self.stats.standardized and not force:
            logger.info("data already standardized")
            return self
        logger.info("standardizing time series")

        # If observations starttime after synthetic, zero pad the front of obs
        dt_st = self.st_obs[0].stats.starttime - self.st_syn[0].stats.starttime
        if dt_st > 0:
            self.st_obs = zero_pad(self.st_obs, dt_st, before=True, after=False)

        # Match sampling rates only if they differ
        if self.st_syn[0].stats.sampling_rate != \
                self.st_obs[0].stats.sampling_rate: 
            if standardize_to == "syn":
                self.st_obs.resample(self.st_syn[0].stats.sampling_rate)
            else:
                self.st_syn.resample(self.st_obs[0].stats.sampling_rate)

        # Match start/endtimes by trimming the 'obs' waveforms to match 'syn'
        self.st_obs, self.st_syn = trim_streams(
            st_a=self.st_obs, st_b=self.st_syn,
            force={"obs": "a", "syn": "b"}[standardize_to]
            )

        # If 'obs' is not long enough, pad end of 'obs' with 0s
        self.st_obs, self.st_syn = match_npts(
            st_a=self.st_obs, st_b=self.st_syn,
            force={"obs": "a", "syn": "b"}[standardize_to]
            )

        # Allow normalization of waveform amplitudes to one another or to
        # a given value
        if normalize_to is not None:
            self.st_obs, self.st_syn = normalize(
                st_a=self.st_obs, st_b=self.st_syn,
                choice={"obs": "a", "syn": "b", "one": "one"}[normalize_to]
            )

        # Determine if synthetics start before the origintime
        if hasattr(self.st_syn[0].stats, "time_offset"):
            self.stats.time_offset_sec = self.st_syn[0].stats.time_offset
        elif self.event is not None:
            self.stats.time_offset_sec = (self.st_syn[0].stats.starttime -
                                          self.event.preferred_origin().time
                                          )
        else:
            logger.warning("cannot find information relating to synthetic time "
                           "offset. Setting to 0")
            self.stats.time_offset_sec = 0
        logger.info(f"syn time offset == {self.stats.time_offset_sec}s")

        # Calculate epicentral distance and backazimuth
        if self.event and self.inv:
            self.gcd, self.baz = gcd_and_baz(event=self.event,
                                             sta=self.inv[0][0])

        self.stats.standardized = True

        return self

    def preprocess(self, which="both", filter_=True, corners=2,
                   remove_response=False, taper_percentage=0.05, zerophase=True,
                   normalize_to=None, convolve_with_stf=True,
                   half_duration=None, **kwargs):
        """
        Apply a simple, default preprocessing scheme to observed and synthetic
        waveforms in place.

        Default preprocessing tasks: Remove response (optional),
        rotate (optional), filter, convolve with source time function (optional)

        User is free to skip this step and perform their own preprocessing on
        `Manager.st_obs` and `Manager.st_syn` if they require their own unique
        processing workflow.

        :type which: str
        :param which: "obs", "syn" or "both" to choose which stream to process
            defaults to "both"
        :type filter_: bool
        :param filter_: filter data using Config.min_period and Config.max_period
            with `corners` filter corners. Apply tapers and demeans before
            and after application of filter.
        :type taper_percentage: float
        :param taper_percentage: percentage [0, 1] of taper to apply to head and
            tail of the data before and after preprocessing
        :type corners: int
        :param corners: number of filter corners to apply if `filter`==True
        :type zerophase: bool
        :param zerophase: apply a zerophase filter (True) or not (False).
            Zerophase filters are run back and forth meaning no phase shift
            is applied, but more waveform distorition may be present.
        :type remove_response: bool
        :param remove_response: flag, remove instrument response from
            'obs' type data using the provided `inv`. Defaults to False.
            Kwargs are passed directly to the the ObsPy `remove_response`
            function. See ObsPy docs for available options.
        :type convolve_with_stf: bool
        :param convolve_with_stf: flag, convolve 'syn' type data with a Gaussian
            source time function to mimic a finite source. Used when half
            half duration in seismic simulations is set to 0.
            Defaults to True and relies on parameters `half_duration`
        :type half_duration: float
        :param half_duration: Source time function half duration in units of
            seconds. Only used if `convolve_with_stf`==True
        :type normalize_to: str
        :param normalize_to: allow for normalizing the amplitudes of the two
            traces. Choices are:
            'obs': normalize synthetic waveforms to the max amplitude of obs
            'syn': normalize observed waveform to the max amplitude of syn
            'one': normalize both waveforms so that their max amplitude is 1
        """
        if which.lower() == "obs":
            preproc_list = {"obs": self.st_obs}
        elif which.lower() == "syn":
            preproc_list = {"syn": self.st_syn}
        else:
            preproc_list = {"obs": self.st_obs, "syn": self.st_syn}

        # Apply preprocessing in-place on streams
        for key, st in preproc_list.items():
            # Remove response from 'obs' type data only
            if remove_response:
                if (key == "obs" and self.config.st_obs_type == "obs") or (
                        key == "syn" and self.config.st_syn_type == "obs"):
                    st.remove_response(inventory=self.inv, plot=False, **kwargs)

            # Set mean to 0 and taper ends to prep for filtering
            st.detrend("demean")
            st.taper(taper_percentage)

            if self.config.rotate_to_rtz and self.baz is not None:
                logger.info(f"rotate {key} NE->RT by {self.baz} degrees")
                st.rotate(method="NE->RT", back_azimuth=self.baz)

            if filter_ and (self.config.min_period or self.config.max_period):
                logger.info(f"filtering {key} {self.config.min_period}--"
                            f"{self.config.max_period}")
                apply_filter(st=st, min_period=self.config.min_period,
                             max_period=self.config.max_period,
                             corners=corners, zerophase=zerophase
                             )
                # Detrend and taper post filter
                st.detrend("simple")
                st.detrend("demean")
                st.taper(taper_percentage)

            # Convolve waveform with source time function for `syn` type data
            if convolve_with_stf and half_duration:
                if (key == "obs" and self.config.st_obs_type == "syn") or (
                        key == "syn" and self.config.st_syn_type == "syn"):
                    stf_convolve(st=st, half_duration=half_duration)

        # Allow normalization of waveform amplitudes to one another or to
        # a given value
        if normalize_to is not None:
            self.st_obs, self.st_syn = normalize(
                st_a=self.st_obs, st_b=self.st_syn,
                choice={"obs": "a", "syn": "b", "one": "one"}[normalize_to]
            )

        # Set stats post preprocessing
        self.stats.obs_processed = is_preprocessed(self.st_obs)
        self.stats.syn_processed = is_preprocessed(self.st_syn)
        self.stats.len_obs = len(self.st_obs)
        self.stats.len_syn = len(self.st_syn)

        return self

    def window(self, windows=None, force=False):
        """
        Evaluate misfit windows using Pyflex. Save windows to ASDFDataSet.
        Allows previously defined windows to be retrieved from ASDFDataSet.

        .. note::
            * Windows are stored as dictionaries of pyflex.Window objects.
            * All windows are saved into the ASDFDataSet, even if retrieved.
            * STA/LTA information is collected and stored internally.

        :type windows: dict
        :param windows: optional argument for User to provide their own windows
            to the window function. This will override the window selection
            process and simply apply the window directly to the class and
            adjust the stats `nwin` for the total number of windows.
        :type force: bool
        :param force: ignore flag checks and run function, useful if e.g.
            external preprocessing is used that doesn't meet flag criteria
        """
        # Pre-check to see if data has already been standardized
        self.check()

        if not self.stats.standardized and not force:
            raise ManagerError("cannot window, waveforms not standardized")

        # Synthetic STA/LTA as Pyflex WindowSelector.calculate_preliminaries()
        for comp in self.config.component_list:
            try:
                self.staltas[comp] = pyflex.stalta.sta_lta(
                    data=envelope(self.st_syn.select(component=comp)[0].data),
                    dt=self.st_syn.select(component=comp)[0].stats.delta,
                    min_period=self.config.min_period
                )
            except IndexError:
                continue

        # If no windows provided, gather windows using waveform data
        if not windows:
            self._select_windows_plus()
        else:
            nwin = 0
            for comp, window in windows.items():
                nwin += len(window)

            self.windows = windows
            self.stats.nwin = sum(len(_) for _ in self.windows.values())

        logger.info(f"{self.stats.nwin} window(s) total found")

        # Print out some window stats for reference
        for comp, windows_ in self.windows.items():
            for w, win in enumerate(windows_):
                logger.debug(f"{comp}_{w}: "
                            f"cc={win.max_cc_value:.2f} / "
                            f"dt={win.cc_shift * win.dt:.2f}s / "
                            f"dlnA={win.dlnA:.2f}")

        return self

    def retrieve_windows_from_dataset(self, ds=None, iteration=None,
                                      step_count=None, components=None, 
                                      revalidate=False):
        """
        Window selection function that retrieves previously saved windows from a
        PyASDF ASDFDataset, recalculates window criteria using the old windows
        with the current data, and attaches window information to Manager.

        :type ds: pyasdf.ASDFDataSet
        :param ds: ASDFDataSet with windows to select. If None given, will
            search for internal definition of `ds`
        :type iteration: int or str
        :param iteration: retrieve windows from the given iteration. If None,
            will search for the previous evaluation to select windows from
        :type step_count: int or str
        :param step_count: retrieve windows from the given step count
            in the given dataset. If None, will search for previous evaluation
            to select windows from
        :type components: list
        :param components: if only windows from certain components should be 
            returned from the dataset. If not given, defaults to Config 
            `component_list`. Should be the inform of a list, e.g., ['N', 'E']
        :type revalidate: bool
        :param revalidate: check acceptability of waveform fit in the 
            retrieved windows, that is, if time shift, dlna or cross correlation
            fall outside of the accepted values defined in the Config object,
            then the window will be rejected directly. IfFalse then
            windows will be returned regardless of their newly assessed misfit.
        """
        if ds is None:
            ds = self.ds
        assert(ds is not None), f"ASDFDataSet `ds` required to retrieve windows"

        # Determine how to treat fixed windows
        if (iteration is None) or (step_count is None):
            # If no iteration/step_count values are given, automatically search
            # the previous step_count for windows in relation to the current
            # iteration/step_count
            iteration = self.config.iteration
            step_count = self.config.step_count
            return_previous = True
        else:
            # If fix windows and iteration/step_count are given, search the
            # dataset for windows under the current iteration/step_count
            return_previous = False

        net, sta, _, _ = self.st_obs[0].get_id().split(".")
        # Function will return empty dictionary if no acceptable windows found
        windows = load_windows(
            ds=ds, net=net, sta=sta, iteration=iteration, step_count=step_count,
            components=components or self.config.component_list,
            return_previous=return_previous
            )

        # Recalculate window criteria for new values for cc, tshift, dlnA etc...
        for comp, windows_ in windows.items():
            # Use Pyflex machinery to re-evaluate the windows based on the
            # current setup of waveforms
            try:
                obs = self.st_obs.select(component=comp)[0]
                syn = self.st_syn.select(component=comp)[0]
                if revalidate:
                    logger.info("revalidating windows against Config criteria")
                    ws = pyflex.WindowSelector(observed=obs, synthetic=syn,
                            config=self.config.pyflex_config, event=self.event,
                            station=self.inv)
                    ws.windows = windows_
                    ws.reject_based_on_data_fit_criteria()
                    windows[comp] = ws.windows

                # If no reject on data fit, simply recalculate the window 
                # criteria and return all windows to User
                else:
                    logger.debug("recalculating window criteria (comp_#):")
                    for w, win in enumerate(windows_):
                        # Log for double check or manual review of new criteria
                        logger.debug(f"{comp}_{w} (old): "
                                    f"cc={win.max_cc_value:.2f} / "
                                    f"dt={win.cc_shift * win.dt:.2f}s / "
                                    f"dlnA={win.dlnA:.2f}")
                        win._calc_criteria(obs.data, syn.data)
                        logger.debug(f"{comp}_{w} (new): "
                                    f"cc={win.max_cc_value:.2f} / "
                                    f"dt={win.cc_shift * win.dt:.2f}s / "
                                    f"dlnA={win.dlnA:.2f}")

            # IndexError thrown when trying to access an empty Stream
            except IndexError:
                continue

        self.windows = windows
        self.stats.nwin = sum(len(_) for _ in self.windows.values())

    def _select_windows_plus(self):
        """
        Mid-level custom window selection function that calls Pyflex select 
        windows, but includes additional window suppression functionality.
        Includes custom Pyflex addition of outputting rejected windows, which
        will be used internally for plotting.

        .. note::

            Pyflex will throw a ValueError if the arrival of the P-wave
            is too close to the initial portion of the waveform, considered the
            'noise' section. This happens for short source-receiver distances
            (< 100km).

            This error becomes a PyflexError if no event/station attributes
            are provided to the WindowSelector

            We could potentially deal with this by zero-padding the
            waveforms, and running select_windows() again, but for now we just
            raise a ManagerError and allow processing to continue
        """
        logger.info(f"windowing waveforms with Pyflex")

        nwin, window_dict, reject_dict = 0, {}, {}
        for comp in self.config.component_list:
            try:
                obs = self.st_obs.select(component=comp)[0]
                syn = self.st_syn.select(component=comp)[0]
            # IndexError thrown when trying to access an empty Stream
            except IndexError:
                continue

            # Pyflex throws a TauP warning from ObsPy #2280, ignore that
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                ws = pyflex.WindowSelector(observed=obs, synthetic=syn,
                                           config=self.config.pyflex_config,
                                           event=self.event,
                                           station=self.inv)
                try:
                    windows = ws.select_windows()
                except (IndexError, pyflex.PyflexError):
                    # see docstring note for why this error is to be addressed
                    raise ManagerError("Cannot window, most likely because "
                                       "the source-receiver distance is too "
                                       "small w.r.t the minimum period")

            # Suppress windows that contain low-amplitude signals
            if self.config.win_amp_ratio > 0:
                windows, ws.rejects["amplitude"] = \
                       reject_on_global_amplitude_ratio(
                                            data=obs.data, windows=windows,
                                            ratio=self.config.win_amp_ratio
                                            )
            # ==================================================================
            # NOTE: Additional windowing criteria may be added here if necessary
            # ==================================================================
            if windows:
                window_dict[comp] = windows
            if ws.rejects:
                reject_dict[comp] = ws.rejects

            # Count windows and tell User
            logger.info(f"{comp}: {len(windows)} window(s)")
            nwin += len(windows)

        self.windows = window_dict
        self.rejwins = reject_dict
        self.stats.nwin = nwin

    def measure(self, force=False):
        """
        Measure misfit and calculate adjoint sources using PyAdjoint.

        Method for caluculating misfit set in Config, Pyadjoint expects
        standardized traces with the same spectral content, so this function
        will not run unless these flags are passed.

        Returns a dictionary of adjoint sources based on component.
        Saves resultant dictionary to a pyasdf dataset if given.

        .. note::

            Pyadjoint returns an unscaled misfit value for an entire set of
            windows. To return a "total misfit" value as defined by 
            Tape (2010) Eq. 6, the total summed misfit will need to be scaled by 
            the number of misfit windows chosen in Manager.window().

        :type force: bool
        :param force: ignore flag checks and run function, useful if e.g.
            external preprocessing is used that doesn't meet flag criteria
        """
        self.check()

        if self.config.adj_src_type is None:
            logger.info("adjoint source type is 'None', will not measure")
            return

        # Check that data has been filtered and standardized
        if not self.stats.standardized and not force:
            raise ManagerError("cannot measure misfit, not standardized")
        elif self.stats.nwin == 0 and not force:
            raise ManagerError("cannot measure misfit, no windows recovered")
        
        logger.debug(f"measuring misfit with adjoint source type: "
                     f"{self.config.adj_src_type}")

        # Create list of windows needed for Pyadjoint
        adjoint_windows = self._format_windows()

        # Run Pyadjoint to retrieve adjoint source objects
        total_misfit, adjoint_sources = 0, {}
        for comp, adj_win in adjoint_windows.items():
            # Streams may not have matching components to given windows, e.g.,
            # during an inversion 
            observed = self.st_obs.select(component=comp)
            synthetic = self.st_syn.select(component=comp)

            if not observed or not synthetic:
                logger.warning(f"no matching observed or synthetic data "
                               f"for component {comp}, cannot measure")
                continue

            # Assuming that only one trace is available per stream
            adj_src = pyadjoint.calculate_adjoint_source(
                config=self.config.pyadjoint_config,
                observed=observed[0], synthetic=synthetic[0],
                windows=adj_win, plot=False
                )

            # Re-format component name to reflect SPECFEM convention
            adj_src.component = f"{channel_code(adj_src.dt)}X{comp}"

            # Save adjoint sources in dictionary object. Sum total misfit
            adjoint_sources[comp] = adj_src
            logger.info(f"{comp} component misfit == {adj_src.misfit:.3f}")
            total_misfit += adj_src.misfit

        # Save adjoint source internally and to dataset
        self.adjsrcs = adjoint_sources

        # Run check to get total misfit
        self.check()
        logger.info(f"total misfit == {self.stats.misfit:.3f}")

        return self

    def _format_windows(self):
        """
        In `pyadjoint.calculate_adjoint_source`, the window needs to be a
        list of lists, with each list containing the [left_window, right_window]
        Each window argument should be given in units of time (seconds).

        :rtype: dict of list of lists
        :return: dictionary with key related to individual components,
            and corresponding to a list of lists containing window start and end
        """
        adjoint_windows = {}

        if self.windows is not None:
            for comp, window in self.windows.items():
                adjoint_windows[comp] = []
                dt = self.st_obs.select(component=comp)[0].stats.delta
                # Prepare Pyflex window indices to give to Pyadjoint
                for win in window:
                    # Window units given in seconds
                    adj_win = [win.left * dt, win.right * dt]
                    adjoint_windows[comp].append(adj_win)
        # If no windows given, calculate adjoint source on whole trace
        else:
            logger.debug("no windows given, adjoint sources will be "
                         "calculated on full trace")
            for comp in self.config.component_list:
                dt = self.st_obs.select(component=comp)[0].stats.delta
                npts = self.st_obs.select(component=comp)[0].stats.npts
                # We offset the bounds of the entire trace by 1s to play nice
                # with PyAdjoints quirky method of generating the adjsrc. 
                # The assumption being the end points will be zero anyway
                adjoint_windows[comp] = [[1, npts * dt-1]]

        return adjoint_windows

    def plot(self, choice="both", save=None, show=True, corners=None,
             figsize=None, dpi=100, **kwargs):
        """
        Plot observed and synthetics waveforms, misfit windows, STA/LTA and
        adjoint sources for all available components. Append information
        about misfit, windows and window selection. Also as subplot create a
        source receiver map which contains annotated information detailing
        src-rcv relationship like distance and BAz. Options to plot either or.

        For valid key word arguments see `visuals.manager_plotter` and
        `visuals.map_maker`

        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :param corners: {lat_min, lat_max, lon_min, lon_max}
            corners to cut the map to, otherwise a global map is provided
        :type choice: str
        :param choice: choice for what to plot:
            * 'wav': plot waveform figure only
            * 'map': plot a source-receiver map only
            * 'both' (default): plot waveform and source-receiver map together
        :type figsize: tuple
        :param figsize: optional size of the figure, set by plot()
        :type dpi: int
        :param dpi: optional dots per inch (resolution) of figure
        """
        self.check()

        # Precheck for correct data to plot
        if choice in ["wav", "both"] and not self.stats.standardized:
            raise ManagerError("cannot plot waveforms, not standardized")

        if choice in ["map", "both"] and (self.inv is None or
                                          self.event is None):
            logger.warning("cannot plot map, no event and/or inv found")
            choice = "wav"

        # Plot only waveform
        if choice == "wav":
            wm = WaveMaker(mgmt=self, **kwargs)
            wm.plot(show=show, save=save)
            plt.close()
        # Plot only map
        elif choice == "map":
            mm = MapMaker(inv=self.inv, cat=self.event, corners=corners,
                          **kwargs)
            mm.plot(show=show, save=save)
            plt.close()
        # Plot waveform and map on the same figure
        elif choice == "both":
            if figsize is None:
                figsize = (1400 / dpi, 600 / dpi)

            # Create an overlying GridSpec that will contain both plots
            gs = mpl.gridspec.GridSpec(1, 2, wspace=0.25, hspace=0.)
            fig = plt.figure(figsize=figsize, dpi=dpi)

            # Plot the waveform on the left
            wm = WaveMaker(mgmt=self)
            wm.plot(figure=fig, subplot_spec=gs[0], show=False, save=False,
                    **kwargs)

            # Plot the map on the right
            mm = MapMaker(cat=self.event, inv=self.inv, figsize=figsize,
                          figure=fig, gridspec=gs, corners=corners,
                          **kwargs)
            mm.plot(figure=fig, gridspec=gs, show=False, save=None)

            if save:
                plt.savefig(save)
            if show:
                plt.show()
            else:
                plt.close()

        # One final shutdown of all figures just incase
        plt.close("all")

