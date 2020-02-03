#!/usr/bin/env python3
"""
Main workflow components of Pyatoa

Manager is the central workflow control object. It calls on mid and low level
classes to gather data; it measures two Obspy stream objects using Pyflex to
generate misfit windows based on parameters set by the Config, and it calculates
adjoint sources using Pyadjoint for misfit quantification.
"""
import warnings

import obspy
import pyflex
import pyadjoint
import traceback
import numpy as np
from os.path import basename
from obspy.signal.filter import envelope

from pyatoa import logger
from pyatoa.core.gatherer import Gatherer
from pyatoa.plugins.pyadjoint_config import src_type

from pyatoa.utils.asdf.additions import write_adj_src_to_asdf
from pyatoa.utils.asdf.extractions import windows_from_ds
from pyatoa.utils.tools.srcrcv import gcd_and_baz, seismogram_length
from pyatoa.utils.tools.format import create_window_dictionary, channel_codes
from pyatoa.utils.tools.calculate import abs_max
from pyatoa.utils.tools.process import preproc, trimstreams, stf_convolve, \
    zero_pad_stream

from pyatoa.utils.visuals.mapping import manager_map
from pyatoa.utils.visuals.waveforms import window_maker


class Manager:
    """
    Core object within Pyatoa.

    Workflow management function that internally calls on all other objects
    within the package in order to gather, process and analyze waveform data.
    """
    def __init__(self, config, ds=None, empty=True, station_code=None,
                 event=None, st_obs=None, st_syn=None, inv=None, windows=None,
                 staltas=None, adj_srcs=None):
        """
        If no pyasdf dataset is given in the initiation of the Manager, all
        data fetching will happen via given pathways in the config file,
        or through external getting via FDSN pathways

        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: ASDF data set from which to read and write data
        :type empty: bool
        :param empty: Do not instantiate Gatherer or look for event.
            Useful for when User provides own event object. if 'event' is given
            then this parameter does not do anything.
        :type station_code: str
        :param station_code: station code for data gather, e.g. 'NZ.BFZ.10.HH?'
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
        :type adj_srcs: dict of pyadjoint.AdjointSource objects
        :param adj_srcs: adjoint source waveforms stored in dictionaries

        """
        # Main workflow requirements
        self.config = config
        self.ds = ds
        self.gatherer = None
        self.station_code = station_code
        self.inv = inv
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
        self.windows = windows
        self.staltas = staltas
        self.adj_srcs = adj_srcs
        # Internal statistics
        self._num_windows = 0
        self._len_obs = 0
        self._len_syn = 0
        self._misfit = 0
        self._time_offset_sec = 0
        self._half_dur = 0
        # Internal flags for workflow status
        self._dataset_id = None
        self._event_name = None
        self._inv_name = None
        self._standardize_flag = False
        self._obs_filter_flag = False
        self._syn_filter_flag = False

        # If event ID set, launch gatherer, gather an event
        if not empty or event is not None:
            self._launch(event=event)
        else:
            # If 'empty' or no event, dont launch gatherer, event is None
            self.event = None

    def __str__(self):
        """
        Print statement shows available data detailing workflow
        """
        self._check()
        return ("DATA\n"
                f"\tdataset (ds):                 {self._dataset_id}\n"
                f"\tevent:                        {self._event_name}\n"
                f"\tmoment tensor (half_dur):     {self._half_dur}\n"
                f"\tinventory (inv):              {self._inv_name}\n"
                f"\tobserved data (st_obs):       {self._len_obs}\n"
                f"\tsynthetic data (st_syn):      {self._len_syn}\n"
                "WORKFLOW\n"
                f"\tstandardized:                 {self._standardize_flag}\n"
                f"\tst_obs filtered:              {self._obs_filter_flag}\n"
                f"\tst_syn filtered:              {self._syn_filter_flag}\n"
                f"\tmisfit windows (windows):     {self._num_windows}\n"
                f"\tmisfit (adj_srcs):            {self._misfit:.2E}\n"
                )
    
    @property
    def st(self):
        """
        Return all streams available in Class
        """
        if isinstance(self.st_syn, obspy.Stream) and \
                isinstance(self.st_obs, obspy.Stream):
            return self.st_syn + self.st_obs
        elif isinstance(self.st_syn, obspy.Stream) and \
                not isinstance(self.st_obs, obspy.Stream):
            return self.st_syn
        elif isinstance(self.st_obs, obspy.Stream) and \
                not isinstance(self.st_syn, obspy.Stream):
            return self.st_obs
        else:
            return None

    @property
    def num_windows(self):
        return self._num_windows

    @property
    def misfit(self):
        return self._misfit

    @property
    def time_offset_sec(self):
        return self._time_offset_sec

    @property
    def half_dur(self):
        return self._half_dur

    def _check(self):
        """
        Update flag information for the User to know location in the workflow.

        NOTE:
        This function rechecks conditions whenever called. This is a safeguard
        incase something has gone awry mid-workflow. Flags could be set and
        forget but that could lead to trouble somewhere down the line.

        Flags should only ever be set by _check() or by the individual functions
        that are allowed to set their own flags, never User.
        """
        # Give dataset filename if available
        if (self.ds is not None) and (self._dataset_id is None):
            self._dataset_id = basename(self.ds.filename)

        # Event as object check, set until reset()
        if (self._event_name is None) and \
                isinstance(self.event, obspy.core.event.Event):
            self._event_name = self.event.resource_id

        # Observed waveforms as Stream objects, and preprocessed
        if isinstance(self.st_obs, obspy.Stream) and len(self.st_obs):
            self._len_obs = len(self.st_obs)
            # Check if a filter has been applied in the processing
            self._obs_filter_flag = (
                    hasattr(self.st_obs[0].stats, "processing") and
                    ("filter(options" in
                     "".join(self.st_obs[0].stats.processing))
            )
        else:
            self._len_obs = 0
            self._obs_filter_flag = False

        # Synthetic waveforms as Stream objects and preprocessed
        if isinstance(self.st_syn, obspy.Stream) and len(self.st_syn):
            self._len_syn = len(self.st_syn)
            # Check if a filter has been applied in the processing
            self._syn_filter_flag = (
                    hasattr(self.st_syn[0].stats, "processing") and
                    ("filter(options" in
                     "".join(self.st_syn[0].stats.processing))
            )
        else:
            self._len_syn = 0
            self._syn_filter_flag = False

        # Standardized waveforms by checking npts and sampling rate.
        # If any of the traces fails the check, the entire flag is False
        # TO DO: also check start and end times?
        if (self.st_obs and self.st_syn) is not None:
            self._standardize_flag = True
            for obs, syn in zip(self.st_obs, self.st_syn):
                if (obs.stats.sampling_rate == syn.stats.sampling_rate) and \
                        (obs.stats.npts == syn.stats.npts):
                    continue
                else:
                    self._standardize_flag = False
                    break
        else:
            self._standardize_flag = False

        # Inventory station check
        if isinstance(self.inv, obspy.Inventory):
            self._inv_name = f"{self.inv[0].code}.{self.inv[0][0].code}"
        else:
            self._inv_name = None

        # Pyflex check if run, return the number of windows made
        self._num_windows = 0
        if isinstance(self.windows, dict):
            for key, win in self.windows.items():
                self._num_windows += len(win)

        # Pyadjoint check if adj_srcs and calculate total misfit
        self._misfit = 0
        if isinstance(self.adj_srcs, dict):
            total_misfit = 0
            for key, adj_src in self.adj_srcs.items():
                total_misfit += adj_src.misfit
            self._misfit = 0.5 * total_misfit / self._num_windows

    def _launch(self, reset=False, event=None, idx=0):
        """
        Appends an event to the Manager class, instantiates low-level Gatherer.

        If an Event or Catalog object is given, passes that to both Manager and
        Gatherer classes.

        If no Event given, queries FDSN via the Gatherer class to search for
        an event.

        :type reset: bool
        :param reset: Reset the Gatherer class for a new run
        :type event: str or obspy.core.event.Event or
            obspy.core.event.catalog.Catalog
        :param event: either an event object to manually set event, or an ID
            to use when searching FDSN for event objects
        :type idx: int
        :param idx: if set event given as Catalog, idx allows user to specify
            which event index in Catalog, defaults to 0.
        """
        # Launch or reset the Gatherer
        if (self.gatherer is None) or reset:
            logger.info("initiating/resetting gatherer")
            self.gatherer = Gatherer(config=self.config, ds=self.ds)

        # If no Event ID is specified in Config, do nothing.
        if self.config.event_id is not None:
            # If the User provides their own Obspy Event or Catalog object
            if event and not isinstance(event, str):
                # If catalog given, take the `idx` entry
                if isinstance(event, obspy.core.event.catalog.Catalog):
                    logger.info(f"event given as catalog, taking entry {idx}")
                    event = event[idx]
                # Populate the Gatherer and Manager with event object
                self.gatherer.set_event(event)
                self.event = self.gatherer.event
            # If no event given by User, turn to Gatherer
            else:
                # If Gatherer has an event, have Manager copy it
                if self.gatherer.event:
                    self.event = self.gatherer.event
                # If Gatherer empty, try gather event
                else:
                    self.event = self.gatherer.gather_event()

    def reset(self, hard_reset=False):
        """
        Delete all collected data in the Manager, to restart workflow.
        Retains Config, ds.

        Soft reset retains event information so that another Station can be
        gathered for the same event, without needing to relaunch Gatherer and
        gather the same event.

        :type hard_reset: bool
        :param hard_reset: hard or soft reset, soft reset doesnt re-instantiate
            gatherer class, and leaves the same event, useful for repeating
            workflow. Hard reset re-initites. default soft
        """
        self.station_code = None
        self.st_obs = None
        self.st_syn = None
        self.inv = None
        self.windows = None
        self.staltas = None
        self.adj_srcs = None
        self._num_windows = 0
        self._len_obs = 0
        self._len_syn = 0
        self._misfit = 0
        self._time_offset_sec = 0
        self._half_dur = 0
        self._dataset_id = None
        self._event_name = None
        self._inv_name = None
        self._standardize_flag = False
        self._obs_filter_flag = False
        self._syn_filter_flag = False

        if hard_reset:
            self._launch(reset=True)
        self._check()

    def gather(self, station_code, choice=None):
        """
        Launch a gatherer object and gather event, station and waveform
        information given a station code. Fills the manager based on information
        most likely to be available (we expect an event to be available more
        often than waveform data).
        Catches general exceptions along the way, stops gathering if errors.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type choice: list
        :param choice: allows user to gather individual bits of data, rather
            than gathering all. Allowed: 'inv', 'st_obs', 'st_syn'
        """
        if self.gatherer is None:
            self._launch()
        try:
            self.station_code = station_code
            logger.info(f"gathering {station_code} for {self.config.event_id}")
            # Gather all data
            if choice is None:
                logger.debug("gathering station information")
                self.inv = self.gatherer.gather_station(station_code)
                logger.debug("gathering observation waveforms")
                self.st_obs = self.gatherer.gather_observed(station_code)
                logger.debug("gathering synthetic waveforms")
                self.st_syn = self.gatherer.gather_synthetic(station_code)
            # Gather specific data based on user defined choice
            else:
                if "inv" in choice:
                    logger.debug("gathering station information")
                    self.inv = self.gatherer.gather_station(station_code)
                if "st_obs" in choice:
                    logger.debug("gathering observation waveforms")
                    self.st_obs = self.gatherer.gather_observed(station_code)
                if "st_syn" in choice:
                    logger.debug("gathering synthetic waveforms")
                    self.st_syn = self.gatherer.gather_synthetic(station_code)
        except obspy.clients.fdsn.header.FDSNNoDataException:
            logger.info("No data found internally or externally")
            return
        except Exception as e:
            traceback.print_exc() 
            return

    def write(self, write_to="ds"):
        """
        Write the data collected inside Manager to either a Pyasdf Dataset,
        or to individual files

        :type write_to: str
        :param write_to: choice to write data to, if "ds" writes to
            Pyasdf Dataset

        write_to == "ds"
            If `gather` is skipped but data should still be saved into an
            ASDFDataSet for data storage, this function will fill that dataset
            in the same fashion as the Gatherer class
        write_to == "/path/to/output"
            write out all the internal data of the manager to a path
        """
        if write_to == "ds":
            # Only populate if all requisite parts are given
            if self.event:
                try:
                    self.ds.add_quakeml(self.event)
                except ValueError:
                    logger.warning("event already present, not added")
            if self.inv:
                try:
                    self.ds.add_stationxml(self.inv)
                except TypeError:
                    logger.warning("inv already present, not added")
            # PyASDF has its own warnings if waveform data already present
            if self.st_obs:
                self.ds.add_waveforms(waveform=self.st_obs,
                                      tag=self.config.observed_tag)
            if self.st_syn:
                self.ds.add_waveforms(waveform=self.st_syn,
                                      tag=self.config.synthetic_tag)
            if self.windows:
                self.save_windows()
            if self.adj_srcs:
                self.save_adj_srcs()
        else:
            raise NotImplementedError

    def load(self, station_code, model, ds=None):
        """
        Populate the manager using a given PyASDF Dataset, based on user-defined
        station code. Useful for re-instantiating an existing workflow that
        has already gathered data and saved it to a dataset.

        :type station_code: str
        :param station_code: SEED conv. code, e.g. NZ.BFZ.10.HHZ
        :type model: str
        :param model: model number to search for data, e.g. 'm00'
        :type ds: None or pyasdf.ASDFDataSet
        :param ds: dataset can be given to load from, will not set the ds
        """
        # Allows a ds to be provided outside the attribute
        if self.ds and ds is None:
            ds = self.ds
        else:
            raise AttributeError("load requires a Dataset")

        assert len(station_code.split('.')) == 4,\
            "station_code must be in form 'NN.SSS.LL.CCC'"

        self.event = ds.events[0]

        net, sta, _, _ = station_code.split('.')
        sta_tag = f"{net}.{sta}"
        if sta_tag in ds.waveforms.list():
            self.st_obs = ds.waveforms[sta_tag][self.config.observed_tag]
            self.st_syn = ds.waveforms[sta_tag][
                self.config.synthetic_tag.format(model)]
            self.inv = ds.waveforms[sta_tag].StationXML

    def standardize(self, force=False, standardize_to="syn"):
        """
        Standardize the observed and synthetic traces in place. Ensure that the
        data streams have the same start and endtimes, and sampling rate, so
        that data comparisons can be made; all preprocessing related to timing
        and sampling rate.

        :type force: bool
        :param force: allow the User to force the functino to run even if checks
            say that the two Streams are already standardized
        :type standardize_to: str
        :param standardize_to: allows User to set which Stream conforms to which
            by default the Observed traces should conform to the Synthetic ones
            because exports to Specfem should be controlled by the Synthetic
            sampling rate, npts, etc.
        """
        self._check()
        if min(self._len_obs, self._len_syn) == 0:
            logger.warning("cannot standardize, not enough waveform data")
            return
        elif self._standardize_flag and not force:
            logger.warning("already standardized")
            return
        logger.info("standardizing streams")

        # Zero pad the data if set by Config
        if self.config.zero_pad:
            self.st_obs = zero_pad_stream(self.st_obs, self.config.zero_pad)
            self.st_syn = zero_pad_stream(self.st_syn, self.config.zero_pad)
            logger.debug(f"zero padding front, back by {self.config.zero_pad}s")

        # Resample one Stream to match the other
        if standardize_to == "syn":
            self.st_obs.resample(self.st_syn[0].stats.sampling_rate)
        else:
            self.st_syn.resample(self.st_obs[0].stats.sampling_rate)

        # Trim observations and synthetics to the length of chosen
        trim_to = {"obs": "a", "syn": "b"}
        self.st_obs, self.st_syn = trimstreams(
            st_a=self.st_obs, st_b=self.st_syn, force=trim_to[standardize_to])

        # Retrieve the first timestamp in the .sem? file from Specfem
        if self.event is not None:
            self._time_offset_sec = (self.st_syn[0].stats.starttime -
                                     self.event.preferred_origin().time
                                     )
        else:
            self._time_offset_sec = 0
        logger.debug(f"time offset set to {self._time_offset_sec}s")

        self._standardize_flag = True

    def preprocess(self, which="both", overwrite=None):
        """
        Standard preprocessing of observed and synthetic data in place.
        Called in identical manner for observation and synthetic waveforms.
        Synthetic waveform is convolved with a source time function.

        This function can of course be overwritten by a User defined function

        :type which: str
        :param which: "obs", "syn" or "both" to choose which stream to process
            defaults to both
        :type overwrite: function
        :param overwrite: If a function is provided, it will overwrite the 
            standard preprocessing function. All arguments that are given
            to the standard preprocessing function will be passed as kwargs to
            the new function. This allows for customized preprocessing
        """
        self._check()
        # Make sure an instrument response is available for removal, or that
        # this is a synthetic-synthetic case
        if (not isinstance(self.inv, obspy.core.inventory.Inventory)) \
                and (not self.config.synthetics_only):
            logger.warning("cannot preprocess, no inventory")
            return
        if overwrite:
            # Ensure the overwrite call is a function
            assert(hasattr(overwrite, '__call__')), "overwrite must be function"
            preproc_fx = overwrite
        else:
            preproc_fx = preproc

        # If required, rotate based on source receiver lat/lon values
        baz = None
        if self.config.rotate_to_rtz:
            _, baz = gcd_and_baz(event=self.event, sta=self.inv[0][0])

        # Preprocess observation and synthetic data the same
        if self.st_obs is not None and not self._obs_filter_flag and \
                which.lower() in ["obs", "both"]:
            # Determine if different preprocessing required for syn-syn case
            if self.config.synthetics_only:
                obs_inv = None
                obs_synthetic_unit = self.config.synthetic_unit
            else:
                obs_inv = self.inv
                obs_synthetic_unit = None
            logger.info("preprocessing observation data")
            self.st_obs = preproc_fx(st_original=self.st_obs, inv=obs_inv,
                                     synthetic_unit=obs_synthetic_unit,
                                     back_azimuth=baz,
                                     unit_output=self.config.unit_output,
                                     corners=self.config.filter_corners,
                                     filter_bounds=[self.config.min_period,
                                                    self.config.max_period],
                                     )

        if self.st_syn is not None and not self._syn_filter_flag and \
                which.lower() in ["syn", "both"]:
            logger.info("preprocessing synthetic data")
            self.st_syn = preproc_fx(st_original=self.st_syn, inv=None,
                                     synthetic_unit=self.config.synthetic_unit,
                                     back_azimuth=baz,
                                     unit_output=self.config.unit_output,
                                     corners=self.config.filter_corners,
                                     filter_bounds=[self.config.min_period,
                                                    self.config.max_period]
                                     )
        
        # Check to see if preprocessing failed
        self._check()
        if not self._obs_filter_flag or not self._syn_filter_flag:
            logger.warning("preprocessing failed")
            return

        # Convolve synthetic data with a gaussian source-time-function
        self._convolve_source_time_function(which)

    def _convolve_source_time_function(self, which="both"):
        """
        Convolve synthetic data with a gaussian source time function, time
        shift by a given half duration.

        TO DO:
            check if time_offset is doing what I want it to do

        :type which: str
        :param which: "obs", "syn" or "both" to choose which stream to process
            defaults to both
        """
        try:
            moment_tensor = self.event.focal_mechanisms[0].moment_tensor
            self._half_dur = moment_tensor.source_time_function.duration / 2
            if which.lower() in ["syn", "both"]:
                self.st_syn = stf_convolve(st=self.st_syn,
                                           half_duration=self._half_dur,
                                           time_shift=False,
                                           # time_offset=self._time_offset_sec
                                           )
            # If a synthetic-synthetic case, convolve observations too
            if self.config.synthetics_only and which in ["obs", "both"]:
                self.st_obs = stf_convolve(st=self.st_obs,
                                           half_duration=self._half_dur,
                                           time_shift=False,
                                           )
        except (AttributeError, IndexError):
            logger.info("moment tensor not found for event, cannot convolve")

    def window(self, fix_windows=False, force=False):
        """
        Call Pyflex to calculate best fitting misfit windows given observation
        and synthetic data. Data must be standardized.

        Save ouputs as dictionaries of window objects, as well as STA/LTAs.
        If a pyasdf dataset is present, save misfit windows as auxiliary data

        :type fix_windows: bool
        :param fix_windows: do not pick new windows, but load windows from the
            given dataset
        :type force: bool
        :param force: ignore flag checks and run function, useful if e.g.
            external preprocessing is used that doesn't meet flag criteria
        """
        # Pre-check to see if data has already been standardized
        self._check()
        if not self._standardize_flag and not force:
            logger.warning("cannot window, waveforms not standardized")
            return
        if fix_windows and not self.ds:
            logger.warning("cannot fix window, no windows in dataset")
            fix_windows = False

        # Get STA/LTA information
        staltas = {}
        for comp in self.config.component_list:
            try:
                staltas[comp] = pyflex.stalta.sta_lta(
                    dt=self.st_syn.select(component=comp)[0].stats.delta,
                    min_period=self.config.min_period,
                    data=envelope(self.st_syn.select(component=comp)[0].data)
                )
            except IndexError:
                continue
        self.staltas = staltas

        # Get misfit windows
        # If windows are to be fixed, ensure that there are windows in dataset
        if fix_windows and (hasattr(self.ds, "auxiliary_data") and
                            hasattr(self.ds.auxiliary_data, "MisfitWindows")):
            net, sta, _, _ = self.st_obs[0].get_id().split(".")
            try:
                self.windows, self._num_windows = windows_from_ds(self.ds, net, 
                                                                  sta)
            except AttributeError:
                self.windows, self._num_windows = self.select_windows()
        # If not fixed windows, calculate windows using Pyflex
        else:
            # Windows and staltas saved as dictionary objects by component name
            self.windows, self._num_windows = self.select_windows()

        # Save to dataset only if new windows are picked
        if self.ds is not None and (self._num_windows != 0) and not fix_windows:
            self.save_windows()

        # Let the User know the outcomes of Pyflex
        logger.info(f"{self._num_windows} window(s) total found")

    def select_windows(self):
        """
        Custom window selection function to include suppression by amplitude

        :type comp: str
        :param comp: component to select waveform data by
        :rtype window: pyflex.Window
        :return window: the windows calculated by Pyflex and filtered by amp rat
        """
        logger.info(f"running Pyflex w/ map: {self.config.pyflex_map}")

        num_windows, windows = 0, {}
        for comp in self.config.component_list:
            try:
                # Pyflex throws a TauP warning from ObsPy #2280, ignore that
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning) 
                    # Run Pyflex to select misfit windows as list of Window obj
                    window = pyflex.select_windows(
                        observed=self.st_obs.select(component=comp),
                        synthetic=self.st_syn.select(component=comp),
                        config=self.config.pyflex_config, event=self.event,
                        station=self.inv)

                # Suppress windows that contain signals smaller than some
                # fraction of the peak amplitude contained in the synthetic 
                if window and self.config.window_amplitude_ratio > 0:
                    windows_by_amplitude = []
                    for win_ in window:
                        waveform_peak = abs_max(
                            self.st_syn.select(component=comp)[0].data
                        )
                        window_peak = abs_max(
                            self.st_syn.select(
                                component=comp)[0].data[win_.left:win_.right]
                        )
                        if (abs(window_peak / waveform_peak) >
                                self.config.window_amplitude_ratio):
                            windows_by_amplitude.append(win_)
                        else:
                            logger.info(
                                "removing window due to global amplitude ratio: "
                                f"{ abs(window_peak / waveform_peak)} < "
                                f"{self.config.window_amplitude_ratio}")
                            continue
                    window = windows_by_amplitude
                # Check if amplitude windowing removed windows
                if window:
                    windows[comp] = window

                _nwin = len(window)
            except IndexError:
                _nwin = 0

                # Count windows and tell User
                num_windows += _nwin
                logger.info(f"{_nwin} window(s) for comp {comp}")

        return windows, num_windows

    def save_windows(self, data_type="MisfitWindows"):
        """
        Save the misfit windows that are calculated by Pyflex into a Dataset
        """
        logger.debug("saving misfit windows to PyASDF")
        for comp in self.windows.keys():
            for i, window in enumerate(self.windows[comp]):
                tag = (f"{self.config.model_number or 'Default'}/"
                       f"{self.st_obs[0].stats.network}_"
                       f"{self.st_obs[0].stats.station}_{comp}_{i}"
                       )

                # ASDF auxiliary_data subgroups don't play nice with nested
                # dictionaries, which the window parameters are. Format them
                # a bit simpler for saving into the dataset
                window_dict = create_window_dictionary(window)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    self.ds.add_auxiliary_data(data=np.array([True]),
                                               data_type=data_type,
                                               parameters=window_dict,
                                               path=tag
                                               )

    def measure(self, force=False):
        """
        Run Pyadjoint on Obs and Syn data for misfit windows or on full trace.

        Method for caluculating misfit set in Config, Pyadjoint expects
        standardized traces with the same spectral content, so this function
        will not run unless these flags are passed.

        Returns a dictionary of adjoint sources based on component.
        Saves resultant dictionary to a pyasdf dataset if given.

        :type force: bool
        :param force: ignore flag checks and run function, useful if e.g.
            external preprocessing is used that doesn't meet flag criteria
        """
        self._check()

        # Check that data has been filtered and standardized
        if not self._standardize_flag and not force:
            logger.warning("cannot measure misfit, traces not standardized")
            return
        elif not (self._obs_filter_flag and self._syn_filter_flag) \
                and not force:
            logger.warning(
                "cannot measure misfit, waveforms not filtered")
            return
        elif self._num_windows == 0 and not force:
            logger.warning("cannot measure misfit, no windows")
            return
        logger.info(
            f"running Pyadjoint w/ adj_src_type: {self.config.adj_src_type}")

        # Create list of windows needed for Pyadjoint
        adjoint_windows = self._format_windows()

        # Run Pyadjoint to retrieve adjoint source objects
        total_misfit, adjoint_sources = 0, {}
        for comp, adj_win in adjoint_windows.items():
            adj_src = pyadjoint.calculate_adjoint_source(
                adj_src_type=src_type(self.config.adj_src_type),
                config=self.config.pyadjoint_config,
                observed=self.st_obs.select(component=comp)[0],
                synthetic=self.st_syn.select(component=comp)[0],
                window=adj_win, plot=False
                )

            # Save adjoint sources in dictionary object. Sum total misfit
            adjoint_sources[comp] = adj_src
            logger.info(f"{adj_src.misfit:.3f} misfit for comp {comp}")
            total_misfit += adj_src.misfit

        # Save adjoint source internally for plotting
        self.adj_srcs = adjoint_sources

        # Save adjoint source to dataset
        if self.ds:
            self.save_adj_srcs()
        
        # Save total misfit, calculated a la Tape (2010) Eq. 6
        if self._num_windows:
            self._misfit = 0.5 * total_misfit/self._num_windows
        else:
            self._misfit = total_misfit

        # Let the User know the outcome of Pyadjoint
        logger.info(f"total misfit {self._misfit:.3f}")

    def _format_windows(self):
        """
        In `pyadjoint.calculate_adjoint_source`, the window needs to be a list
        of lists, with each list containing the [left_window, right_window];
        each window argument should be given in units of time (seconds)
        Note:
            This is not in the PyAdjoint docs

        :rtype: dict of list of lists
        :return: dictionary with key related to individual components,
            and corresponding to a list of lists containing window start and end
        """
        adjoint_windows = {}
        if self.windows is not None:
            for comp, window in self.windows.items():
                adjoint_windows[comp] = []
                # Prepare Pyflex window indices to give to Pyadjoint
                for win in window:
                    adj_win = [win.left * self.st_obs[0].stats.delta,
                               win.right * self.st_obs[0].stats.delta]
                    adjoint_windows[comp].append(adj_win)
        # If no windows given, calculate adjoint source on whole trace
        else:
            logger.debug("no windows given, adjoint sources will be "
                         "calculated on full trace")
            for comp in self.config.component_list:
                adjoint_windows[comp] = [[0, self.st_obs.select(
                                                component=comp)[0].stats.npts]]

        return adjoint_windows

    def save_adj_srcs(self):
        """
        Save adjoint sources to Pyasdf Dataset
        """
        for key, adj_src in self.adj_srcs.items():
            logger.debug(f"saving adjoint sources {key} to PyASDF")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # The tag hardcodes an X as the second channel index
                # to signify that these are synthetic waveforms
                # This is required by Specfem3D
                tag = "{mod}/{net}_{sta}_{ban}X{cmp}".format(
                    mod=self.config.model_number or "Default",
                    net=adj_src.network, sta=adj_src.station,
                    ban=channel_codes(self.st_syn[0].stats.delta),
                    cmp=adj_src.component[-1]
                )
                write_adj_src_to_asdf(adj_src=adj_src, ds=self.ds, tag=tag,
                                      time_offset=self._time_offset_sec)

    def plot(self, append_title='', length_sec=None,
             figsize=(11.69, 8.27), dpi=100, show=True, save=None,
             return_figure=False, **kwargs):
        """
        Waveform plots for all given components.
        If specific components are not given they are omitted from the plot.
        Plotting should be dynamic, i.e. if only 2 components are present in the
        streams, only two subplots should be generated in the figure.

        :type append_title: str
        :param append_title: any string to append the title of the plot
        :type length_sec: str or int or float or None
        :param length_sec: if "dynamic", dynamically determine length of wav 
            based on src rcv distance. If type float or type int, User defined
            waveform length in seconds. If None, default to full trace length
        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :type figsize: tuple of floats
        :param figsize: length and width of the figure
        :type dpi: int
        :param dpi: dots per inch of the figure
        :type return_figure: bool
        :param return_figure: option to get the figure back
        """
        # Precheck for waveform data
        self._check()
        if not self._standardize_flag:
            logger.warning("cannot plot, waveforms not standardized")
            return
        logger.info("plotting waveform")

        # Calculate the seismogram length based on given options
        # Dynamically find length based on src rcv distanc
        if isinstance(length_sec, str) and (length_sec == "dynamic"):
            length_sec = seismogram_length(
                    slow_wavespeed_km_s=1, binsize=50, minimum_length=100,
                    distance_km=gcd_and_baz(event=self.event,
                                            sta=self.inv[0][0])[0])
        # If not int or not float, then default to showing full trace
        elif not (isinstance(length_sec, int) or isinstance(length_sec, float)):
            length_sec = None

        # Call on window making function to produce waveform plots
        fig_window = window_maker(
            st_obs=self.st_obs, st_syn=self.st_syn, config=self.config,
            time_offset_sec=self._time_offset_sec, windows=self.windows,
            staltas=self.staltas, adj_srcs=self.adj_srcs, length_sec=length_sec,
            append_title=append_title, figsize=figsize, dpi=dpi, show=show,
            save=save, **kwargs
        )
        if return_figure:
            return fig_window

    def srcrcvmap(self, map_corners=None, stations=None, show_nz_faults=False,
                  annotate_names=False, color_by_network=False,
                  figsize=(8, 8.27), dpi=100, show=True, save=None, **kwargs):
        """
        Map plot showing a map of the given target region. All stations that
        show data availability (according to the station master list) are
        plotted as open markers. Event is plotted as a beachball if a moment
        tensor is given, station of interest highlighted, both are connected
        with a dashed line.
        Source receier information plotted in lower right hand corner of map.

        :type map_corners: dict
        :param map_corners: {lat_min, lat_max, lon_min, lon_max}
        :type stations: obspy.core.inventory.Inventory
        :param stations: background stations to plot on the map that will
        not interact with the source
        :type annotate_names: bool
        :param annotate_names: annotate station names
        :type color_by_network: bool
        :param color_by_network: color based on network name
        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :type show_nz_faults: bool
        :param show_nz_faults: plot active faults and hikurangi trench from
            internally saved coordinate files. takes extra time over simple plot
        :type figsize: tuple of floats
        :param figsize: length and width of the figure
        :type dpi: int
        :param dpi: dots per inch of the figure
        """
        self._check()
        logger.info("plotting map")
 
        # Warn user if no inventory is given
        if not isinstance(self.inv, obspy.Inventory):
            logger.warning("no inventory given, plotting blank map")

        if map_corners is None:
            map_corners = self.config.map_corners

        # Call external function to generate map
        manager_map(map_corners=map_corners, inv=self.inv, event=self.event,
                    stations=stations, show_nz_faults=show_nz_faults,
                    color_by_network=color_by_network,
                    annotate_names=annotate_names, show=show, figsize=figsize,
                    dpi=dpi, save=save, **kwargs
                    )
