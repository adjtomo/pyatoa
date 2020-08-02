#!/usr/bin/env python3
"""
A class to control workflow and temporarily store and manipulate data
"""
import warnings

import obspy
import pyflex
import pyadjoint
from os.path import basename
from obspy.signal.filter import envelope

from pyatoa import logger
from pyatoa.core.config import Config
from pyatoa.core.gatherer import Gatherer, GathererNoDataException
from pyatoa.utils.process import is_preprocessed
from pyatoa.utils.asdf.fetch import windows_from_dataset
from pyatoa.utils.window import reject_on_global_amplitude_ratio
from pyatoa.utils.srcrcv import gcd_and_baz
from pyatoa.utils.asdf.add import add_misfit_windows, add_adjoint_sources
from pyatoa.utils.process import (preproc, trim_streams, stf_convolve, zero_pad, 
                                  match_npts)

from pyatoa.visuals.manager_plotter import ManagerPlotter


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
        self.nwin = 0
        self.len_obs = 0
        self.len_syn = 0 
        self.misfit = 0
        self.half_dur = 0
        self.time_offset_sec = 0
        self.standardized = False 
        self.obs_filtered = False 
        self.syn_filtered = False

    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]    

    def __str__(self):
        str_ = ""
        for key, value in self.items():
            str_ += f"{key:>15}: {value}\n"
        return str_[:-1]


class Manager:
    """
    Pyatoas core workflow object.

    Manager is the central workflow control object. It calls on mid and
    low level classes to gather data, standardize and preprocess stream objects,
    generate misfit windows, and calculate adjoint sources. Has a variety of
    internal sanity checks to ensure that the workflow stays on the rails.
    """
    def __init__(self, config=None, ds=None, empty=True, station_code=None,
                 event=None, st_obs=None, st_syn=None, inv=None, windows=None,
                 staltas=None, adjsrcs=None, gcd=None, baz=None,
                 gatherer=None):
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
        :type adjsrcs: dict of pyadjoint.AdjointSource objects
        :param adjsrcs: adjoint source waveforms stored in dictionaries
        :type gcd: float
        :param gcd: great circle distance between source and receiver in km
        :type baz: float
        :param baz: Backazimuth between source and receiver in units of degrees
        :type gatherer: pyatoa.core.gatherer.Gatherer
        :param gatherer: A previously instantiated Gatherer class.
            Should not have to be passed in by User, but is used for reset()
        """
        # Main workflow requirements
        if config:
            self.config = config
        else:
            self.config = None
        self.ds = ds
        self.gatherer = gatherer
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
        self.gcd = gcd
        self.baz = baz
        self.windows = windows
        self.staltas = staltas or {}
        self.adjsrcs = adjsrcs
        self.rejwins = {}

        # Internal statistics to keep track of the workflow progress
        self.stats = ManagerStats()

        # If event ID set, launch gatherer, gather an event
        if not empty or event is not None:
            self.setup(event=event)
        else:
            # If 'empty' or no event, dont launch gatherer, event is None
            self.event = None

        # Run internal checks on data
        self._check()

    def __str__(self):
        """
        Print statement shows available data detailing workflow
        """
        self._check()
        return ("Manager Data\n"
                f"    dataset   [ds]:        {self.stats.dataset_id}\n"
                f"    quakeml   [event]:     {self.stats.event_id}\n"
                f"    station   [inv]:       {self.stats.inv_name}\n"
                f"    observed  [st_obs]:    {self.stats.len_obs}\n"
                f"    synthetic [st_syn]:    {self.stats.len_syn}\n"
                "Stats and Status\n"  
                f"    half_dur:              {self.stats.half_dur}\n"
                f"    time_offset_sec:       {self.stats.time_offset_sec}\n"
                f"    standardized:          {self.stats.standardized}\n"
                f"    obs_filtered:          {self.stats.obs_filtered}\n"
                f"    syn_filtered:          {self.stats.syn_filtered}\n"
                f"    nwin (windows):        {self.stats.nwin}\n"
                f"    misfit (adjsrcs):      {self.stats.misfit:.2E}\n"
                )

    def __repr__(self):
        return self.__str__()

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

    def _check(self):
        """
        (Re)check the stats of the workflow and data within the Manager.

        Rechecks conditions whenever called, incase something has gone awry
        mid-workflow. Stats should only be set by this function.
        """
        # Give dataset filename if available
        if self.stats.dataset_id is None and self.ds is not None:
            self.stats.dataset_id = basename(self.ds.filename)

        # Determine the resource identifier for the Event object
        if self.stats.event_id is None and self.event is not None:
            self.stats.event_id = self.event.resource_id

        # Get the network and station name from the Inventory object
        if self.stats.inv_name is None and self.inv is not None:
            self.stats.inv_name = ".".join([self.inv[0].code,
                                            self.inv[0][0].code])

        # Check if waveforms are Stream objects, and if preprocessed
        if self.st_obs is not None:
            self.stats.len_obs = len(self.st_obs)
            self.stats.obs_filtered = is_preprocessed(self.st_obs)

        if self.st_syn is not None:
            self.stats.len_syn = len(self.st_syn)
            self.stats.syn_filtered = is_preprocessed(self.st_syn)

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
        if self.stats.half_dur is None and self.event is not None:
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

    def setup(self, event=None, idx=0, append_focal_mechanism=True):
        """
        One-time setup of the Manager class.

        Assigns or instantiates required auxiliary classes necessary for the
        Manger to function, including the Config, Gatherer and Event attributes.

        :type event: obspy.core.event.Event or obspy.core.event.catalog.Catalog
        :param event: event or Catalog to use for the central event attribute
        :type idx: int
        :param idx: if `event` is a Catalog, idx allows user to specify
            which index in Catalog to choose event from. Default is 0.
        """
        if self.config is None:
            logger.info("no Config found, initiating default")
            self.config = Config()

        # Launch or reset the Gatherer
        if self.gatherer is None:
            self.gatherer = Gatherer(config=self.config, ds=self.ds)

        # Determine event information
        if event is not None:
            # User provided Event/Catalog object should be distributed
            if isinstance(event, obspy.core.event.catalog.Catalog):
                logger.info(f"event given as catalog, taking entry {idx}")
                event = event[idx]
            self.event = event
            self.gatherer.origintime = event.preferred_origin().time
        else:
            # No User provided event, gather based on event id
            if self.config.event_id is not None:
                self.event = self.gatherer.gather_event(
                    append_focal_mechanism=append_focal_mechanism
                    )

    def reset(self):
        """
        Restart workflow by deleting all collected data in the Manager, but
        retain dataset, event, config, and gatherer so a new station can be
        processed with the same configuration as the previous workflow.
        """
        self.__init__(ds=self.ds, event=self.event, config=self.config,
                      gatherer=self.gatherer)
        self._check()

    def write(self, write_to="ds"):
        """
        Write the data collected inside Manager to either a Pyasdf Dataset,
        or to individual files (not implemented).

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
            if self.event:
                try:
                    self.ds.add_quakeml(self.event)
                except ValueError:
                    logger.warning("Event already present, not added")
            if self.inv:
                try:
                    self.ds.add_stationxml(self.inv)
                except TypeError:
                    logger.warning("StationXML already present, not added")
            # PyASDF has its own warnings if waveform data already present
            if self.st_obs:
                self.ds.add_waveforms(waveform=self.st_obs,
                                      tag=self.config.observed_tag)
            if self.st_syn:
                self.ds.add_waveforms(waveform=self.st_syn,
                                      tag=self.config.synthetic_tag)
            if self.windows:
                self.save_windows()
            if self.adjsrcs:
                self.save_adjsrcs()
        else:
            raise NotImplementedError

    def load(self, station_code, path=None, ds=None, synthetic_tag=None,
             observed_tag=None):
        """
        Populate the manager using a previously populated ASDFDataSet.
        Useful for re-instantiating an existing workflow that has already 
        gathered data and saved it to an ASDFDataSet.

        :type station_code: str
        :param station_code: SEED conv. code, e.g. NZ.BFZ.10.HHZ
        :type path: str
        :param path: if no Config object is given during init, the User
            can specify the config path here to load data from the dataset.
            This skips the need to initiate a separate Config object.
        :type ds: None or pyasdf.ASDFDataSet
        :param ds: dataset can be given to load from, will not set the ds
        :type synthetic_tag: str
        :param synthetic_tag: waveform tag of the synthetic data in the dataset
            e.g. 'synthetic_m00s00'
        :type observed_tag: str
        :param observed_tag: waveform tag of the observed data in the dataset
            e.g. 'observed'
        """
        # Allows a ds to be provided outside the attribute
        if self.ds and ds is None:
            ds = self.ds
        else:
            raise TypeError("load requires a Dataset")

        # If no Config object in Manager, load from dataset 
        if self.config is None:
            if path is None:
                raise TypeError("load requires Config or 'path'")
            else:
                self.config = Config(ds=ds, path=path)
                logger.info(f"loading config from dataset {path}")

        assert len(station_code.split('.')) == 2, \
            "station_code must be in form 'NN.SSS'"

        # Reset and populate using the dataset
        self.__init__(config=self.config, ds=ds, empty=True)
        self.event = ds.events[0]
        net, sta = station_code.split('.')
        sta_tag = f"{net}.{sta}"
        if sta_tag in ds.waveforms.list():
            self.inv = ds.waveforms[sta_tag].StationXML
            self.st_syn = ds.waveforms[sta_tag][synthetic_tag or
                                                self.config.synthetic_tag]
            self.st_obs = ds.waveforms[sta_tag][observed_tag or
                                                self.config.observed_tag]
        else:
            logger.warning(f"no data for {sta_tag} found in dataset")

        self._check()
        return self

    def flow(self, fix_windows=False, preprocess_overwrite=None, iteration=None,
             step_count=None):
        """
        A convenience function to run the full workflow with a single command.
        Does not include gathering as that could be from gather() or load()

        Will raise ManagerError for controlled exceptions meaning workflow 
        cannot continue.

        :type fix_windows: bool
        :param fix_windows: do not pick new windows, but load windows from the
            given dataset
        :type preprocess_overwrite: function
        :param preprocess_overwrite: overwrite the core preprocess functionality
        """
        self.standardize()
        self.preprocess(overwrite=preprocess_overwrite)
        self.window(fix_windows=fix_windows, iteration=iteration,
                    step_count=step_count)
        self.measure()

    def gather(self, station_code, choice=None, query_fdsn=True):
        """
        Gather station dataless and waveform data using the Gatherer class.
        In order collect observed waveforms, dataless, and finally synthetics.

        If any part of gathering fails, raise ManagerError

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type choice: list
        :param choice: allows user to gather individual bits of data, rather
            than gathering all. Allowed: 'inv', 'st_obs', 'st_syn'
        :rtype: bool
        :return: status of the function, 1: successful / 0: failed
        """
        # Default to gathering all data.
        if choice is None:
            choice = ["st_obs", "inv", "st_syn"]
        try:
            if self.gatherer is None:
                # Instantiate Gatherer and retrieve Event information
                self.setup()
            net, sta, loc, cha = station_code.split(".")
            logger.info(f"gathering data for {station_code}")
            if "st_obs" in choice:
                # Ensure observed waveforms gathered first, as if this fails
                # then there is no point to gathering the rest
                self.st_obs = self.gatherer.gather_observed(station_code)
            if "inv" in choice:
                self.inv = self.gatherer.gather_station(station_code)
            if "st_syn" in choice:
                self.st_syn = self.gatherer.gather_synthetic(station_code)

            return self
        except GathererNoDataException as e:
            # Catch the Gatherer exception and redirect as ManagerError 
            # so that it can be caught by flow()
            raise ManagerError("data gatherer could not find some data") from e
        except Exception as e:
            # Gathering should be robust, but if something slips through, dont
            # let it kill a workflow, display and raise ManagerError
            logger.warning(e, exc_info=True)
            raise ManagerError("Uncontrolled error in data gathering") from e

    def standardize(self, force=False, standardize_to="syn"):
        """
        Standardize the observed and synthetic traces in place. 
        Ensures Streams have the same starttime, endtime, sampling rate, npts.

        :type force: bool
        :param force: allow the User to force the functino to run even if checks
            say that the two Streams are already standardized
        :type standardize_to: str
        :param standardize_to: allows User to set which Stream conforms to which
            by default the Observed traces should conform to the Synthetic ones
            because exports to Specfem should be controlled by the Synthetic
            sampling rate, npts, etc.
        :rtype: bool
        :return: status of the function, 1: successful / 0: failed
        """
        self._check()
        if not self.stats.len_obs or not self.stats.len_syn:
            raise ManagerError("cannot standardize, not enough waveform data")
        elif self.stats.standardized and not force:
            logger.info("data already standardized")
            return self
        logger.info("standardizing streams")

        # If observations starttime after synthetic, zero pad the front of obs
        dt_st = self.st_obs[0].stats.starttime - self.st_syn[0].stats.starttime
        if dt_st > 0:
            self.st_obs = zero_pad(self.st_obs, dt_st, before=True, after=False)

        # Match sampling rates
        if standardize_to == "syn":
            self.st_obs.resample(self.st_syn[0].stats.sampling_rate)
        else:
            self.st_syn.resample(self.st_obs[0].stats.sampling_rate)

        # Match start and endtimes
        self.st_obs, self.st_syn = trim_streams(
            st_a=self.st_obs, st_b=self.st_syn,
            force={"obs": "a", "syn": "b"}[standardize_to]
            )

        # Match the number of samples 
        self.st_obs, self.st_syn = match_npts(
            st_a=self.st_obs, st_b=self.st_syn,
            force={"obs": "a", "syn": "b"}[standardize_to]
            )

        # Determine if synthetics start before the origintime
        if self.event is not None:
            self.stats.time_offset_sec = (self.st_syn[0].stats.starttime -
                                          self.event.preferred_origin().time
                                          )
            logger.debug(f"time offset is {self.stats.time_offset_sec}s")
        else:
            self.stats.time_offset_sec = 0

        self.stats.standardized = True

        return self

    def preprocess(self, which="both", overwrite=None):
        """
        Preprocessing of observed and synthetic data in place.

        Remove response (observed), rotate, filter, convolve with source time
        functin (synthetic).

        Can be overwritten by a User defined function that takes the Manager
        as an input variable.

        :type which: str
        :param which: "obs", "syn" or "both" to choose which stream to process
            defaults to both
        :type overwrite: function
        :param overwrite: If a function is provided, it will overwrite the 
            standard preprocessing function. All arguments that are given
            to the standard preprocessing function will be passed as kwargs to
            the new function. This allows for customized preprocessing
        :rtype: bool
        :return: status of the function, 1: successful / 0: failed
        """
        if not self.inv and not self.config.synthetics_only:
            raise ManagerError("cannot preprocess, no inventory")
        if overwrite:
            assert(hasattr(overwrite, '__call__')), "overwrite must be function"
            preproc_fx = overwrite
        else:
            preproc_fx = preproc

        # If required, will rotate based on source receiver lat/lon values
        if self.config.rotate_to_rtz:
            if not self.inv:
                logger.warning("cannot rotate components, no inventory")
            else:
                self.gcd, self.baz = gcd_and_baz(event=self.event,
                                                 sta=self.inv[0][0])

        # Preprocess observation waveforms
        if self.st_obs is not None and not self.stats.obs_filtered and \
                which.lower() in ["obs", "both"]:
            logger.info("preprocessing observation data")
            self.st_obs = preproc_fx(self, choice="obs")
            if self.config.synthetics_only and self.stats.half_dur:
                # A synthetic-synthetic case means observed needs STF too
                self.st_obs = stf_convolve(st=self.st_obs,
                                           half_duration=self.stats.half_dur
                                           )

        # Preprocess synthetic waveforms
        if self.st_syn is not None and not self.stats.syn_filtered and \
                which.lower() in ["syn", "both"]:
            logger.info("preprocessing synthetic data")
            self.st_syn = preproc_fx(self, choice="syn")
            if self.stats.half_dur:
                # Convolve synthetics with a Gaussian source time function
                self.st_syn = stf_convolve(st=self.st_syn,
                                           half_duration=self.stats.half_dur
                                           )

        # Set stats
        self.stats.len_obs = len(self.st_obs)
        self.stats.len_syn = len(self.st_syn)
        self.stats.obs_filtered = True
        self.stats.syn_filtered = True

        return self

    def window(self, fix_windows=False, iteration=None, step_count=None,
               force=False, save=True):
        """
        Evaluate misfit windows using Pyflex. Save windows to ASDFDataSet.
        Allows previously defined windows to be retrieved from ASDFDataSet.

        Note:
            -Windows are stored as dictionaries of pyflex.Window objects.
            -All windows are saved into the ASDFDataSet, even if retrieved.
            -STA/LTA information is collected and stored internally.

        :type fix_windows: bool
        :param fix_windows: do not pick new windows, but load windows from the
            given dataset from 'iteration' and 'step_count'
        :type iteration: int or str
        :param iteration: if 'fix_windows' is True, look for windows in this
            iteration. If None, will check the latest iteration/step_count
            in the given dataset
        :type step_count: int or str
        :param step_count: if 'fix_windows' is True, look for windows in this
            step_count. If None, will check the latest iteration/step_count
            in the given dataset
        :type force: bool
        :param force: ignore flag checks and run function, useful if e.g.
            external preprocessing is used that doesn't meet flag criteria
        :type save: bool
        :param save: save the gathered windows to an ASDF Dataset
        """
        # Pre-check to see if data has already been standardized
        self._check()

        if not self.stats.standardized and not force:
            raise ManagerError("cannot window, waveforms not standardized")
        if fix_windows and not self.ds:
            logger.warning("cannot fix window, no dataset")
            fix_windows = False
        elif fix_windows and (iteration is None or step_count is None):
            # If no iteration/step_count values are given, automatically search
            # the previous step_count for windows in relation to the current
            # iteration/step_count
            iteration = self.config.iteration
            step_count = self.config.step_count
            return_previous_step = True
        else:
            # If fix windows and iteration/step_count are given, search the
            # dataset for windows under the current iteration/step_count
            return_previous_step = False

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

        # Find misfit windows, from a dataset or through window selection
        if fix_windows:
            self.retrieve_windows(iteration, step_count, return_previous_step)
        else:
            self.select_windows_plus()

        if save:
            self.save_windows()
        logger.info(f"{self.stats.nwin} window(s) total found")

        return self

    def retrieve_windows(self, iteration, step_count, return_previous_step):
        """
        Mid-level window selection function that retrieves windows from a 
        PyASDF Dataset, recalculates window criteria, and attaches window 
        information to Manager. No access to rejected window information.

        :type iteration: int or str
        :param iteration: retrieve windows from the given iteration
        :type step_count: int or str
        :param step_count: retrieve windows from the given step count
            in the given dataset
        :type return_previous_step: bool
        :param return_previous_step: if True: return windows from the previous
            step count in relation to the given iteration/step_count.
            if False: return windows from the given iteration/step_count
        """
        logger.info(f"retrieving windows from dataset")

        net, sta, _, _ = self.st_obs[0].get_id().split(".")
        # Function will return empty dictionary if no acceptable windows found
        windows = windows_from_dataset(ds=self.ds, net=net, sta=sta,
                                       iteration=iteration,
                                       step_count=step_count,
                                       return_previous_step=return_previous_step
                                       )

        # Recalculate window criteria for new values for cc, tshift, dlnA etc...
        logger.debug("recalculating window criteria")
        for comp, windows_ in windows.items():
            try:
                d = self.st_obs.select(component=comp)[0].data
                s = self.st_syn.select(component=comp)[0].data
                for w, win in enumerate(windows_):
                    # Post the old and new values to the logger for sanity check
                    logger.debug(f"{comp}{w}_old - "
                                 f"cc:{win.max_cc_value:.2f} / "
                                 f"dt:{win.cc_shift:.1f} / "
                                 f"dlnA:{win.dlnA:.2f}")
                    win._calc_criteria(d, s)
                    logger.debug(f"{comp}{w}_new - "
                                 f"cc:{win.max_cc_value:.2f} / "
                                 f"dt:{win.cc_shift:.1f} / "
                                 f"dlnA:{win.dlnA:.2f}")
            # IndexError thrown when trying to access an empty Stream
            except IndexError:
                continue

        self.windows = windows
        self.stats.nwin = sum(len(_) for _ in self.windows.values())

    def select_windows_plus(self):
        """
        Mid-level custom window selection function that calls Pyflex select 
        windows, but includes additional window suppression functionality.
        Includes custom Pyflex addition of outputting rejected windows, which
        will be used internally for plotting.
        """
        logger.info(f"running Pyflex w/ map: {self.config.pyflex_preset}")

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
                windows = ws.select_windows()

            # vvv Further window filtering can be applied here vvv
            # Suppress windows that contain low-amplitude signals
            if self.config.win_amp_ratio > 0:
                windows, ws.rejects["amplitude"] = \
                       reject_on_global_amplitude_ratio(
                                            data=obs.data, windows=windows,
                                            ratio=self.config.win_amp_ratio
                                            )
            # ^^^ Further window filtering can be applied here ^^^
            if windows:
                window_dict[comp] = windows
            if ws.rejects:
                reject_dict[comp] = ws.rejects

            # Count windows and tell User
            logger.info(f"{len(windows)} window(s) selected for comp {comp}")
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

        Note:
            Pyadjoint returns an unscaled misfit value for an entire set of
            windows. To return a "total misfit" value as defined by 
            Tape (2010) Eq. 6, the total summed misfit will need to be scaled by 
            the number of misfit windows chosen in Manager.window().

        :type force: bool
        :param force: ignore flag checks and run function, useful if e.g.
            external preprocessing is used that doesn't meet flag criteria
        :rtype: bool
        :return: status of the function, 1: successful / 0: failed
        """
        self._check()

        # Check that data has been filtered and standardized
        if not self.stats.standardized and not force:
            raise ManagerError("cannot measure misfit, not standardized")
        elif not (self.stats.obs_filtered and self.stats.syn_filtered) \
                and not force:
            raise ManagerError("cannot measure misfit, not filtered")
        elif not self.stats.nwin and not force:
            raise ManagerError("cannot measure misfit, no windows")
        logger.debug(f"running Pyadjoint w/ type: {self.config.adj_src_type}")

        # Create list of windows needed for Pyadjoint
        adjoint_windows = self._format_windows()

        # Run Pyadjoint to retrieve adjoint source objects
        total_misfit, adjoint_sources = 0, {}
        for comp, adj_win in adjoint_windows.items():
            try:
                adj_src = pyadjoint.calculate_adjoint_source(
                    adj_src_type=self.config.adj_src_type,
                    config=self.config.pyadjoint_config,
                    observed=self.st_obs.select(component=comp)[0],
                    synthetic=self.st_syn.select(component=comp)[0],
                    window=adj_win, plot=False
                    )
                # Save adjoint sources in dictionary object. Sum total misfit
                adjoint_sources[comp] = adj_src
                logger.info(f"{adj_src.misfit:.3f} misfit for comp {comp}")
                total_misfit += adj_src.misfit
            except IndexError:
                continue

        # Save adjoint source internally and to dataset
        self.adjsrcs = adjoint_sources
        self.save_adjsrcs()

        # Run check to get total misfit
        self._check()
        logger.info(f"total misfit {self.stats.misfit:.3f}")

        return self

    def save_windows(self):
        """
        Convenience function to save collected misfit windows into an 
        ASDFDataSet with some preliminary checks

        Auxiliary data tag is hardcoded as 'MisfitWindows'
        """
        if self.ds is None:
            logger.warning("Manager has no ASDFDataSet, cannot save windows")
        elif not self.windows:
            logger.warning("Manager has no windows to save")
        elif not self.config.save_to_ds:
            logger.warning("config parameter save_to_ds is set False, "
                           "will not save windows")
        else:
            logger.debug("saving misfit windows to ASDFDataSet")
            add_misfit_windows(self.windows, self.ds, path=self.config.aux_path)

    def save_adjsrcs(self):
        """
        Convenience function to save collected adjoint sources into an 
        ASDFDataSet with some preliminary checks

        Auxiliary data tag is hardcoded as 'AdjointSources'        
        """
        if self.ds is None:
            logger.warning("Manager has no ASDFDataSet, cannot save "
                           "adjoint sources")
        elif not self.adjsrcs:
            logger.warning("Manager has no adjoint sources to save")
        elif not self.config.save_to_ds:
            logger.warning("config parameter save_to_ds is set False, "
                           "will not save adjoint sources")
        else:
            logger.debug("saving adjoint sources to ASDFDataSet")
            add_adjoint_sources(adjsrcs=self.adjsrcs, ds=self.ds,
                                path=self.config.aux_path,
                                time_offset=self.stats.time_offset_sec)

    def _format_windows(self):
        """
        In `pyadjoint.calculate_adjoint_source`, the window needs to be a list 
        of lists, with each list containing the [left_window, right_window]; 
        each window argument should be given in units of time (seconds). 

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

    def plot(self, choice="both", save=None, show=True, corners=None, **kwargs):
        """
        Plot observed and synthetics waveforms, misfit windows, STA/LTA and
        adjoint sources for all available components. Append information
        about misfit, windows and window selection.

        For valid key word arguments see ManagerPlotter

        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :param corners: {lat_min, lat_max, lon_min, lon_max}
            corners to cut the map to, otherwise a global map is provided
        :type choice: str
        :param choice: choice for what to plot:
            'wav': plot waveform figure only
            'map': plot a source-receiver map only
            'both' (default): plot waveform and source-receiver map together
        """
        self._check()
        # Precheck for correct data to plot
        if choice in ["wav", "both"] and not self.stats.standardized:
            raise ManagerError("cannot plot, waveforms not standardized")

        if choice in ["map", "both"] and (self.inv is None or
                                          self.event is None):
            raise ManagerError("cannot plot map, no event and/or inv found")

        mp = ManagerPlotter(mgmt=self)
        if choice == "wav":
            mp.plot_wav(show=show, save=save, **kwargs)
        elif choice == "map":
            mp.plot_map(corners=corners, show=show, save=save,  **kwargs)
        elif choice == "both":
            mp.plot(corners=corners, show=show, save=save, **kwargs)


