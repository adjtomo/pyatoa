#!/usr/bin/env python3
"""
A class to control workflow and temporarily store and manipulate data
"""
import os
import obspy
import pyflex
import warnings
import pyadjoint
from obspy.signal.filter import envelope
from pyatoa import logger
from pyatoa.core.config import Config
from pyatoa.core.gatherer import Gatherer, GathererNoDataException
from pyatoa.utils.form import channel_code
from pyatoa.utils.process import is_preprocessed
from pyatoa.utils.asdf.load import load_windows, load_adjsrcs
from pyatoa.utils.window import reject_on_global_amplitude_ratio
from pyatoa.utils.srcrcv import gcd_and_baz
from pyatoa.utils.asdf.add import add_misfit_windows, add_adjoint_sources
from pyatoa.utils.process import (default_process, trim_streams, zero_pad,
                                  match_npts)

from pyatoa.visuals.mgmt_plot import ManagerPlotter


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
    def __init__(self, config=None, ds=None, event=None, st_obs=None,
                 st_syn=None, inv=None, windows=None, staltas=None,
                 adjsrcs=None, gcd=None, baz=None, gatherer=None):
        """
        Initiate the Manager class with or without pre-defined attributes.

        .. note::
            If `ds` is not given in data can only be gathered via the
            config.paths attribute or using the ObsPy client service.
            Data will also not be saved.

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
        :type gatherer: pyatoa.core.gatherer.Gatherer
        :param gatherer: A previously instantiated Gatherer class.
            Should not have to be passed in by User, but is used for reset()
        """
        self.ds = ds
        self.inv = inv

        # Instantiate a Config object
        if config is not None:
            self.config = config
        else:
            logger.info("no config provided, initiating default")
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

        # Instantiate a Gatherer object and pass along info
        if gatherer is None:
            self.gatherer = Gatherer(config=self.config, ds=self.ds,
                                     origintime=origintime)
        else:
            self.gatherer = gatherer

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
        retain dataset, event, config, and gatherer so a new station can be
        processed with the same configuration as the previous workflow.
        """
        self.__init__(ds=self.ds, event=self.event, config=self.config,
                      gatherer=self.gatherer)

    def write(self, write_to="ds"):
        """
        Write the data collected inside Manager to either a Pyasdf Dataset,
        or to individual files (not implemented).

        :type write_to: str
        :param write_to: choice to write data to, if "ds" writes to a
            pyasdf.asdf_data_set.ASDFDataSet

            * write_to == "ds":
                If gather is skipped but data should still be saved into an
                ASDFDataSet for data storage, this function will
                fill that dataset in the same fashion as the Gatherer class
            * write_to == "/path/to/output":
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

    def write_adjsrcs(self, path="./", write_blanks=True):
        """
        Write internally stored adjoint source traces into SPECFEM3D defined
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
        from copy import deepcopy

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

    def load(self, code, path=None, ds=None, synthetic_tag=None,
             observed_tag=None, config=True, windows=False,
             adjsrcs=False):
        """
        Populate the manager using a previously populated ASDFDataSet.
        Useful for re-instantiating an existing workflow that has already 
        gathered data and saved it to an ASDFDataSet.

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
                logger.warning(f"No Config object in dataset for path {path}")

        assert(self.config is not None), "Config object required for load"
        assert len(code.split('.')) == 2, "'code' must be in form 'NN.SSS'"
        if windows or adjsrcs:
            assert(path is not None), "'path' required to load auxiliary data"
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
                self.windows = load_windows(ds, net, sta, iter_, step, False)
            if adjsrcs:
                self.adjsrcs = load_adjsrcs(ds, net, sta, iter_, step)
        else:
            logger.warning(f"no data for {sta_tag} found in dataset")

        self.check()
        return self

    def flow(self, **kwargs):
        """
        A convenience function to run the full workflow with a single command.
        Does not include gathering. Takes kwargs related to all underlying
        functions.

        .. code:: python

            mgmt = Manager()
            mgmt.flow() == mgmt.standardize().preprocess().window().measure()

        :raises ManagerError: for any controlled exceptions
        """
        force = kwargs.get("force", False)
        standardize_to = kwargs.get("standardize_to", "syn")
        fix_windows = kwargs.get("fix_windows", False)
        iteration = kwargs.get("iteration", None)
        step_count = kwargs.get("step_count", None)
        overwrite = kwargs.get("overwrite", None)
        which = kwargs.get("which", "both")
        save = kwargs.get("save", True)

        self.standardize(standardize_to=standardize_to, force=force)
        self.preprocess(overwrite=overwrite, which=which, **kwargs)
        self.window(fix_windows=fix_windows, iteration=iteration,
                    step_count=step_count, force=force, save=save)
        self.measure(force=force, save=save)

    def gather(self, code=None, choice=None, event_id=None, **kwargs):
        """
        Gather station dataless and waveform data using the Gatherer class.
        In order collect observed waveforms, dataless, and finally synthetics.

        For valid kwargs see methods in :doc:`core.gatherer`

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type choice: list
        :param choice: allows user to gather individual bits of data, rather
            than gathering all. Allowed: 'inv', 'st_obs', 'st_syn'
        :raises ManagerError: if any part of the gathering fails.

        Keyword Arguments
        ::
            bool try_fm:
                Try to retrieve and append focal mechanism information to the
                Event object.
            str prefix:
                Prefix for event id when searching for event information,
                can be used to search ordered files e.g., CMTSOLUTION_001
            str suffix:
                Suffix for event id when searching for event information
            str station_level:
                The level of the station metadata if retrieved using the ObsPy
                Client. Defaults to 'response'
            str resp_dir_template:
                Directory structure template to search for response files.
                By default follows the SEED convention:
                'path/to/RESPONSE/{sta}.{net}/'
            str resp_fid_template:
                Response file naming template to search for station dataless.
                By default, follows the SEED convention
                'RESP.{net}.{sta}.{loc}.{cha}'
            str obs_dir_template:
                directory structure to search for observation data. Follows the
                SEED convention: 'path/to/obs_data/{year}/{net}/{sta}/{cha}'
            str obs_fid_template:
                File naming template to search for observation data. Follows the
                SEED convention: '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
            str syn_cfgpath:
                Config.cfgpaths key to search for synthetic data. Defaults to
                'synthetics', but for the may need to be set to 'waveforms' in
                certain use-cases, e.g. synthetics-synthetic inversions.
            str syn_unit:
                Optional argument to specify the letter used to identify the
                units of the synthetic data: For Specfem3D: ["d", "v", "a", "?"]
                'd' for displacement, 'v' for velocity,  'a' for acceleration.
                Wildcards okay. Defaults to '?'
            str syn_dir_template:
                Directory structure template to search for synthetic waveforms.
                Defaults to empty string
            str syn_fid_template:
                The naming template of synthetic waveforms defaults to:
                "{net}.{sta}.*{cmp}.sem{syn_unit}"
        """
        # Default to gathering all data
        if choice is None:
            choice = ["event", "inv", "st_obs", "st_syn"]
        try:
            # Attempt to gather event information before waveforms/metadata
            if "event" in choice and self.event is None:
                if event_id is None:
                    event_id = self.config.event_id
                self.event = self.gatherer.gather_event(event_id, **kwargs)
            if code is not None:
                logger.info(f"gathering data for {code}")
                if "st_obs" in choice:
                    # Ensure observed waveforms gathered before synthetics and
                    # metadata. If this fails, no point to gathering the rest
                    self.st_obs = self.gatherer.gather_observed(code, **kwargs)
                if "inv" in choice:
                    self.inv = self.gatherer.gather_station(code, **kwargs)
                if "st_syn" in choice:
                    self.st_syn = self.gatherer.gather_synthetic(code, **kwargs)

            return self
        except GathererNoDataException as e:
            # Catch the Gatherer exception and redirect as ManagerError 
            # so that it can be caught by flow()
            raise ManagerError("Data Gatherer could not find some data") from e
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
        :param force: allow the User to force the function to run even if checks
            say that the two Streams are already standardized
        :type standardize_to: str
        :param standardize_to: allows User to set which Stream conforms to which
            by default the Observed traces should conform to the Synthetic ones
            because exports to Specfem should be controlled by the Synthetic
            sampling rate, npts, etc.
        """
        self.check()
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

    def preprocess(self, which="both", overwrite=None, **kwargs):
        """
        Preprocess observed and synthetic waveforms in place.
        Default preprocessing tasks: Remove response (observed), rotate, filter,
        convolve with source time function (synthetic).

        .. note::
            Default preprocessing can be overwritten using a
            user-defined function that takes Manager and choice as inputs
            and outputs an ObsPy Stream object.

        .. note::
            Documented kwargs only apply to default preprocessing.

        :type which: str
        :param which: "obs", "syn" or "both" to choose which stream to process
            defaults to both
        :type overwrite: function
        :param overwrite: If a function is provided, it will overwrite the 
            standard preprocessing function. All arguments that are given
            to the standard preprocessing function will be passed as kwargs to
            the new function. This allows for customized preprocessing

        Keyword Arguments
        ::
            int water_level:
                water level for response removal
            float taper_percentage:
                amount to taper ends of waveform
            bool remove_response:
                remove instrument response using the Manager's inventory object.
                Defaults to True
            bool apply_filter:
                filter the waveforms using the Config's min_period and
                max_period parameters. Defaults to True
            bool convolve_with_stf:
                Convolve synthetic data with a Gaussian source time function if
                a half duration is provided.
        """
        if not self.inv and not self.config.synthetics_only:
            raise ManagerError("cannot preprocess, no inventory")
        if overwrite:
            assert(hasattr(overwrite, '__call__')), "overwrite must be function"
            preproc_fx = overwrite
        else:
            preproc_fx = default_process

        # If required, will rotate based on source receiver lat/lon values
        if self.config.rotate_to_rtz:
            if not self.inv:
                logger.warning("cannot rotate components, no inventory")
            else:
                self.gcd, self.baz = gcd_and_baz(event=self.event,
                                                 sta=self.inv[0][0])

        # Preprocess observation waveforms
        if self.st_obs is not None and not self.stats.obs_processed and \
                which.lower() in ["obs", "both"]:
            logger.info("preprocessing observation data")
            self.st_obs = preproc_fx(self, choice="obs", **kwargs)
            self.stats.obs_processed = True

        # Preprocess synthetic waveforms
        if self.st_syn is not None and not self.stats.syn_processed and \
                which.lower() in ["syn", "both"]:
            logger.info("preprocessing synthetic data")
            self.st_syn = preproc_fx(self, choice="syn", **kwargs)
            self.stats.syn_processed = True

        # Set stats
        self.stats.len_obs = len(self.st_obs)
        self.stats.len_syn = len(self.st_syn)

        return self

    def window(self, fix_windows=False, iteration=None, step_count=None,
               force=False, save=True):
        """
        Evaluate misfit windows using Pyflex. Save windows to ASDFDataSet.
        Allows previously defined windows to be retrieved from ASDFDataSet.

        .. note::
            * Windows are stored as dictionaries of pyflex.Window objects.
            * All windows are saved into the ASDFDataSet, even if retrieved.
            * STA/LTA information is collected and stored internally.

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
        self.check()

        if self.config.pyflex_preset is None:
            logger.info("pyflex preset is set to 'None', will not window")
            return

        if not self.stats.standardized and not force:
            raise ManagerError("cannot window, waveforms not standardized")

        # Determine how to treat fixed windows
        if fix_windows and not self.ds:
            logger.warning("cannot fix window, no dataset")
            fix_windows = False
        elif fix_windows and (iteration is None or step_count is None):
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
            self.retrieve_windows(iteration, step_count, return_previous)
        else:
            self.select_windows_plus()

        if save:
            self.save_windows()
        logger.info(f"{self.stats.nwin} window(s) total found")

        return self

    def retrieve_windows(self, iteration, step_count, return_previous):
        """
        Mid-level window selection function that retrieves windows from a 
        PyASDF Dataset, recalculates window criteria, and attaches window 
        information to Manager. No access to rejected window information.

        :type iteration: int or str
        :param iteration: retrieve windows from the given iteration
        :type step_count: int or str
        :param step_count: retrieve windows from the given step count
            in the given dataset
        :type return_previous: bool
        :param return_previous: if True: return windows from the previous
            step count in relation to the given iteration/step_count.
            if False: return windows from the given iteration/step_count
        """
        logger.info(f"retrieving windows from dataset")

        net, sta, _, _ = self.st_obs[0].get_id().split(".")
        # Function will return empty dictionary if no acceptable windows found
        windows = load_windows(ds=self.ds, net=net, sta=sta,
                               iteration=iteration, step_count=step_count,
                               return_previous=return_previous
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
            logger.info(f"{len(windows)} window(s) selected for comp {comp}")
            nwin += len(windows)

        self.windows = window_dict
        self.rejwins = reject_dict
        self.stats.nwin = nwin

    def measure(self, force=False, save=True):
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
        :type save: bool
        :param save: save adjoint sources to ASDFDataSet
        """
        self.check()

        if self.config.adj_src_type is None:
            logger.info("adjoint source type is 'None', will not measure")
            return

        # Check that data has been filtered and standardized
        if not self.stats.standardized and not force:
            raise ManagerError("cannot measure misfit, not standardized")
        elif not (self.stats.obs_processed and self.stats.syn_processed) \
                and not force:
            raise ManagerError("cannot measure misfit, not filtered")
        elif self.stats.nwin == 0 and not force:
            raise ManagerError("cannot measure misfit, no windows recovered")
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

                # Re-format component name to reflect SPECFEM convention
                adj_src.component = f"{channel_code(adj_src.dt)}X{comp}"

                # Save adjoint sources in dictionary object. Sum total misfit
                adjoint_sources[comp] = adj_src
                logger.info(f"{adj_src.misfit:.3f} misfit for comp {comp}")
                total_misfit += adj_src.misfit
            except IndexError:
                continue

        # Save adjoint source internally and to dataset
        self.adjsrcs = adjoint_sources
        if save:
            self.save_adjsrcs()

        # Run check to get total misfit
        self.check()
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
        .. note::
            In `pyadjoint.calculate_adjoint_source`, the window needs to be a
            list of lists, with each list containing the
            [left_window, right_window]; each window argument should be given in
            units of time (seconds). This is not in the PyAdjoint docs.

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

    def plot(self, choice="both", save=None, show=True, corners=None, **kwargs):
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
        """
        self.check()
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


