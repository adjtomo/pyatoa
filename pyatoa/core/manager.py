#!/usr/bin/env python3
"""
Main workflow components of Pyatoa.
Manager is the central workflow control object. It calls on mid and low level
classes to gather data, and then runs these through Pyflex for misfit window
identification, and then into Pyadjoint for misfit quantification. Config class
required to set the necessary parameters

Crate class is a simple data storage object which is easily emptied and filled
such that the manager can remain relatively high level and not get bogged down
by excess storage requirements. Crate also feeds flags to the manager to
signal which processes have already occurred in the workflow.

TODO: create a moment tensor object that can be attached to the event object
"""
import copy
import warnings

import obspy
import pyflex
import pyadjoint
import numpy as np
from obspy.signal.filter import envelope

from pyatoa import logger
from pyatoa.utils.gathering.data_gatherer import Gatherer
from pyatoa.utils.operations.source_receiver import gcd_and_baz
from pyatoa.utils.operations.formatting import create_window_dictionary, \
    write_adj_src_to_asdf
from pyatoa.utils.processing.preproc import preproc, trimstreams
from pyatoa.utils.processing.synpreproc import stf_convolve_gaussian

from pyatoa.utils.configurations.external_configurations import \
    set_pyflex_configuration, set_pyadjoint_configuration


class Crate:
    """
    An internal storage class for clutter-free rentention of data for individual
    stations. Simple flagging system for quick glances at workflow progress.

    Mid-level object that is called and manipulated by the Manager class.
    """
    def __init__(self, station_code=None):
        """
        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type st_obs: obspy.core.stream.Stream
        :param st_obs: Stream object containing waveforms of observations
        :type st_syn: obspy.core.stream.Stream
        :param st_syn: Stream object containing waveforms of observations
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory that should only contain the station of interest,
            it's relevant channels, and response information
        :type event: obspy.core.event.Event
        :param event: An event object containing relevant earthquake information
        :type windows: dict of pyflex.Window objects
        :param windows: misfit windows calculated by Pyflex, stored in a
            dictionary based on component naming
        :type adj_srcs: dict of pyadjoint.AdjointSource objects
        :param adj_srcs: adjoint source waveforms stored in dictionaries

            (i.e. first column in the .sem? file) are shifted compared to the
            CMTSOLUTION, in units of seconds.
        :type *flag: bool
        :param *_flag: if * is present in the Crate, or if * has been processed
        """
        self.station_code = station_code
        self.st_obs = None
        self.st_syn = None
        self.inv = None
        self.event = None
        self.windows = None
        self.staltas = None
        self.adj_srcs = None

        self.event_flag = False
        self.st_obs_flag = False
        self.st_syn_flag = False
        self.inv_flag = False
        self.obs_process_flag = False
        self.syn_process_flag = False
        self.syn_shift_flag = False
        self.pyflex_flag = False
        self.pyadjoint_flag = False

    def _checkflags(self):
        """
        Update flags based on what is available in the crate. The 3 in the
        stream process flags comes from the 2 steps taken before preprocessing,
        downsampling and trimming
        """
        if isinstance(self.st_obs, obspy.Stream):
            self.st_obs_flag = len(self.st_obs)
            self.obs_process_flag = (
                    hasattr(self.st_obs[0].stats, "processing") and
                    len(self.st_obs[0].stats.processing) >= 3
            )
        else:
            self.st_obs_flag, self.obs_process_flag = False, False
        if isinstance(self.st_syn, obspy.Stream):
            self.st_syn_flag = len(self.st_syn)
            self.syn_process_flag = (
                    hasattr(self.st_syn[0].stats, "processing") and
                    len(self.st_syn[0].stats.processing) >= 3
            )
        else:
            self.st_syn_flag, self.syn_process_flag = False, False
        if isinstance(self.inv, obspy.Inventory):
            self.inv_flag = "{net}.{sta}".format(net=self.inv[0].code,
                                                 sta=self.inv[0][0].code)
        else:
            self.inv_flag = False
        if isinstance(self.event, obspy.core.event.Event):
            self.event_flag = self.event.resource_id
        else:
            self.event_flag = False
        self.pyflex_flag = isinstance(self.windows, dict)
        self.pyadjoint_flag = isinstance(self.adj_srcs, dict)


class Manager:
    """
    Core object within Pyatoa.

    Workflow management function that internally calls on all other objects
    within the package in order to gather, process and analyze waveform data.
    """
    def __init__(self, config, ds=None, empty=False):
        """
        If no pyasdf dataset is given in the initiation of the Manager, all
        data fetching will happen via given pathways in the config file,
        or through external getting via FDSN pathways

        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: ASDF data set from which to read and write data
        :type gatherer: pyatoa.utils.gathering.data_gatherer.Gatherer
        :param gatherer: gathering function used to get and fetch data
        :type crate: pyatoa.core.Manager.Crate
        :param crate: Crate to hold all your information
        """
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.gatherer = None
        self.crate = Crate()
        if not empty:
            self.launch()

    def __str__(self):
        """
        Print statement shows available information inside the workflow.
        """
        self.crate._checkflags()
        return ("CRATE\n"
                "\tEvent:                     {event}\n"
                "\tInventory:                 {inventory}\n"
                "\tObserved Stream(s):        {obsstream}\n"
                "\tSynthetic Stream(s):       {synstream}\n"
                "MANAGER\n"
                "\tObs Data Preprocessed:       {obsproc}\n"
                "\tSyn Data Preprocessed:       {synproc}\n"
                "\tSynthetic Data Shifted:      {synshift}\n"
                "\tPyflex runned:               {pyflex}\n"
                "\tPyadjoint runned:            {pyadjoint}\n"
                ).format(event=self.crate.event_flag,
                         obsstream=self.crate.st_obs_flag,
                         synstream=self.crate.st_syn_flag,
                         inventory=self.crate.inv_flag,
                         obsproc=self.crate.obs_process_flag,
                         synproc=self.crate.syn_process_flag,
                         synshift=self.crate.syn_shift_flag,
                         pyflex=self.crate.pyflex_flag,
                         pyadjoint=self.crate.pyadjoint_flag
                         )

    def _launch_gatherer(self, reset=False):
        """
        Initiates gatherer class, should only need to be done once per asdf
        dataset. Option given to reset and reinstantiate gatherer object.

        :param reset: bool
        :return: if True, launch a new gatherer regardless of whether one exists
        """
        if (self.gatherer is None) or reset:
            logger.info("initiating gatherer")
            gatherer = Gatherer(config=self.config, ds=self.ds)
            self.gatherer = gatherer

    def reset(self):
        """
        To avoid user interaction with the Crate class.
        Convenience function to instantiate a new Crate, and hence start the
        workflow from the start without losing your event or gatherer.
        """
        self.crate = Crate()
        self.launch()

    @property
    def event(self):
        return self.gatherer.event

    @property
    def st(self):
        if isinstance(self.crate.st_syn, obspy.Stream) and \
                isinstance(self.crate.st_obs, obspy.Stream):
            return self.crate.st_syn + self.crate.st_obs
        elif isinstance(self.crate.st_syn, obspy.Stream) and \
                not isinstance(self.crate.st_obs, obspy.Stream):
            return self.crate.st_syn
        elif isinstance(self.crate.st_obs, obspy.Stream) and \
                not isinstance(self.crate.st_syn, obspy.Stream):
            return self.crate.st_obs
        else:
            return None

    @property
    def st_obs(self):
        return self.crate.st_obs

    @property
    def st_syn(self):
        return self.crate.st_syn

    @property
    def inv(self):
        return self.crate.inv

    @property
    def windows(self):
        return self.crate.windows

    @property
    def adj_srcs(self):
        return self.crate.adj_srcs

    def launch(self):
        """
        Initiate the prerequisite parts of the Manager class. Populate with
        an obspy event object
        """
        self._launch_gatherer()
        if self.gatherer.event is not None:
            self.crate.event = self.gatherer.event
        else:
            logger.info("gathering event information")
            self.crate.event = self.gatherer.gather_event()

    def gather_data(self, station_code):
        """
        Launch a gatherer object and gather event, station and waveform
        information given a station code. Fills the crate based on information
        most likely to be available (we expect an event to be available more
        often than waveform data).
        Catches general exceptions along the way, stops gathering if errors.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        """
        try:
            self.crate.station_code = station_code
            logger.info("GATHERING {station} for {event}".format(
                station=station_code, event=self.config.event_id)
            )
            logger.info("gathering station information")
            self.crate.inv = self.gatherer.gather_station(station_code)
            logger.info("gathering observation waveforms")
            self.crate.st_obs = self.gatherer.gather_observed(station_code)
            logger.info("gathering synthetic waveforms")
            self.crate.st_syn = self.gatherer.gather_synthetic(station_code)
        except Exception as e:
            print(e)
            return

    def preprocess(self):
        """
        Preprocess observed and synthetic data in place on waveforms in crate.
        """
        if not (isinstance(self.crate.st_obs, obspy.Stream) and
                isinstance(self.crate.st_syn, obspy.Stream)
                ):
            warnings.warn("cannot preprocess, missing waveform data",
                          UserWarning)
            return

        logger.info("preprocessing observation data")
        if self.config.rotate_to_rtz:
            _, baz = gcd_and_baz(self.crate.event, self.crate.inv)
        else:
            baz = None
        # adjoint sources require the same sampling_rate as the synthetics
        sampling_rate = self.crate.st_syn[0].stats.sampling_rate
        self.crate.st_obs = preproc(self.crate.st_obs, inv=self.crate.inv,
                                    resample=sampling_rate,
                                    pad_length_in_seconds=20, back_azimuth=baz,
                                    output=self.config.unit_output,
                                    filterbounds=[self.config.min_period,
                                                  self.config.max_period],
                                    corners=4
                                    )
        logger.info("preprocessing synthetic data")
        self.crate.st_syn = preproc(self.crate.st_syn, resample=None,
                                    pad_length_in_seconds=20,
                                    output=self.config.unit_output,
                                    back_azimuth=baz, corners=4,
                                    filterbounds=[self.config.min_period,
                                                  self.config.max_period],

                                    )
        self.crate.st_obs, self.crate.st_syn = trimstreams(
            st_a=self.crate.st_obs, st_b=self.crate.st_syn, force="b")
        # vv essentially the first timestamp in the .sem? file from specfem
        self.crate.time_offset = (self.crate.st_syn[0].stats.starttime -
                                  self.crate.event.preferred_origin().time
                                  )
        try:
            # convolve synthetic data with a gaussian source-time-function
            half_duration = (self.crate.event.focal_mechanisms[0].
                             moment_tensor.source_time_function.duration) / 2

            self.crate.st_syn = stf_convolve_gaussian(
                st=self.crate.st_syn, half_duration=half_duration,
                time_shift=False
            )
            self.crate.syn_shift_flag = True  # TODO: make this flag smarter
        except AttributeError:
            print("half duration value not found in event")

    def run_pyflex(self):
        """
        Call Pyflex to calculate best fitting misfit windows given observation
        and synthetic data in the crate. Return dictionaries of window objects,
        as well as STA/LTA traces, to the crate. If a pyasdf dataset is present,
        save misfit windows in as auxiliary data.
        If no misfit windows are found for a given station, throw a warning
        because pyadjoint won't run.
        Pyflex configuration is given by the config as a list of values with
        the following descriptions:

        i  Standard Tuning Parameters:
        0: water level for STA/LTA (short term average/long term average)
        1: time lag acceptance level
        2: amplitude ratio acceptance level (dlna)
        3: normalized cross correlation acceptance level
        i  Fine Tuning Parameters
        4: c_0 = for rejection of internal minima
        5: c_1 = for rejection of short windows
        6: c_2 = for rejection of un-prominent windows
        7: c_3a = for rejection of multiple distinct arrivals
        8: c_3b = for rejection of multiple distinct arrivals
        9: c_4a = for curtailing windows w/ emergent starts and/or codas
        10:c_4b = for curtailing windows w/ emergent starts and/or codas
        """
        if not (isinstance(self.crate.st_obs, obspy.Stream) and
                isinstance(self.crate.st_syn, obspy.Stream)
                ):
            warnings.warn("cannot run Pyflex, no waveform data")
            return

        pf_config, pf_event, pf_station = set_pyflex_configuration(
            config=self.config, inv=self.crate.inv, event=self.crate.event
            )
        empties = 0
        windows, staltas = {}, {}
        for COMP in self.config.component_list:
            window = pyflex.select_windows(
                observed=self.crate.st_obs.select(component=COMP),
                synthetic=self.crate.st_syn.select(component=COMP),
                config=pf_config, event=pf_event, station=pf_station,
                )
            stalta = pyflex.stalta.sta_lta(
                data=envelope(
                    self.crate.st_syn.select(component=COMP)[0].data),
                dt=self.crate.st_syn.select(component=COMP)[0].stats.delta,
                min_period=self.config.min_period
                )
            staltas[COMP] = stalta
            logger.info("{0} window(s) found for component {1}".format(
                len(window), COMP)
                )
            if not window:
                empties += 1
                continue
            else:
                windows[COMP] = window
        # TODO: should this be an exception?
        if empties == len(self.config.component_list):
            warnings.warn("Empty windows", UserWarning)

        self.crate.windows = windows
        self.crate.staltas = staltas

        if self.ds is not None:
            logger.info("Saving misfit windows to PyASDF")
            for COMP in windows.keys():
                for i, window in enumerate(windows[COMP]):
                    tag = "{mod}/{net}_{sta}_{cmp}_{num}".format(
                        net=self.crate.st_obs[0].stats.network,
                        sta=self.crate.st_obs[0].stats.station,
                        cmp=COMP, mod=self.config.model_number, num=i)
                    wind_dict = create_window_dictionary(window)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        # auxiliary needs data, give it a bool; auxis love bools
                        self.ds.add_auxiliary_data(data=np.array([True]),
                                                   data_type="MisfitWindows",
                                                   path=tag,
                                                   parameters=wind_dict
                                                   )

    def run_pyadjoint(self):
        """
        Run pyadjoint on observation and synthetic data given misfit windows
        calculated by pyflex. Method for caluculating misfit set in config,
        pyadjoint config set in external configurations. Returns a dictionary
        of adjoint sources based on component. Saves resultant dictionary into
        the crate, as well as to a pyasdf dataset if given.

        NOTE: This is not in the PyAdjoint docs, but in
        pyadjoint.calculate_adjoint_source, the window needs to be a list of
        lists, with each list containing the [left_window,right_window];
        each window argument should be given in units of time (seconds)
        """
        if (self.crate.windows is None) or (isinstance(self.crate.windows, dict)
                                            and not len(self.crate.windows)
                                            ):
            warnings.warn("cannot run Pyadjoint, no Pyflex outputs",
                          UserWarning)
            return

        logger.info("running pyAdjoint for type {} ".format(
            self.config.adj_src_type)
            )
        pa_config = set_pyadjoint_configuration(config=self.config)

        adjoint_sources = {}
        for key in self.crate.windows:
            adjoint_windows = []
            for win in self.crate.windows[key]:
                adj_win = [win.left * self.crate.st_obs[0].stats.delta,
                           win.right * self.crate.st_obs[0].stats.delta]
                adjoint_windows.append(adj_win)

            adj_src = pyadjoint.calculate_adjoint_source(
                adj_src_type=self.config.adj_src_type,
                observed=self.crate.st_obs.select(component=key)[0],
                synthetic=self.crate.st_syn.select(component=key)[0],
                config=pa_config, window=adjoint_windows,
                plot=False
                )
            adjoint_sources[key] = adj_src
            if self.ds is not None:
                logger.info("saving adjoint sources {} to PyASDF".format(key))
                with warnings.catch_warnings():
                    tag = "{mod}/{net}_{sta}_BX{cmp}".format(
                        mod=self.config.model_number, net=adj_src.network,
                        sta=adj_src.station, cmp=adj_src.component[-1]
                        )
                    warnings.simplefilter("ignore")
                    write_adj_src_to_asdf(adj_src, self.ds, tag,
                                          time_offset=self.crate.time_offset)
        self.crate.adj_srcs = adjoint_sources

    def plot_wav(self, **kwargs):
        """
        Waveform plots for all given components of the crate.
        If specific components are not given (e.g. adjoint source waveform),
        they are omitted from the final plot. Plotting should be dynamic, i.e.
        if only 2 components are present in the streams, only two subplots
        should be generated in the figure.
        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :type figsize: tuple of floats
        :param figsize: length and width of the figure
        :type dpi: int
        :param dpi: dots per inch of the figure
        """
        if not (isinstance(self.crate.st_obs, obspy.Stream) and
                isinstance(self.crate.st_syn, obspy.Stream)
                ):
            warnings.warn("cannot plot waveforms, no waveform data",
                          UserWarning)
            return

        from pyatoa.visuals.plot_waveforms import window_maker
        show = kwargs.get("show", True)
        save = kwargs.get("save", None)
        figsize = kwargs.get("figsize", (11.69, 8.27))
        dpi = kwargs.get("dpi", 100)

        window_maker(st_obs=self.crate.st_obs, st_syn=self.crate.st_syn,
                     windows=self.crate.windows, staltas=self.crate.staltas,
                     adj_srcs=self.crate.adj_srcs,
                     time_offset=self.crate.time_offset,
                     stalta_wl=self.config.pyflex_config[0],
                     unit_output=self.config.unit_output,
                     config=self.config, figsize=figsize, dpi=dpi, show=show,
                     save=save
                     )

    def plot_map(self, **kwargs):
        """
        Map plot showing a map of the given target region. All stations that
        show data availability (according to the station master list) are
        plotted as open markers. Event is plotted as a beachball if a moment
        tensor is given, station of interest highlighted, both are connected
        with a dashed line.
        Source receier information plotted in lower right hand corner of map.
        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :type show_faults: bool
        :param show_faults: plot active faults and hikurangi trench from
            internally saved coordinate files. takes extra time over simple plot
        :type figsize: tuple of floats
        :param figsize: length and width of the figure
        :type dpi: int
        :param dpi: dots per inch of the figure
        """
        if not isinstance(self.crate.inv, obspy.Inventory):
            warnings.warn("no inventory given, plotting blank map", UserWarning)
        from pyatoa.visuals.plot_map import generate_map
        show = kwargs.get("show", True)
        save = kwargs.get("save", None)
        show_faults = kwargs.get("show_faults", False)
        annotate_names = kwargs.get("annotate_names", False)

        figsize = kwargs.get("figsize", (8, 8.27))
        dpi = kwargs.get("dpi", 100)

        generate_map(config=self.config, event=self.crate.event,
                     inv=self.crate.inv, show_faults=show_faults,
                     annotate_names=annotate_names,
                     show=show, figsize=figsize, dpi=dpi, save=save
                     )






