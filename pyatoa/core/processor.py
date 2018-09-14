#!/usr/bin/env python3
"""
Invoke Pyadjoint and Pyflex on observed and synthetic data to generate misfit
windows and adjoint sources
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
from pyatoa.utils.operations.formatting import create_window_dictionary
from pyatoa.utils.processing.preproc import preproc, trimstreams
from pyatoa.utils.processing.synpreproc import stf_convolve, \
    half_duration_from_m0
from pyatoa.utils.configurations.external_configurations import \
    set_pyflex_configuration, set_pyadjoint_configuration
from pyatoa.utils.gathering.grab_auxiliaries import grab_geonet_moment_tensor


class Crate:
    """
    An internal storage class for clutter-free rentention of data for individual
    stations. Simple flagging system for quick glances at workflow progress.

    Mid-level object that is called and manipulated by the Processor class.
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
        :type *flag: bool
        :param *_flag: if * is present in the Crate, or if * has been processed
        """
        self.station_code = station_code
        self.st_obs = None
        self.st_syn = None
        self.inv = None
        self.event = None
        self.moment_tensor = None
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
        Update flags based on what is available in the crate.
        """
        self.st_obs_flag = isinstance(self.st_obs, obspy.Stream)
        if self.st_obs_flag:
            self.obsprocess_flag = hasattr(self.st_obs[0].stats, "processing")
        self.st_syn_flag = isinstance(self.st_syn, obspy.Stream)
        if self.st_syn_flag:
            self.syn_process_flag = hasattr(self.st_syn[0].stats, "processing")
        self.inv_flag = isinstance(self.inv, obspy.Inventory)
        self.event_flag = isinstance(self.event, obspy.core.event.Event)
        self.pyflex_flag = isinstance(self.windows, dict)
        self.pyadjoint_flag = isinstance(self.adj_srcs, dict)


class Processor:
    """
    Core object within Pyatoa.

    Workflow management function that internally calls on all other objects
    within the package in order to gather, process and analyze waveform data.

    """
    def __init__(self, config, ds=None):
        """
        If no pyasdf dataset is given in the initiation of the processor, all
        data fetching will happen via given pathways in the config file,
        or through external getting via FDSN pathways

        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: ASDF data set from which to read and write data
        :type gatherer: pyatoa.utils.gathering.data_gatherer.Gatherer
        :param gatherer: gathering function used to get and fetch data
        :type crate: pyatoa.core.processor.Crate
        :param crate: Crate to hold all your information
        """
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.gatherer = None
        self.crate = Crate()

    def __str__(self):
        """
        Print statement shows available information inside the workflow.
        """
        self.crate._checkflags()
        return ("Processor class\n"
                "\tData Gathering:\n"
                "\t\tEvent:                     {event}\n"
                "\t\tInventory:                 {inventory}\n"
                "\t\tObserved Stream:           {obsstream}\n"
                "\t\tSynthetic Stream:          {synstream}\n"
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

    def _reset(self):
        """
        To avoid user interaction with the Crate class.
        Convenience function to instantiate a new Crate, and hence start the
        workflow from the start without losing your event or gatherer.
        """
        self.crate = Crate()

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
            self._launch_gatherer()
            self.crate.station_code = station_code
            logger.info("gathering event information")
            if self.gatherer.event is not None:
                self.crate.event = self.gatherer.event
            else:
                self.crate.event = self.gatherer.gather_event()
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
        logger.info("preprocessing observation data")
        if self.config.rotate_to_rtz:
            _, baz = gcd_and_baz(self.crate.event, self.crate.inv)
        self.crate.st_obs = preproc(self.crate.st_obs, inv=self.crate.inv,
                                    resample=5, pad_length_in_seconds=20,
                                    output="VEL", back_azimuth=baz,
                                    filterbounds=[self.config.min_period,
                                                  self.config.max_period]
                                    )
        logger.info("preprocessing synthetic data")
        self.crate.st_syn = preproc(self.crate.st_syn, resample=5,
                                    pad_length_in_seconds=20, output="VEL",
                                    back_azimuth=baz,
                                    filterbounds=[self.config.min_period,
                                                  self.config.max_period]
                                    )
        self.crate.st_obs, self.crate.st_syn = trimstreams(self.crate.st_obs,
                                                           self.crate.st_syn)

    def shift_synthetic(self):
        """
        TODO: determine how important half duration and time shift are
        TODO: move the guts of this function into synpreproc.py (?)

        Put synthetic data into the correct origin time given a GeoNet moment
        tensor list retrieved from the GeoNet moment tensor CSV file.

        Convolve the synthetic waveforms in place with a shape function
        (default is a bartlett, which is basically a triangle function)
        """
        self.crate.moment_tensor = grab_geonet_moment_tensor(
            self.config.event_id)
        half_duration = half_duration_from_m0(self.crate.moment_tensor["Mo"])
        self.crate.st_syn = stf_convolve(st=self.crate.st_syn,
                                         half_duration=half_duration,
                                         window="bartlett",
                                         time_shift=False
                                         )
        self.crate.syn_shift_flag = True  # TODO: make this flag change smarter

    def run_pyflex(self):
        """
        Call Pyflex to calculate best fitting misfit windows given observation
        and synthetic data in the crate. Return dictionaries of window objects,
        as well as STA/LTA traces, to the crate. If a pyasdf dataset is 
        """
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
            for COMP in windows.keys():
                for i, window in enumerate(windows[COMP]):
                    internalpath = "{mod}/{net}_{sta}_{comp}_{num}".format(
                        evid=self.config.component_list,
                        net=self.crate.st_obs[0].stats.network,
                        sta=self.crate.st_obs[0].stats.station,
                        comp=COMP, mod=self.config.model_number, num=i)
                    wind_dict = create_window_dictionary(window)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        # auxiliary needs data, give it a bool; auxis love bools
                        self.ds.add_auxiliary_data(data=np.array([True]),
                                                   data_type="MisfitWindows",
                                                   path=internalpath,
                                                   parameters=wind_dict
                                                   )

    def run_pyadjoint(self):
        """
        !!! Not in the PyAdjoint docs:
        in pyadjoint.calculate_adjoint_source: window needs to be a list of
        lists, with each list containing the [left_window,right_window] where
        each window argument is given in seconds
        """
        if (self.crate.windows is None) or (isinstance(self.crate.windows, dict)
                                            and not len(self.crate.windows)
                                            ):
            warnings.warn("windows must be collected before pyadjoint can run")
            return

        logger.info("Running pyAdjoint for type {} ".format(
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
            if self.ds:
                logger.info("Saving adjoint source {} to PyASDF".format(key))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    adj_src.write_to_asdf(self.ds, time_offset=0)
        self.crate.adj_srcs = adjoint_sources

    def populate(self, station_code):
        """
        convenience function to run through the pyatoa process
        :return:
        """
        self.gather_data(station_code)
        self.preprocess()
        self.shift_synthetic()
        self.run_pyflex()
        self.run_pyadjoint()

    def plot_wav(self, show=True, *args, **kwargs):
        """

        :param args:
        :param kwargs:
        :return:
        """
        from pyatoa.visuals.plot_waveforms import window_maker
        window_maker(st_obs=self.crate.st_obs, st_syn=self.crate.st_syn,
                     windows=self.crate.windows, staltas=self.crate.staltas,
                     adj_srcs=self.crate.adj_srcs,
                     stalta_wl=self.config.pyflex_config[0],
                     unit_output=self.config.unit_output,
                     config=self.config, figsize=(11.69, 8.27), dpi=100,
                     show=show
                     )

    def plot_map(self, show=True, save=False, *args, **kwargs):
        """

        :param args:
        :param kwargs:
        :return:
        """
        from pyatoa.visuals.plot_map import generate_map
        generate_map(config=self.config, event=self.crate.event,
                     inv=self.crate.inv, show_faults=True,
                     show=show, figsize=(11.69, 8.27), dpi=100
                     )





