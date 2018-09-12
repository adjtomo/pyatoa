#!/usr/bin/env python3
"""
Invoke Pyadjoint and Pyflex on observed and synthetic data to generate misfit
windows and adjoint sources
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
    an internal storage class that Processor can dump information into to keep
    clutter out of the Processor class
    """
    def __init__(self, station_code=None):
        self.station_code = station_code
        self.st_obs = None
        self.st_syn = None
        self.inv = None
        self.event = None
        self.moment_tensor = None
        self.windows = None
        self.staltas = None
        self.adj_srcs = None

        self.eventflag = False
        self.stobsflag = False
        self.stsynflag = False
        self.invflag = False
        self.obsprocessflag = False
        self.synprocessflag = False
        self.syntheticdataflag = False
        self.pyflexflag = False
        self.pyadjointflag = False

    def _checkflags(self):
        """
        sort through the crate and figure out what's in there
        :return:
        """
        self.stobsflag = isinstance(self.st_obs, obspy.Stream)
        if self.stobsflag:
            self.obsprocessflag = hasattr(self.st_obs[0].stats, "processing")
        self.stsynflag = isinstance(self.st_syn, obspy.Stream)
        if self.stsynflag:
            self.synprocessflag = hasattr(self.st_syn[0].stats, "processing")
        self.invflag = isinstance(self.inv, obspy.Inventory)
        self.eventflag = isinstance(self.event, obspy.core.event.Event)
        self.pyflexflag = isinstance(self.windows, dict)
        self.pyadjointflag = isinstance(self.adj_srcs, dict)


class Processor:
    """
    wrapper to contain all the adjoint source creation functionalities
    """
    def __init__(self, config, ds=None):
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.gatherer = None
        self.crate = Crate()

    def __str__(self):
        """
        print statement for processor class
        :return:
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
                ).format(event=self.crate.eventflag,
                         obsstream=self.crate.stobsflag,
                         synstream=self.crate.stsynflag,
                         inventory=self.crate.invflag,
                         obsproc=self.crate.obsprocessflag,
                         synproc=self.crate.synprocessflag,
                         synshift=self.crate.syntheticdataflag,
                         pyflex=self.crate.pyflexflag,
                         pyadjoint=self.crate.pyadjointflag
                         )

    def gather_data(self, station_code):
        """
        call on data_gatherer to populate data space, catch on exception,
        possibility to check flags to see where data gathering failed
        """
        try:
            self.crate.station_code = station_code
            if self.gatherer is None:
                logger.info("initiating gatherer")
                gatherer = Gatherer(config=self.config, ds=self.ds)
                self.gatherer = gatherer

            logger.info("gathering source, receiver and waveform data")
            if self.gatherer.event is not None:
                self.crate.event = self.gatherer.event
            else:
                self.crate.event = self.gatherer.gather_event()
            self.crate.inv = self.gatherer.gather_station(station_code)
            self.crate.st_obs = self.gatherer.gather_observed(station_code)
            self.crate.st_syn = self.gatherer.gather_synthetic(station_code)
        except Exception as e:
            print(e)
            return

    def preprocess(self):
        """
        preprocess observed and synthetic data
        :return:
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
        put synthetic data into the correct origin time.
        !!!TODO: determine how important half duration and time shift are
        :return:
        """
        self.crate.moment_tensor = grab_geonet_moment_tensor(
            self.config.event_id)
        half_duration = half_duration_from_m0(self.crate.moment_tensor["Mo"])
        self.crate.st_syn = stf_convolve(st=self.crate.st_syn,
                                         half_duration=half_duration,
                                         window="bartlett",
                                         time_shift=False
                                         )

    def run_pyflex(self):
        """
        call pyflex module to calculate best fitting misfit windows
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
        if self.crate.windows is None:
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
        f = window_maker(st_obs=self.crate.st_obs, st_syn=self.crate.st_syn,
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
        m = generate_map(config=self.config, event=self.crate.event,
                         inv=self.crate.inv, show_faults=True,
                         show=show, figsize=(11.69, 8.27), dpi=100
                         )





