#!/usr/bin/env python3
"""
Invoke Pyadjoint and Pyflex on observed and synthetic data to generate misfit
windows and adjoint sources
"""
import copy
import warnings

import pyflex
import numpy as np
from obspy.signal.filter import envelope

from pyatoa import logger
from pyatoa.utils.gathering.data_gatherer import Gatherer
from pyatoa.utils.operations.formatting import create_window_dictionary
from pyatoa.utils.processing.preproc import preproc, trimstreams
from pyatoa.utils.processing.synpreproc import stf_convolve, half_duration_from_m0
from pyatoa.utils.configurations.external_configurations import \
    set_pyflex_configuration, set_pyadjoint_configuration
from pyatoa.utils.gathering.grab_auxiliaries import grab_geonet_moment_tensor

class Crate:
    """
    an internal storage class that Processor can dump information into
    """
    def __init__(self):
        self.station_code = None
        self.st_obs = None
        self.st_syn = None
        self.inv = None
        self.event = None
        self.moment_tensor = None
        self.windows = None
        self.staltas = None
        self.adj_srcs = None

    def set_station(self, station_code):
        """
        reset all parameters in the crate for a new station
        :param station_code:
        :return:
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


class Processor:
    """
    wrapper to contain all the adjoint source creation functionalities
    """
    def __init__(self, config, ds=None):
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.gatherer = None
        self.crate = Crate()

        self.eventflag = False
        self.stobsflag = False
        self.stsynflag = False
        self.invflag = False
        self.dataprocessflag = False
        self.syntheticdataflag = False
        self.pyflexflag = False
        self.pyadjointflag = False

    def __str__(self):
        """
        print statement for processor class
        :return:
        """
        template = ("Processor class\n"
                    "\tData Gathering:\n"
                    "\t\tEvent:                     {ev}\n"
                    "\t\tInventory:                 {iv}\n"
                    "\t\tStation Code:              {sc}\n"
                    "\t\tObserved Stream:           {os}\n"
                    "\t\tSynthetic Stream:          {ss}\n"
                    "\tData Preprocessed:           {pr}\n"
                    "\tSynthetic Data Shifted:      {sy}\n"
                    "\tPyflex run:                  {pf}\n"
                    "\tPyadjoint run:               {pa}\n"
                    )
        return template.format(ev=self.eventflag, sc=self.crate.station_code,
                               os=self.stobsflag, ss=self.stsynflag,
                               iv=self.invflag, pr=self.dataprocessflag,
                               sy=self.syntheticdataflag, pf=self.pyflexflag,
                               pa=self.pyadjointflag
                               )

    def gather_data(self, station_code):
        """
        call on data_gatherer to populate data space, catch on exception,
        possibility to check flags to see where data gathering failed
        """
        try:
            self.crate.set_station(station_code)
            if self.gatherer is None:
                logger.info("initiating gatherer")
                gatherer = Gatherer(config=self.config, ds=self.ds)
                self.gatherer = gatherer

            logger.info("gathering source, receiver and waveform data")
            if self.gatherer.event is not None:
                self.crate.event = self.gatherer.event
            else:
                self.crate.event = self.gatherer.gather_event()
            self.eventflag = True
            self.crate.inv = self.gatherer.gather_station(station_code)
            self.invflag = True
            self.crate.st_obs = self.gatherer.gather_observed(station_code)
            self.stobsflag = True
            self.crate.st_syn = self.gatherer.gather_synthetic(station_code)
            self.stsynflag = True
        except Exception:
            return

    def preprocess(self):
        """
        preprocess observed and synthetic data
        :return:
        """
        self.crate.st_obs = preproc(self.crate.st_obs, inv=self.crate.inv,
                                    resample=5, pad_length_in_seconds=20,
                                    output="VEL",
                                    filter=[self.config.min_period,
                                            self.config.max_period]
                                    )
        self.crate.st_syn = preproc(self.crate.st_syn, resample=5,
                                    pad_length_in_seconds=20, output="VEL",
                                    filter=[self.config.min_period,
                                            self.config.max_period]
                                    )
        self.crate.st_obs, self.crate.st_syn = trimstreams(self.crate.st_obs,
                                                           self.crate.st_syn)
        self.dataprocessflag = True

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
        self.syntheticdataflag = True

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
        if empties == len(self.config.component_list):
            # TODO: should this be an exception?
            warnings.warn("Empty windows", UserWarning)

        if self.ds is not None:
            for COMP in self.config.component_list.keys():
                for i, window in enumerate(windows[COMP]):
                    internalpath = "{net}_{sta}_{comp}_{num}_{mod}".format(
                        evid=self.config.component_list,
                        net=self.crate.st_obs[0].stats.network,
                        sta=self.crate.st_obs[0].stats.station,
                        comp=COMP, mod=self.config.model_number, num=i)
                    win_dict = create_window_dictionary(window)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        # auxiliary data requires a data object, give it a bool;
                        # auxi's love bools
                        self.ds.add_auxiliary_data(data=np.array([True]),
                                                   data_type="MisfitWindows",
                                                   path=internalpath,
                                                   parameters=win_dict
                                                   )

        self.crate.windows = windows
        self.crate.staltas = staltas

    def run_pyadjoint(self):
        """
        !!! Not in the PyAdjoint docs:
        in pyadjoint.calculate_adjoint_source: window needs to be a list of
        lists, with each list containing the [left_window,right_window] where
        each window argument is given in seconds
        """
        logger.info("Running pyAdjoint for type [{}] ".format(
            self.config.adj_src_type)
            )
        pa_config = set_pyadjoint_configuration(config=self.config)
        delta = self.crate.st_obs[0].stats.delta

        adjoint_sources = {}
        for key in self.crate.windows:
            # collect all windows into a single list object
            adjoint_windows = []
            for win in self.crate.windows[key]:
                adj_win = [win.left * delta, win.right * delta]
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
                logger.info("Saving adj src [{}] to PyASDF".format(key))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    adj_src.write_to_asdf(PD["dataset"], time_offset=0)

        self.crate.adj_srcs = adjoint_sources

    def process(self, station_code):
        """
        run through the pyatoa process for a single station
        :return:
        """
        self.crate.set_station(station_code)
        self.gather_data()






    def plot_waveforms(self, *args, **kwargs):
        """

        :param args:
        :param kwargs:
        :return:
        """

    def plot_map(self, *args, **kwargs):
        """

        :param args:
        :param kwargs:
        :return:
        """




