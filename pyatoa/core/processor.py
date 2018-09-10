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
from pyatoa.utils.processing.preproc import preproc
from pyatoa.utils.configurations.external_configuration import \
    set_pyflex_configuration, set_pyadjoint_configuration

class Processer():
    """
    wrapper to contain all the adjoint source creation functionalities
    """
    def __init__(self,config,ds):
        self.cfg = copy.deepcopy(config)
        self.ds = ds
        self.st_obs = None
        self.st_syn = None
        self.inv = None
        self.event = None
        self.gatherer = None
        self.windows = None
        self.staltas = None
        self.adj_srcs = None

    def gather_data(self,station_code):
        """
        call on data_gatherer to populate data space
        """
        if self.gatherer is None:
            gatherer = Gatherer(cfg=self.cfg,ds=self.ds)
            self.gatherer = gatherer
        self.st_obs, self.st_syn, self.inv, self.event = \
                                        self.gatherer.gather_all(station_code)

    def preprocess(self):
        """
        preprocess observed and synthetic data
        :return:
        """
        self.st_obs = preproc(self.st_obs, inv=self.inv, resample=50,
                              pad_length_in_seconds=20, output="VEL",
                              filter=[self.cfg.min_period,self.cfg.max_period]
                              )
        self.st_syn = preproc(self.st_syn, resample=50,
                              pad_length_in_seconds=20, output="VEL",
                              filter=[self.cfg.min_period,self.cfg.max_period]
                              )

    def run_pyflex(self):
        """
        call pyflex module to calculate best fitting misfit windows
        """
        pf_config, pf_event, pf_station = set_pyflex_configuration(
            config=self.cfg, station=self.inv, event=self.event
            )
        empties = 0
        windows, staltas = {}, {}
        for COMP in self.cfg.component_list:
            window = pyflex.select_windows(
                observed=self.st_obs.select(component=COMP),
                synthetic=self.st_syn.select(component=COMP),
                config=pf_config, event=pf_event, station=pf_station,
                )
            stalta = pyflex.stalta.sta_lta(
                data=envelope(
                    self.st_syn.select(component=COMP)[0].data),
                dt=self.st_syn.select(component=COMP)[0].stats.delta,
                min_period=self.cfg.min_period
                )
            staltas[COMP] = stalta

            # TODO: add log statement to print out how many windows found
            logger.info("{0} window(s) found for component {1}".format(
                len(window), COMP)
                )
            if not window:
                empties += 1
                continue
            else:
                windows[COMP] = window
        if empties == len(self.cfg.component_list):
            # TODO: should this be an exception?
            warnings.warn("Empty windows", UserWarning)

        if self.ds is not None:
            for COMP in self.cfg.component_list.keys():
                for i, window in enumerate(windows[COMP]):
                    internalpath = "{net}_{sta}_{comp}_{num}_{mod}".format(
                        evid=self.cfg.component_list,
                        net=self.st_obs[0].stats.network,
                        sta=self.st_obs[0].stats.station,
                        comp=COMP, mod=self.cfg.model_number, num=i)
                    win_dict = create_window_dictionary(window)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        # auxiliary data requires a data object, give it a bool
                        self.ds.add_auxiliary_data(data=np.array([True]),
                                                   data_type="MisfitWindows",
                                                   path=internalpath,
                                                   parameters=win_dict
                                                   )

        self.windows = windows
        self.staltas = staltas

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
        pa_config = set_pyadjoint_configuration(config=self.cfg)

        delta = st[0].stats.delta
        cfg = choose_config("pyadjoint", PD)

        adjoint_sources = {}
        for key in self.windows:
            # collect all windows into a single list object
            adjoint_windows = []
            for win in self.windows[key]:
                adj_win = [win.left * delta, win.right * delta]
                adjoint_windows.append(adj_win)

            adj_src = pyadjoint.calculate_adjoint_source(
                adj_src_type=self.cfg.adj_src_type,
                observed=self.st_obs.select(component=key)[0],
                synthetic=self.st_syn.select(component=key)[0],
                config=pa_config, window=adjoint_windows,
                plot=False
                )
            adjoint_sources[key] = adj_src

            if self.ds:
                logger.info("Saving adj src [{}] to PyASDF".format(key))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    adj_src.write_to_asdf(PD["dataset"], time_offset=0)

        self.adj_srcs = adjoint_sources

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




