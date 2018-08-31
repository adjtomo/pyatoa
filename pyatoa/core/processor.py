#!/usr/bin/env python3
"""
Invoke Pyadjoint and Pyflex on observed and synthetic data to generate misfit
windows and adjoint sources
"""
import copy

from pyatoa.core.config import Config
from pyatoa.utisl.gathering.data_gatherer import Gatherer

class Processer():
    """
    wrapper to contain all the adjoint source creation functionalities
    """
    def __init__(self,config,ds,event=None):
        self.cfg = copy.deepcopy(config)
        self.ds = ds
        self.event = None
        self.gatherer = None

    def gather_data(self,station_code):
        """
        call on data_gatherer to populate data space
        """
        gatherer = Gatherer(cfg=self.cfg,ds=self.ds)
        gatherer.gather_all(station_code)

    def run_pyflex
