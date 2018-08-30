#!/usr/bin/env python3
"""
Invoke Pyadjoint and Pyflex on observed and synthetic data to generate misfit
windows and adjoint sources
"""
import copy

from .config import Config
from .data_gather import Gatherer

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
        gatherer = Gatherer(cfg=self.cfg,ds=None)
        gatherer.gather_all(station_code)
        self.gatherer = gatherer
        if isinstance(ds,pyasdf.ASDFDataSet):
