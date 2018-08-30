#!/usr/bin/env python
"""
Data fetching class
"""
import copy
import obspy

from pyatoa.utils.data_gather.internal_fetcher import Fetcher
from pyatoa.utils.data_gather.external_getter import Getter


class Gatherer():
    """
    A class used to fetch data via internal_fetcher and external_getter
    dependent on data availability.
    """
    def __init__(self, config, ds=None, startpad=20, endpad=200):
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.startpad = startpad
        self.endpad = endpad
        self.fetcher = Fetcher(config=self.config,ds=self.ds)
        self.getter = Getter(config=self.config)
        self.event = event
        self.origintime = None
        self.st_obs = None
        self.st_syn = None
        self.inv = None

    def _gather_event(self):
        """
        get event information, check internally first
        """
        if self.ds is not None:
            try:
                self.event = self.fetcher.asdf_event_fetch()
            except IndexError:
                self.event = self.getter.event_get()
        else:
            self.event = self.getter.event_get()
        self.origintime = self.event.origins[0].time

    def _gather_station(self, station_code):
        """
        get station information, check internally first.
        """
        if self.ds is not None:
            try:
                self.inv = self.fetcher.asdf_station_fetch(station_code)
            except KeyError:
                self.inv = self.getter.station_get(station_code)
        else:
            self.inv = self.getter.station_get(station_code)

    def _gather_waveforms(self, station_code):
        """
        get waveform as obspy streams, check internally first.
        """
        if self.ds is not None:
            try:
                self.st_obs = self.fetcher.obs_waveform_fetch(station_code)
            except FetcherException:
                self.st_obs = self.getter.waveform_get(station_code)
        else:
            self.st_obs = self.getter.waveform_get(station_code)
        self.st_syn = self.fetcher.syn_waveform_fetch(station_code)

    def gather_all(self, station_code):
        """
        gather event, station and waveform data in one go
        """
        if not isinstance(self.event, obspy.core.event.event.Event):
            self._gather_event()
        self._gather_station(station_code)
        self._gather_waveforms(station_code)

    def write_to_pyasdf(self,ds):
        """
        save station, event information and raw waveforms into pyasdf
        """
        ds.add_quakeml(self.event)
        ds.add_stationxml(self.inv)
        ds.add_waveforms(waveform=self.st_obs,tag='observed',
            event_id=self.event)
        ds.add_waveforms(waveform=self.st_syn,tag='synthetic_{}'.format(
            self.config.model_number),event_id=self.event
            )





















