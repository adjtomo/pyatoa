#!/usr/bin/env python
"""
Mid level data gathering class that calls on the internal fetcher and the
external getter to grab data
"""
import copy
import obspy

from pyatoa.utils.gathering.internal_fetcher import Fetcher
from pyatoa.utils.gathering.external_getter import Getter


class Gatherer():
    """
    A class used to fetch data via internal_fetcher and external_getter
    dependent on data availability.
    """
    def __init__(self, config, ds=None, startpad=20, endpad=200):
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.event = None
        self.fetcher = Fetcher(config=self.config,ds=self.ds,startpad=startpad,
                               endpad=endpad)
        self.getter = Getter(config=self.config, startpad=startpad,
                             endpad=endpad)

    def _gather_event(self):
        """
        get event information, check internally first
        """
        if self.ds is not None:
            try:
                self.event = self.fetcher.asdf_event_fetch()
            except IndexError:
                self.event = self.getter.event_get()
                ds.add_quakeml(self.event)
        else:
            self.event = self.getter.event_get()
        self.fetcher.origintime = self.event.origins[0].time
        self.getter.origintime = self.event.origins[0].time

    def _gather_station(self, station_code):
        """
        get station information, check internally first.
        """
        try:
            return self.fetcher.station_fetch(station_code)
        except FileNotFoundError:
            return self.getter.station_get(station_code)

    def _gather_waveforms(self, station_code):
        """
        get waveform as obspy streams, check internally first.
        """
        try:
            st_obs = self.fetcher.obs_waveform_fetch(station_code)
        except FileNotFoundError:
            st_obs = self.getter.waveform_get(station_code)
        st_syn = self.fetcher.syn_waveform_fetch(station_code)
        return st_obs, st_syn

    def gather_all(self, station_code):
        """
        convenience function gather event, station and waveform data in one go
        """
        if self.event is None:
            self._gather_event()
        inv = self._gather_station(station_code)
        st_obs, st_syn = self._gather_waveforms(station_code)
        return st_obs, st_syn, inv, self.event

    def gather_moment_tensor(self):
        """
        get moment tensor information
        :return:
        """
        moment_tensor = self.fetcher.geonet_moment_tensor_fetch()
        gcmt = self.gcmt_fetch()

    # def write_to_pyasdf(self,ds):
    #     """
    #     save station, event information and raw waveforms into pyasdf
    #     """
    #     if self.event is not None:
    #         ds.add_quakeml(self.event)
    #     if self.inv is not None:
    #         ds.add_stationxml(self.inv)
    #     if self.st_obs is not None:
    #
    #     if self.st_syn is not None:
    #





















