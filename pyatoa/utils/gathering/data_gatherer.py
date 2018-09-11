#!/usr/bin/env python
"""
Mid level data gathering class that calls on the internal fetcher and the
external getter to grab data
"""
import copy
import warnings

import obspy

from pyatoa import logger
from pyatoa.utils.gathering.internal_fetcher import Fetcher
from pyatoa.utils.gathering.external_getter import Getter
from pyatoa.utils.gathering.grab_auxiliaries import grab_geonet_moment_tensor,\
    grab_gcmt_moment_tensor, timeshift_halfduration


class Gatherer:
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

    def gather_event(self):
        """
        get event information, check internally first
        """
        if self.event is not None:
            return
        if self.ds is not None:
            try:
                self.event = self.fetcher.asdf_event_fetch()
                logger.info("event fetched from pyasdf dataset")
            except IndexError:
                self.event = self.getter.event_get()
                ds.add_quakeml(self.event)
                logger.info("event got from external, added to pyasdf dataset")
        else:
            self.event = self.getter.event_get()
            logger.info("event got from external")
        self.fetcher.origintime = self.event.origins[0].time
        self.getter.origintime = self.event.origins[0].time
        return self.event

    def gather_moment_tensor(self):
        """
        get moment tensor information
        TODO: this probably needs to be changed once we figure out how to get
        the empircal scaling law for moment magnitude and half duration
        :return:
        """
        warnings.warn("Old method using GCMT which is not always reliable,"
                      "instead use gather_mt", DeprecationWarning)
        try:
            geonet_moment_tensor = grab_geonet_moment_tensor(
                self.config.event_id)
        except AttributeError:
                warnings.warn("Geonet moment tensor doesn't exist", UserWarning)
                return None, None, None
        try:
            gcmt_moment_tensor = grab_gcmt_moment_tensor(
                obspy.UTCDateTime(geonet_moment_tensor["Date"]),
                geonet_moment_tensor["Mw"]
                )
        except FileNotFoundError:
            warnings.warn("GCMT moment tensor doesn't exist", UserWarning)
            return geonet_moment_tensor, None, None
        time_shift, half_duration = timeshift_halfduration(gcmt_moment_tensor,
                                                           geonet_moment_tensor
                                                           )
        return geonet_moment_tensor, time_shift, half_duration

    def gather_momtens(self):
        """
        gather geonet moment tensor from internal pathway, calculate the half
        duration using empirical scaling
        :return:
        """

    def gather_station(self, station_code):
        """
        get station information, check internally first.
        """
        try:
            return self.fetcher.station_fetch(station_code)
        except FileNotFoundError:
            return self.getter.station_get(station_code)

    def gather_observed(self, station_code):
        """
        get waveform as obspy streams, check internally first.
        """
        try:
            st_obs = self.fetcher.obs_waveform_fetch(station_code)
        except FileNotFoundError:
            st_obs = self.getter.waveform_get(station_code)
        return st_obs

    def gather_synthetic(self, station_code):
        """
        gather synthetic data
        :param station_code:
        :return:
        """
        return self.fetcher.syn_waveform_fetch(station_code)

    def gather_waveforms(self, station_code):
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
            self.gather_event()
        inv = self.gather_station(station_code)
        st_obs, st_syn = self.gather_waveforms(station_code)
        return st_obs, st_syn, inv, self.event


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





















