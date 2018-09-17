#!/usr/bin/env python
"""
Mid level data gathering class that calls on the internal fetcher and the
external getter to grab event, station and waveform information.

Directly called by the processor class in the main Pyatoa workflow
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
    dependent on data availability. Preferentially searches internal pathways
    and will fall back onto external pathways if no data availability.
    """
    def __init__(self, config, ds=None, startpad=20, endpad=200):
        """
        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset for internal data searching and saving
        :type startpad: int
        :param startpad: padding in seconds before the origin time of an event
            for waveform fetching, to be fed into lower level functions.
        :type endpad: int
        :param endpad: padding in seconds after the origin time of an event
            for wavefomr fetching.
        """
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.event = None
        self.fetcher = Fetcher(config=self.config, ds=self.ds,
                               startpad=startpad, endpad=endpad)
        self.getter = Getter(config=self.config, startpad=startpad,
                             endpad=endpad)

    def gather_event(self):
        """
        Get event information, check internally first. Keep the event as an
        attribute and check if the event has already be retrieved as this
        only needs to be done once per Pyatoa workflow.
        :rtype: obspy.core.event.Event
        :return: event retrieved either via internal or external methods
        """
        if self.event is not None:
            return self.event
        elif self.ds is not None:
            try:
                self.event = self.fetcher._asdf_event_fetch()
            except IndexError:
                self.event = self.getter.event_get()
                self.ds.add_quakeml(self.event)
                logger.info("event got from external, added to pyasdf dataset")
        else:
            self.event = self.getter.event_get()
            logger.info("event got from external")
        self.fetcher.origintime = self.event.origins[0].time
        self.getter.origintime = self.event.origins[0].time
        return self.event

    def gather_moment_tensor(self):
        """
        Get moment tensor information via compeltely hardcoded pathways

        Deprecated and will be deleted.

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

    def gather_station(self, station_code):
        """
        Get station information, check internally in pyasdf dataset. If empty
        pyasdf dataset, save station if found

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        try:
            return self.fetcher.station_fetch(station_code)
        except FileNotFoundError:
            inv = self.getter.station_get(station_code)
            if self.ds is not None:
                self.ds.add_stationxml(inv)
            return None

    def gather_observed(self, station_code):
        """
        get waveform as obspy streams, check internally first.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        try:
            st_obs = self.fetcher.obs_waveform_fetch(station_code)
        except FileNotFoundError:
            st_obs = self.getter.waveform_get(station_code)
            if self.ds is not None:
                self.ds.add_waveforms(waveform=st_obs, tag='observed')
        return st_obs

    def gather_synthetic(self, station_code):
        """
        gather synthetic data

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        return self.fetcher.syn_waveform_fetch(station_code)
