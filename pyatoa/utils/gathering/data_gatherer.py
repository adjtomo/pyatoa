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
from pyatoa.utils.gathering.grab_auxiliaries import grab_geonet_moment_tensor
from pyatoa.utils.operations.source_receiver import generate_focal_mechanism


class Gatherer:
    """
    A class used to fetch data via internal_fetcher and external_getter
    dependent on data availability. Preferentially searches internal pathways
    and will fall back onto external pathways if no data availability.
    """
    def __init__(self, config, ds=None):
        """
        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset for internal data searching and saving
        :type start_pad: int
        :param start_pad: padding in seconds before the origin time of an event
            for waveform fetching, to be fed into lower level functions.
        :type end_pad: int
        :param end_pad: padding in seconds after the origin time of an event
            for wavefomr fetching.
        """
        self.config = config
        self.ds = ds
        self.event = None
        self.fetcher = Fetcher(config=self.config, ds=self.ds,
                               start_pad=self.config.start_pad,
                               end_pad=self.config.end_pad)
        self.getter = Getter(config=self.config,
                             start_pad=self.config.start_pad,
                             end_pad=self.config.end_pad)

    def gather_event(self):
        """
        Get event information, check internally first. Keep the event as an
        attribute and check if the event has already be retrieved as this
        only needs to be done once per Pyatoa workflow.
        Set the event time precision equal to 2 points after the decimal to keep
        the origin time the same as that saved into the CMTSOLUTION file.
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
                self.gather_focal_mechanism()
                self.ds.add_quakeml(self.event)
                logger.info("event got from external, added to pyasdf dataset")
        else:
            self.event = self.getter.event_get()
            self.gather_focal_mechanism()
            logger.info("event got from external")

        self.event.preferred_origin().time.precision = 2
        self.fetcher.origintime = self.event.preferred_origin().time
        self.getter.origintime = self.event.preferred_origin().time
        return self.event

    def gather_focal_mechanism(self):
        """
        NOTE: hardcoded paths to a .csv file containing GeoNet moment tensor
        values which have hardcoded variable names.

        FDSN fetched events are devoid of a few bits of information that are
        useful for our applications, e.g. moment tensor, focal mechanisms.
        This function will perform the conversions and append the necessary
        information to the event located in the dataset.

        :return:
        """
        if isinstance(self.event, obspy.core.event.Event):
            geonet_mtlist = grab_geonet_moment_tensor(self.config.event_id)
            generate_focal_mechanism(mtlist=geonet_mtlist, event=self.event)
            logger.info("appending GeoNet moment tensor information to event")

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
            logger.info(
                "internal station information not found, searching ext.")
            inv = self.getter.station_get(station_code)
            if self.ds is not None:
                self.ds.add_stationxml(inv)
            return inv

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
            logger.info("internal observation data unavailable, searching ext.")
            st_obs = self.getter.waveform_get(station_code)
            if self.ds is not None:
                self.ds.add_waveforms(waveform=st_obs, tag='observed')
        # if self.config.raw_sampling_rate is not None:
        #     if self.config.raw_sampling_rate < st_obs[0].stats.sampling_rate:
        #         st_obs.resample(sampling_rate=self.config.raw_sampling_rate)
        #     else:
        #         warnings.warn(
        #           "Raw sampling rate {r} is greater than original {o}".format(
        #                 r=self.config.raw_sampling_rate,
        #                 o=st_obs[0].stats.sampling_rate), UserWarning
        #                 )

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
        st_syn = self.fetcher.syn_waveform_fetch(station_code)
        # if self.config.raw_sampling_rate is not None:
        #     if self.config.raw_sampling_rate < st_syn[0].stats.sampling_rate:
        #         st_syn.resample(sampling_rate=self.config.raw_sampling_rate)
        #     else:
        #         warnings.warn(
        #           "Raw sampling rate {r} is greater than original {o}".format(
        #                 r=self.config.raw_sampling_rate,
        #                 o=st_syn[0].stats.sampling_rate), UserWarning
        #                 )
        return st_syn


