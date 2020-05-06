#!/usr/bin/env python
"""
Mid level data gathering class that calls on the internal fetcher and the
external getter to grab event, station and waveform information.

Directly called by the Manager class in the main Pyatoa workflow, shouldn't
need to be called by the User.
"""
import obspy

from pyatoa import logger
from obspy.clients.fdsn import Client
from pyatoa.utils.gather.internal_fetcher import InternalFetcher
from pyatoa.utils.gather.external_getter import (ExternalGetter,
                                                 get_gcmt_moment_tensor
                                                 )
from plugins.new_zealand.grab_auxiliaries import grab_geonet_moment_tensor
from pyatoa.utils.srcrcv import generate_focal_mechanism


class Gatherer(InternalFetcher, ExternalGetter):
    """
    A class used to fetch data internally and externally. Inherets internal
    and external data gathering from separate classes.

    Fetch: used to denote searching local filesystems for data (internal)
    Get: used to denote querying FDSN webservices via ObsPy (external), naming
         *convention from the obspy.fdsn.client.Client function names
    """
    def __init__(self, config, ds=None, event_id=None, origintime=None):
        """
        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset for internal data searching and saving
        """
        self.ds = ds
        self.event = None
        self.config = config
        self.Client = Client(self.config.client)

        # Event can be externally gathered by event_id or origin_time
        self.event_id = event_id
        if self.event_id is None and self.config.event_id is not None:
            self.event_id = self.config.event_id
        self.origintime = origintime

    def set_event(self, event):
        """
        If event is provided by User, propogate origintime to fetcher and
        getter classes so they can look for the correct data
        :type: obspy.core.event.Event
        :param: event supplied by User
        """
        self.event = event
        self.origintime = self.event.preferred_origin().time

    def gather_event(self, append_focal_mechanism=True):
        """
        Get event information, check internally first. Keep the event as an
        attribute and check if the event has already be retrieved as this
        only needs to be done once per Pyatoa workflow.
        Set the event time precision equal to 2 points after the decimal to keep
        the origin time the same as that saved into the CMTSOLUTION file.
        :type append_focal_mechanism: bool
        :param append_focal_mechanism: will also run append_focal_mechanism()
        :rtype: obspy.core.event.Event
        :return: event retrieved either via internal or external methods
        """
        # If Gatherer has previously found an event, do nothing
        if self.event:
            return self.event
        # If dataset is given, search for event
        elif self.ds:
            try:
                self.event = self.asdf_event_fetch()
            except (AttributeError, IndexError):
                self.event = self.event_get()
                if append_focal_mechanism:
                    self.append_focal_mechanism()
                if self.config.save_to_ds:
                    self.ds.add_quakeml(self.event)
                    logger.debug(
                        "event got from external, added to pyasdf dataset")
                else:
                    logger.debug("event not saved")
        # Else, query FDSN for event information
        else:
            self.event = self.ext_getter.event_get()
            if append_focal_mechanism:
                self.append_focal_mechanism()
            logger.debug("event got from external")

        return self.event

    def append_focal_mechanism(self, overwrite=False):
        """
        Focal mechanisms are unforunately not something that is standard to
        store, as calculations differ between regions and agencies.

        That means this function is not always bulletproof.

        FDSN fetched events are devoid of a few bits of information that are
        useful for our applications, e.g. moment tensor, focal mechanisms.
        This function will perform the conversions and append the necessary
        information to the event located in the dataset.

        :type overwrite: bool
        :param overwrite: If the event already has a focal mechanism, this will
            overwrite that focal mechanism
        """
        if isinstance(self.event, obspy.core.event.Event):
            # If the event already has a focal mechanism attribute, don't gather
            if hasattr(self.event, 'focal_mechanisms') and \
                    self.event.focal_mechanisms and not overwrite:
                return
            if self.config.client == "GEONET":
                # Search GeoNet moment tensor catalog, query GitHub repo
                geonet_mtlist = grab_geonet_moment_tensor(self.config.event_id)
                self.event, _ = generate_focal_mechanism(mtlist=geonet_mtlist,
                                                         event=self.event
                                                         )
                logger.info(
                    "appending GeoNet moment tensor information to event")
            else:
                # Query GCMT .ndk files for matching event
                try:
                    self.event = grab_gcmt_moment_tensor(
                        datetime=self.event.preferred_origin().time,
                        magnitude=self.event.preferred_magnitude().mag
                    )
                except FileNotFoundError:
                    logger.info("No GCMT event found")

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
            return self.int_fetcher.station_fetch(station_code)
        except FileNotFoundError:
            logger.debug(
                "internal station information not found, searching ext.")
            inv = self.ext_getter.station_get(station_code)
            if self.ds is not None and self.config.save_to_ds:
                self.ds.add_stationxml(inv)
            else:
                logger.debug("station information is not being saved")
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
            st_obs = self.int_fetcher.obs_waveform_fetch(station_code)
        except FileNotFoundError:
            logger.debug("internal obs data unavailable, searching external")
            try:
                st_obs = self.ext_getter.waveform_get(station_code)
            # Catch all FDSN Exceptions
            except obspy.clients.fdsn.header.FDSNException:
                logger.warning("external obs data unavailable, no observed "
                               "stream can be returned")
                st_obs = None
            if (self.ds is not None) and self.config.save_to_ds and (
                    st_obs is not None):
                self.ds.add_waveforms(waveform=st_obs,
                                      tag=self.config.observed_tag)
            else:
                logger.debug("observed waveforms are not being saved")

        return st_obs

    def gather_synthetic(self, station_code):
        """
        Gather synthetic data. Can only be provided internally via fetcher

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
            st_syn = self.int_fetcher.syn_waveform_fetch(station_code)
            return st_syn
        except FileNotFoundError:
            return None


