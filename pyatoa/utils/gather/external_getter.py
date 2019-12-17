#!/usr/bin/env python3

"""
Wrapper class used to quickly interact with Obspy FDSN Client, based on Config
"""
from obspy.clients.fdsn import Client

from pyatoa import logger


class Getter:
    def __init__(self, config=None, client="GEONET", event_id=None):
        """
        Initiate the external getter class to retrieve objects from FDSN

        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type client: str
        :param client: FDSN client for data getting
        :type event_id: str
        :param event_id: unique event identifier
        """
        # Priroritize the Config.event_id over a given id, but allow a User
        # to set their own event id to avoid reliance on a config object
        self.config = config
        if self.config is not None:
            self.event_id = self.config.event_id
        else:
            self.event_id = event_id
        self.client = client
        self.event = None
        self.origintime = None

        # Initiate Obspy FDSN client
        if self.client:
            self.Client = Client(self.client)
    
    def event_get(self):
        """
        return event information parameters pertaining to a given event id

        :rtype event: obspy.core.event.Event
        :return event: event object
        """
        logger.debug(f"fetching event from {self.client}")
        event = self.Client.get_events(eventid=self.event_id)[0]
        self.origintime = event.origins[0].time
        return event

    def station_get(self, station_code, level='response'):
        """
        return station information with level dependent on call, defaults to
        retrieving response information

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type level: str
        :param level: level argument to be passed to obspy
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        logger.debug(f"fetching station from {self.client}")
        net, sta, loc, cha = station_code.split('.')
        return self.Client.get_stations(
            network=net, station=sta, location=loc, channel=cha,
            starttime=self.origintime-self.config.start_pad,
            endtime=self.origintime+self.config.end_pad, level=level
        )
 
    def waveform_get(self, station_code):
        """
        Call for obspy to download data. For some reason obspy can return traces
        with varying sample lengths, so 10 second cushion and then trim after
        retrieval to make sure traces are the same length

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stream
        """
        logger.debug(f"fetching observations from {self.client}")
        net, sta, loc, cha = station_code.split('.')
        st = self.Client.get_waveforms(
            network=net, station=sta, location=loc, channel=cha,
            starttime=self.origintime-(self.config.start_pad+5),
            endtime=self.origintime+(self.config.end_pad+5)
        )
        # Sometimes FDSN queries return improperly cut start and end times, so
        # we retrieve +/-10 seconds and then cut down
        st.trim(starttime=self.origintime-self.config.start_pad,
                endtime=self.origintime+self.config.end_pad)

        logger.debug(f"stream got external {station_code}")
        return st

    def get_all(self, station_code):
        """
        Convenience function for retrieving station, inventory and receiver

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype event: obspy.core.event.Event
        :return event: event object
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stream
        """
        event = self.event_get()
        inv = self.station_get(station_code)
        st = self.waveform_get(station_code)
        return event, inv, st



