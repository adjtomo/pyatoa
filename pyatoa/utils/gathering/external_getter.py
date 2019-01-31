#!/usr/bin/env python3

"""
Wrapper class used to quickly interact with obspy fdsn Client, based on configs.
Defaults are set up for local earthquake trace lengths taken from GeoNet fdsn
webservice.
"""
import copy
from obspy.clients.fdsn import Client

from pyatoa import logger


class Getter:
    def __init__(self, config=None, client="GEONET", event_id=None,
                 startpad=20, endpad=200):
        """
        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type client: str
        :param client: FDSN client for data getting
        :type event_id: str
        :param event_id: unique event identifier
        :type startpad: int
        :param startpad: padding in seconds before the origin time of an event
            for waveform fetching, to be fed into lower level functions.
        :type endpad: int
        :param endpad: padding in seconds after the origin time of an event
            for wavefomr fetching.
        """

        self.config = config
        if self.config is not None:
            self.event_id = self.config.event_id
        else:
            self.event_id = event_id
        self.client = client
        self.startpad = startpad
        self.endpad = endpad
        self.Client = None
        self.event = None
        self.origintime = None

    def event_get(self):
        """
        return event information parameters pertaining to a given event id

        :rtype event: obspy.core.event.Event
        :return event: event object
        """
        if self.Client is None:
            self.Client = Client(self.client)
        event = self.Client.get_events(eventid=self.event_id)[0]
        self.origintime = event.origins[0].time
        return event

    def station_get(self, station_code, level='response'):
        """
        return station information with level dependent on call, defaults to
        retrieving response information

        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        if self.Client is None:
            self.Client = Client(self.client)
        net, sta, loc, cha = station_code.split('.')
        return self.Client.get_stations(network=net, station=sta, location=loc,
                                        channel=cha,
                                        starttime=self.origintime-self.startpad,
                                        endtime=self.origintime+self.endpad,
                                        level=level
                                        )

    def waveform_get(self, station_code):
        """
        Call for obspy to download data. For some reason obspy can return traces
        with varying sample lengths, so 10 second cushion and then trim after
        retrieval to make sure traces are the same length

        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stre
        """
        if self.Client is None:
            self.Client = Client(self.client)
        net, sta, loc, cha = station_code.split('.')
        st = self.Client.get_waveforms(network=net, station=sta, location=loc,
                                       channel=cha, starttime=self.origintime-(
                                                            self.startpad+10),
                                       endtime=self.origintime+(self.endpad+10)
                                       )
        st.trim(starttime=self.origintime-self.startpad,
                endtime=self.origintime+self.endpad)
        logger.info("stream got external {}".format(station_code))
        return st

    def get_all(self, station_code):
        """
        convenience function for retrieving station, inventory and receiver
        :param station_code:
        :return:
        """
        event = self.event_get()
        inv = self.station_get(station_code)
        st = self.waveform_get(station_code)
        return event, inv, st



