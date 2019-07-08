#!/usr/bin/env python3

"""
Wrapper class used to quickly interact with obspy fdsn Client, based on configs.
Defaults are set up for local earthquake trace lengths taken from GeoNet fdsn
webservice.
"""
from obspy.clients.fdsn import Client

from pyatoa import logger


class Getter:
    def __init__(self, config=None, client="GEONET", event_id=None):
        """
        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type client: str
        :param client: FDSN client for data getting
        :type event_id: str
        :param event_id: unique event identifier
        :type start_pad: int
        :param start_pad: padding in seconds before the origin time of an event
            for waveform fetching, to be fed into lower level functions.
        :type end_pad: int
        :param end_pad: padding in seconds after the origin time of an event
            for wavefomr fetching.
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

    def _parse_channel(value):
        """
        fdsn clients cannot handle [xyz] unix style wildcards, although it can
        handle * and ?, small function to address this so that these wildcard
        styles can still be used in pyatoa
        """
        raise NotImplementedError

        if "[" not in value:
             return value
        else:
            code, components = value.split('[')  
            components = components.split(']')[0]
            # If only one component given in wildcard, return 
            if len(components) == 1:
                return code + components
            # If three components given in wildcard, return unix '?' wildcard
            elif len(components) == 3:
                return code + "?"
            # Else, return a list of components to add together
            else:
                return_channels = []
                for component in components:
                    return_channels.append(code + component)
                return return_channels 
        
    
    def event_get(self):
        """
        return event information parameters pertaining to a given event id

        :rtype event: obspy.core.event.Event
        :return event: event object
        """
        logger.debug("fetching event from {}".format(self.client))
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
        logger.debug("fetching station from {}".format(self.client))
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

        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stre
        """
        logger.debug("fetching observations from {}".format(self.client))
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

        logger.debug("stream got external {}".format(station_code))
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



