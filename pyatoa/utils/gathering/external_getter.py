"""
Wrapper class used to quickly interact with obspy fdsn Client, based on configs.
Defaults are set up for local earthquake trace lengths taken from GeoNet fdsn
webservice.
"""
import copy
from obspy.clients.fdsn import Client


class Getter():
    def __init__(self, config=None, client="GEONET", event_id=None,
                 startpad=20, endpad=200):
        self.config = copy.deepcopy(config)
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
        return event information parameters pertaining to a given eventid
        :return:
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
        """
        if self.Client is None:
            self.Client = Client(self.client)
        net,sta,loc,cha = station_code.split('.')
        return self.Client.get_stations(network=net, station=sta, location=loc,
            channel=cha, starttime=self.origintime-self.startpad,
            endtime=self.origintime+self.endpad, level=level
            )

    def waveform_get(self, station_code):
        """
        call for obspy to download data. for some reason obspy can return traces
        with varying sample lengths, so 10 second cushion and then trim after
        retrieval to make sure traces are the same length
        :return:
        """
        if self.Client is None:
            self.Client = Client(self.client)
        net,sta,loc,cha = station_code.split('.')
        st = self.Client.get_waveforms(network=net, station=sta, location=loc,
            channel=cha, starttime=self.origintime-(self.startpad+10),
            endtime=self.origintime+(self.endpad+10)
            )
        st.trim(starttime=self.origintime-self.startpad,
                endtime=self.origintime+self.endpad)
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
        return event,inv,st



