#!/usr/bin/env python
"""
Data fetching class
"""
import copy
from obspy.clients.fdsn import Client

class ExternalDataFetcher():
    """
    A class used to fetch data via internal pathing for geonet
    """
    def __init__(self,cfg,client="GEONET",event=None,station=None,
                 stream=None)
        self.cfg = copy.deepcopy(cfg)
        self.internal = internal
        self.event = event
        self.station = station
        self.stream = stream
        self.client = Client(client)

    def _get_event(self):
        """
        return event information parameters pertaining to a given eventid
        :return:
        """
        cat = self.client.get_events(eventid=self.config.event_id)
        self.event = cat[0]


    def _get_waveform(self):
        """
        call for obspy to download data
        :return:
        """
        st = self.client.get_waveforms(network=self.cfg.)


    def gather_data(self):
        """
        catch all function to grab all data available
        :return:
        """

