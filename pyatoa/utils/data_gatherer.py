#!/usr/bin/env python
"""
Data fetching class
"""
import copy
from obspy import Stream, read, read_events, read_inventory
from obspy.clients.fdsn import Client

from .source_receiver import gcd_and_baz, theoretical_p_arrival
from .conv_and_calc import calculate_julian_day

class Gatherer():
    """
    A class used to fetch data via multiple pathways
    """
    def __init__(self,cfg,ds,client="GEONET",startpad=20,endpad=200,event=None,
                 station=None,stream=None):
        self.cfg = copy.deepcopy(cfg)
        self.event = event
        self.station = station
        self.stream = stream
        self.client = Client(client)
        self.startpad = startpad
        self.endpad = endpad
        self.gcd = None
        self.baz = None
        self.p_arrival = None
        self.origintime = None

    def _srcrcv_theoretical_info(self):
        """
        return theoretical information pertaining to source receiver information
        """
        self.gcd, self.baz = gcd_and_baz(self.inventory,self.event)
        self.p_arrival = theoretical_p_arrival(self.event.origins[0].depth*1E-3)

    def _gather_event(self):
        """
        get event information, check internally first
        """
        quakeml_path = os.path.join(self.cfg.internal_paths["quakeml"],
                                    self.cfg.event_id+'.xml')
        if os.path.exists(quakeml_path):
            self.event = read_events(quakeml_path)
        else:
            getter = Getter(gatherer=self)
            self.event = getter._get_event()
        self.origintime = self.event.origins[0].time

    def gather_all(self):
        """
        get waveform information
        """
        try:
            fetcher = Fetcher(gatherer=gatherer)
            self.st_obs, self.st_syn, self_inv = fetcher.fetch_all()
        except Exception as e:
            pass

        try:
            getter = Getter(gatherer=gatherer)
            self.st_obs, self.inv = getter.get_all()

    def _write_to_pyasdf(self,ds):
        """
        save station, event information and raw waveforms into pyasdf
        """
        ds.add_quakeml(self.event)
        ds.add_stationxml(inv)
        ds.add_waveforms(waveform)


class Getter():
    """data getter through fdsn webservies
    """
    def __init__(self,gatherer):
        self.G = gatherer

    def _get_observation_waveform(self):
        """
        call for obspy to download data
        :return:
        """
        return self.G.client.get_waveforms(network=self.G.cfg.network,
            station=self.G.cfg.station,location=self.G.cfg.location,
            channel=self.G.cfg.channel,starttime=self.G.event.-self.G.startpad,
            endtime=self.G.event+self.G.endpad
            )

    def _get_event(self):
        """
        return event information parameters pertaining to a given eventid
        :return:
        """
        return self.G.client.get_events(eventid=self.G.config.event_id)[0]

    def _get_station(self):
        """
        return station information with response
        """
        return self.G.client.get_stations(network=self.G.cfg.network,
            station=self.G.cfg.station,location=self.G.cfg.location,
            component=self.G.cfg.component,level="response"
            )

    def get_all(self):
        """
        call internal functions to get requisite data
        """
        st_obs = self._get_observation_waveform()
        inv = self._get_station()
        return st_obs, inv


class Fetcher():
    """data fetcher for internal pathing
    """
    def __init__(self,gatherer):
        self.G = gatherer

    def _fetch_synthetic_waveforms(self):
        """
        retrieve synthetic data from internal pathing

        NOTE: assumes that files have been converted from native ascii
        into mseed. In the future this may be a step that we want to incorporate
        directly into the workflow, and therefore this function will need to
        change
        """
        syn = Stream()
        synfiles = glob.glob(os.path.join(
            self.G.cfg.internal_paths["synthetic"],
            "{n}.{s}.BX?.semv.mseed".format(n=self.G.cfg.network,
                                            s=self.G.cfg.station)
            ))
        for synfile in synfiles:
            syn += read(synfile)
        return syn

    def _fetch_wav_from_geonet(self):
        """
        deprecated
        for internal use: fetch waveform data directly from GeoNet directories
        """
        mseed_template = '{net}.{sta}.{loc}.{cha}.D.{yea}.{jdy}'
        for comp in self.G.cfg.components:
            mseed_path = '/geonet/seismic/{yr}/NZ/{sta}/{ch}.D/'.format(
                year=self.G.origintime.year,sta=self.G.station,
                ch=self.G.channel)
            juldays = calculate_julian_day(origin_time=self.G.origintime,
                startpad=self.G.startpad,endpad=self.G.endpad
                )
            st = Stream()
            for j in juldays:
                mseed_file = glob.glob(os.path.join(mseed_path,'*{}'.format(j)))
                st += read(mseed_file[0])
        st.trim(starttime=self.G.origintime-self.G.startpad,
                endtime=self.G.origintime+self.G.endpad)
        return st

    def _fetch_response_from_geonet(self):
        """
        for internal use: fetch response data directly from GeoNet directories
        """
        inv_paths = glob.glob('/geonet/seed/RESPONSE/{sta}.NZ/*{cha}'.format(
            sta=self.G.station,cha=self.G.channel)
            )
        inv = read_inventory(inv_paths[0])
        for i in inv_paths[1:]:
            inv += read_inventory(i)
        return inv

    def _fetch_fathom_response(self):
        """
        for grabbing internal fathom response information
        """
        inv = read_inventory('{base}/FATHOM/DATALESS/XX.RDF.DATALESS'.format(
            base=self.G.basepath))
        return inv.select(station=self.G.cfg.station)

    def _fetch_internal_waveforms(self):
        """
        for internal use: fetch waveform data directly from GeoNet directories
        """
        path_template = '{base}/{year}/{net}/{sta}/{cha}.D/'
        mseed_template = '{net}.{sta}.{loc}.{cha}.D.{yea}.{jdy}'

        for comp in self.G.cfg.components:
            mseed_path = path_template.format(base=self.G.basepath,
                year=self.G.origintime.year,sta=self.G.station,
                ch=self.G.channel)
            juldays = calculate_julian_day(origin_time=self.G.origintime,
                startpad=self.G.startpad,endpad=self.G.endpad
                )
            st = Stream()
            for j in juldays:
                mseed_file = glob.glob(os.path.join(mseed_path,'*{}'.format(j)))
                st += read(mseed_file[0])
        st.trim(starttime=self.G.origintime-self.G.startpad,
                endtime=self.G.origintime+self.G.endpad)
        return st

    def fetch_all(self):
        """
        attempt to internally grab all data through internal pathing
        """
        st_obs = _fetch_wav_internal(self)
        st_syn = _fetch_synthetic_waveforms(self)
        if self.G.network == "NZ":
            inv = _fetch_response_from_geonet(self)
        elif self.G.network == "XX":
            inv = _fetch_fathom_response(self)
        else:
            print("Inventory fetching for networks outside"
                                                    "NZ and XX not implemented")
        return st_obs, st_syn, inv

















