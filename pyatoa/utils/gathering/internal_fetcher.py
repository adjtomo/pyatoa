"""
Class used to fetch data via internal pathways. Smart enough to know if an
event sits too close to a separation in files and can accomodate accordingly.
Hardcoded directory structure and synthetic file name format.
"""
import os
import glob
import copy

from obspy import Stream, read, read_inventory

from pyatoa import logger
from pyatoa.utils.operations.calculations import overlapping_days
from pyatoa.utils.operations.conversions import ascii_to_mseed


class Fetcher:
    def __init__(self, config, ds=None, origintime=None, startpad=20,
                 endpad=200):
        self.config = copy.deepcopy(config)
        self.ds = ds
        self.startpad = startpad
        self.endpad = endpad
        self.origintime = origintime

    def asdf_event_fetch(self):
        """
        return event information from pyasdf
        :return:
        """
        event = self.ds.events[0]
        self.origintime = event.origins[0].time
        return event

    def asdf_station_fetch(self, station_code):
        """
        return station information from pyasdf
        :return:
        """
        net, sta, _, _ = station_code.split('.')
        return self.ds.waveforms['{n}_{s}'.format(n=net, s=sta)].StationXML

    def asdf_waveform_fetch(self, station_code, tag):
        """
        return stream based on tags from pyasdf
        :param station_code:
        :param tag:
        :return:
        """
        net, sta, _, _ = station_code.split('.')
        return self.ds.waveforms['{n}_{s}'.format(n=net, s=sta)][tag]

    def fetch_response(self, station_code, paths_to_responses=None):
        """
        fetch station xml from given internal pathing
        """
        if not paths_to_responses:
            paths_to_responses = self.config.paths['responses']

        net, sta, loc, cha = station_code.split('.')
        dir_structure = '{sta}.{net}'
        file_template = 'RESP.{net}.{sta}.{loc}.{cha}'

        for path_ in paths_to_responses:
            if not os.path.exists(path_):
                continue
            fid = os.path.join(path_, dir_structure, file_template)
            inv = None
            for filepath in glob.glob(fid.format(
                                            net=net, sta=sta, cha=cha,loc=loc)):
                if inv is None:
                    inv = read_inventory(filepath)
                    logger.info("response found at {}".format(filepath))
                inv += read_inventory(filepath)
            if inv is not None:
                return inv
        else:
            raise FileNotFoundError("No response found for given paths")

    def fetch_by_directory(self, station_code, paths_to_waveforms=None):
        """
        waveform directory structure is hardcoded in here, we assume that
        data is saved as miniseeds in the following directory structure

        path/to/data/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/{FID}
        e.g. path/to/data/2017/NZ/OPRZ/HHZ.D/NZ.OPRZ.10.HHZ.D

        :return:
        """
        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")
        if not paths_to_waveforms:
            paths_to_waveforms = self.config.paths['waveforms']

        net, sta, loc, cha = station_code.split('.')
        dir_structure = '{year}/{net}/{sta}/{cha}*'
        file_template = '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
        jdays = overlapping_days(origin_time=self.origintime,
                                 startpad=self.startpad,
                                 endpad=self.endpad
                                 )
        for path_ in paths_to_waveforms:
            if not os.path.exists(path_):
                continue
            full_path = os.path.join(path_, dir_structure, file_template)
            pathlist = []
            for jday in jdays:
                pathlist.append(full_path.format(net=net, sta=sta, cha=cha,
                                                 loc=loc, jday=jday,
                                                 year=self.origintime.year)
                                                 )
            st = Stream()
            for fid in pathlist:
                for filepath in glob.glob(fid):
                    st += read(filepath)
                    logger.info("stream fetched from directory {}".format(
                        filepath))
            if len(st) > 0: # is this necessary?
                st.merge()
                st.trim(starttime=self.origintime-self.startpad,
                        endtime=self.origintime+self.endpad
                        )
                return st
        else:
            raise FileNotFoundError("No waveforms found for given directories")

    def fetch_by_event(self, station_code, paths_to_waveforms=None):
        """
        Synthetic data is saved by event code, give a hardcoded search path
        to get to the data. If data hasn't been converted to miniseed format,
        call on a conversion to convert from SPECFEM3D native ascii format
        :param station_code:
        :param paths_to_waveforms:
        :return:
        """
        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")
        if paths_to_waveforms is None:
            paths_to_waveforms = []
            for path_ in self.config.paths['waveforms']:
                paths_to_waveforms.append(os.path.join(path_, 'SPECFEM3D'))
        net, sta, _, cha = station_code.split('.')
        comp = cha[-1]

        specfem_fid_template = '{net}.{sta}.*{cmp}.semv*'
        for path_ in paths_to_waveforms:
            if not os.path.exists(path_):
                continue
            full_path = os.path.join(path_, self.config.model_number,
                                     self.config.event_id, specfem_fid_template)
            st = Stream()
            for filepath in glob.glob(
                                full_path.format(net=net, sta=sta, cmp=comp)):
                try:
                    st += ascii_to_mseed(filepath,self.origintime)
                except UnicodeDecodeError:
                    st += read(filepath)
                logger.info("stream fetched by event {}".format(filepath))
            if len(st) > 0:
                st.merge()
                st.trim(starttime=self.origintime - self.startpad,
                        endtime=self.origintime + self.endpad
                        )
                return st
        else:
            raise FileNotFoundError("No event found for given paths")

    def obs_waveform_fetch(self, station_code):
        """
        grab observation data, search internally via pyasdf, and then via
        local pathways
        :return:
        """
        if self.ds is not None:
            tag = "observed"
            try:
                return self.asdf_waveform_fetch(station_code, tag=tag)
            except KeyError:
                st_obs = self.fetch_by_directory(station_code)
                self.ds.add_waveforms(waveform=st_obs, tag=tag)
                return st_obs
        else:
            return self.fetch_by_directory(station_code)

    def syn_waveform_fetch(self, station_code):
        """
        if synthetics are freshly generated from specfem, they should
        be placed into a folder separated by event id. check if they are
        already saved into a pyasdf dataset first
        :param station_code:
        :return:
        """
        if self.ds is not None:
            tag = "synthetic_{}".format(self.config.model_number)
            try:
                return self.asdf_waveform_fetch(station_code, tag=tag)
            except KeyError:
                st_syn = self.fetch_by_event(station_code)
                self.ds.add_waveforms(waveform=st_syn, tag=tag)
                return st_syn
        else:
            return self.fetch_by_event(station_code)

    def station_fetch(self, station_code):
        """
        grab station response information, search internally first
        """
        if self.ds is not None:
            try:
                return self.asdf_station_fetch(station_code)
            except KeyError:
                inv = self.fetch_response(station_code)
                self.ds.add_stationxml(inv)
                return inv
        else:
            return self.fetch_response(station_code)










