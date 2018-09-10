"""
Class used to fetch data via internal pathways. Mostly a container for functions
small number of attributes to ease in the fetching process
"""
import os
import glob
import copy

from obspy import Stream, read

from pyatoa.utils.operations.calculations import overlapping_days
from pyatoa.utils.operations.conversions import ascii_to_mseed


class Fetcher():
    def __init__(self, config, ds=None, origintime = None, startpad=20,
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
        net,sta,_,_ = station_code.split('.')
        return self.ds.waveforms['{n}_{s}'.format(n=net, s=sta)].StationXML

    def asdf_waveform_fetch(self, station_code, tag):
        """
        return stream based on tags from pyasdf
        :param station_code:
        :param tag:
        :return:
        """
        net,sta,_,_ = station_code.split('.')
        return self.ds.waveforms['{n}_{s}'.format(n=net, s=sta)][tag]


    def fetch_by_directory(self, station_code, paths_to_waveforms=None):
        """
        waveform directory structure is hardcoded in here, we assume that
        data is saved as miniseeds in the following directory structure

        path/to/data/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/*.mseed
        e.g. path/to/data/2017/NZ/OPRZ/HHZ.D/NZ.OPRZ.10.HHZ.D.mseed

        :return:
        """
        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")
        if not paths_to_waveforms:
            paths_to_waveforms = self.config.paths['waveforms']
        net, sta, loc, cha = station_code.split('.')

        # hardcoded directory structure and file naming schema
        dir_structure = '{year}/{net}/{sta}/{cha}*'
        file_template = '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'

        for path_ in paths_to_waveforms:
            if not os.path.exists(path_):
                continue
            st = Stream()
            full_path = os.path.join(path_, dir_structure, file_template)

            # check if the event origin is close to midnight
            jdays = overlapping_days(origin_time=self.origintime,
                                     startpad=self.startpad,
                                     endpad=self.endpad
                                     )
            pathlist = []
            for jday in jdays:
                pathlist.append(full_path.format(net=net, sta=sta, cha=cha,
                                                 loc=loc, jday=jday,
                                                 year=self.origintime.year
                                                 ))
            for fid in pathlist:
                for filepath in glob.glob(fid):
                    st += read(filepath)
            if len(st) > 0:
                st.merge()
                st.trim(starttime=self.origintime-self.startpad,
                        endtime=self.origintime+self.endpad
                    )
                break
        else:
            raise FileNotFoundError("No data found for given paths")

        return st

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
        if not paths_to_waveforms:
            paths_to_waveforms = []
            for path_ in self.config.paths['waveforms']:
                paths_to_waveforms.append(os.path.join(path_, 'SPECFEM3D'))
        net, sta, _, cha = station_code.split('.')
        cmp = cha[-1]

        # hardcoded SPECFEM3D filename output template
        filename_template = '{net}.{sta}.*.sem{cmp}*'
        for path_ in paths_to_waveforms:
            if not os.path.exists(path_):
                continue
            st = Stream()
            import ipdb;ipdb.set_trace()
            full_path = os.path.join(path_, self.config.model_number,
                                     self.config.event_id, filename_template)
            for filepath in glob.glob(
                    full_path.format(net=net, sta=sta, cmp=cmp)):
                try:
                    st += read(filepath)
                except TypeError:
                    st += ascii_to_mseed(filepath,self.origintime)
            if len(st) > 0:
                st.merge()
                st.trim(starttime=self.origintime - self.startpad,
                        endtime=self.origintime + self.endpad
                        )
                break
        else:
            raise FileNotFoundError("No data found for given paths")

        return st

    def obs_waveform_fetch(self, station_code):
        """
        grab observation data, search internally via pyasdf, and then via
        local pathways
        :return:
        """
        net, sta, _, _ = station_code.split('.')
        try:
            return asdf_waveform_fetch(station_code, tag='observed')
        except KeyError:
            return fetch_by_directory(station_code)

    def syn_waveform_fetch(self, station_code):
        """
        if synthetics are freshly generated from specfem, they should
        be placed into a folder separated by event id. check if they are
        already saved into a pyasdf dataset first
        :param station_code:
        :return:
        """
        net, sta, _, _ = station_code.split('.')
        try:
            return asdf_waveform_fetch(station_code,
                                        tag='synthetics_{}'.format(
                                            self.config.model_number))
        except KeyError:
            return fetch_by_event(station_code)
