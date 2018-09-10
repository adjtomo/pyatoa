"""
Class used to fetch data via internal pathways. Smart enough to know if an
event sits too close to a separation in files and can accomodate accordingly.
Hardcoded directory structure and synthetic file name format.
"""
import os
import csv
import glob
import copy

from obspy import Stream, Inventory, read, read_events, read_inventory, \
    UTCDateTime

from pyatoa import logger
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
                logger.info("response found at {}".format(filepath))
                if inv is None:
                    inv = read_inventory(filepath)
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
                    logger.info("stream from directory from {}".format(
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

        st = None
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
                logger.info("stream by event from {}".format(filepath))
            if (st is not None) and (len(st) > 0):
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
            try:
                return self.asdf_waveform_fetch(station_code, tag='observed')
            except KeyError:
                st_obs = self.fetch_by_directory(station_code)
                self.ds.add_waveforms(waveform=st_obs,tag='observed')
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
            try:
                return self.asdf_waveform_fetch(
                    station_code,tag='synthetics_{}'.format(
                    self.config.model_number))
            except KeyError:
                st_syn = self.fetch_by_event(station_code)
                self.ds.add_waveforms(waveform=st_syn,tag='synthetic_{}'.format(
                                                    self.config.model_number))
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
                ds.add_stationxml(self.inv)
                return inv
        else:
            return self.fetch_response(station_code)

    def geonet_moment_tensor_fetch(self):
        """
        fetch moment tensor information from an internal csv file, only relevant
        to my new zealand tomography problem
        :return:
        """
        if self.config.paths["moment_tensors"] is None:
            return
        with open(self.config.paths["moment_tensors"]) as f:
            reader = csv.reader(f)
            for i, row in enumerate(reader):
                if i == 0:
                    tags = row
                if self.config.event_id == row[0]:
                    values = []
                    for t, v in zip(tags, row):
                        if (t == "Date") or (t == "PublicID"):
                            values.append(v)
                        else:
                            values.append(float(v))

                    moment_tensor = dict(zip(tags, values))
                    return moment_tensor

    def gcmt_fetch(self):
        """
        fetch global centroid moment tensor information from internal ndk files,

        :return:
        """
        if self.config.paths["gcmt_ndk"] is None:
            return

        month_dict = {1: "jan", 2: "feb", 3: "mar", 4: "apr", 5: "may",
                      6: "jun", 7: "jul", 8: "aug",9: "sep", 10: "oct",
                      11: "nov", 12: "dec"}
        moment_tensor = self.geonet_moment_tensor_fetch()
        mw = moment_tensor["Mw"]
        date = UTCDateTime(moment_tensor["Date"])
        year = str(date.year)
        month = month_dict[date.month]

        fid = "{m}{y}.ndk".format(m=month, y=year[2:])
        fpath = os.path.join(self.config.paths["gcmt_ndk"], "GCMT", year, fid)
        try:
            cat = read_events(fpath)
        except FileNotFoundError:
            # import HTTPerror to be able to catch it
            from urllib.error import HTTPError
            try:
                cat = read_events(
                    "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/""
                    "catalog/NEW_MONTHLY/{y}/{fid}".format(y=year, fid=fid)
                    )
            except HTTPError:
                cat = read_events(
                    "http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
                    "catalog/NEW_QUICK/qcmt.ndk"
                )

        cat_filt = cat.filter("time > {}".format(str(date - 60)),
                              "time < {}".format(str(date + 60)),
                              "magnitude >= {}".format(mw - .5),
                              "magnitude <= {}".format(mw + .5)
                              )
        if not len(cat_filt):
            raise FileNotFoundError("No events found")
        elif len(cat_filt) > 1:
            print("{} events found, choose from list:".format(len(cat_filt)))
            print("{0}\n{1}".format(moment_tensor, cat_filt))
            choice = int(input("Event number (index from 0): "))
            return cat_filt[choice]
        else:
            return cat_filt[0]








