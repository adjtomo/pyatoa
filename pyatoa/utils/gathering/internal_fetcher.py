#!/usr/bin/env python3
"""
Class used to fetch data via internal pathways. Smart enough to know if an
event sits too close to a separation in files (i.e. midnight for waveforms)
and accomodates accordingly. Also will save any fetched data into a Pyasdf
dataset if given.

Hardcoded directory structure and synthetic file name format.
"""
import os
import glob
import copy

from obspy import Stream, Inventory, read, read_inventory

from pyatoa import logger
from pyatoa.utils.operations.calculations import overlapping_days
from pyatoa.utils.operations.conversions import ascii_to_mseed


class Fetcher:
    def __init__(self, config, ds=None, origintime=None):
        """
        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset for internal data searching and saving
        :type origintime: obspy.core.UTCDateTime
        :param origintime: event origin time for event picking
        :type start_pad: int
        :param start_pad: padding in seconds before the origin time of an event
            for waveform fetching, to be fed into lower level functions.
        :type end_pad: int
        :param end_pad: padding in seconds after the origin time of an event
            for wavefomr fetching.
        """
        self.config = config
        self.ds = ds
        self.origintime = origintime

    def asdf_event_fetch(self):
        """
        return event information from pyasdf

        :rtype event: obspy.core.event.Event
        :return event: event object
        """
        event = self.ds.events[0]
        self.origintime = event.origins[0].time
        return event

    def _asdf_station_fetch(self, station_code):
        """
        return station information from pyasdf based on station tag

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.network.Network
        :return: network containing relevant station information
        """
        net, sta, _, _ = station_code.split('.')
        return self.ds.waveforms['{n}_{s}'.format(n=net, s=sta)].StationXML

    def _asdf_waveform_fetch(self, station_code, tag):
        """
        return stream based on tags from pyasdf

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type tag: str
        :param tag: internal asdf tag labelling waveforms
        :rtype: obspy.core.stream.Stream
        :return: waveform contained in a stream
        """
        net, sta, _, _ = station_code.split('.')
        return self.ds.waveforms['{n}_{s}'.format(n=net, s=sta)][tag]

    def _fetch_response(self, station_code, paths_to_responses=None):
        """
        Fetch station xml from given internal pathing. Search through all
        paths given until corresponding inventories are found or until nothing
        is found

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type paths_to_responses: list of str
        :param paths_to_responses: absolute pathways for response file locations
        :rtype inv: obspy.core.inventory.Inventory
        :return inv: inventory containing relevant network and stations
        """
        if not paths_to_responses:
            paths_to_responses = self.config.paths['responses']

        net, sta, loc, cha = station_code.split('.')
        dir_structure = '{sta}.{net}'
        file_template = 'RESP.{net}.{sta}.{loc}.{cha}'

        inv = None
        for path_ in paths_to_responses:
            if not os.path.exists(path_):
                continue
            # Inventory() requires some positional arguements, use None to skip
            inv = None
            fid = os.path.join(path_, dir_structure, file_template).format(
                net=net, sta=sta, cha=cha, loc=loc)
            for filepath in glob.glob(fid):
                if inv is None:
                    inv = read_inventory(filepath)
                    logger.info("response found at {}".format(filepath))
                inv += read_inventory(filepath)

        if inv:
            return inv
        else:
            logger.info("No response found for given paths")
            raise FileNotFoundError()

    def _fetch_by_directory(self, station_code, paths_to_waveforms=None):
        """
        waveform directory structure formatting is hardcoded in here, it is
        assumed that data is saved as miniseeds in the following structure:

        path/to/data/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/{FID}
        e.g. path/to/data/2017/NZ/OPRZ/HHZ.D/NZ.OPRZ.10.HHZ.D

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type paths_to_waveforms: list of str
        :param paths_to_waveforms: absolute pathways for mseed file locations
        :rtype stream: obspy.core.stream.Stream
        :return stream: stream object containing relevant waveforms
        """
        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")
        if not paths_to_waveforms:
            paths_to_waveforms = self.config.paths['waveforms']

        net, sta, loc, cha = station_code.split('.')
        dir_structure = '{year}/{net}/{sta}/{cha}*'
        file_template = '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
        jdays = overlapping_days(origin_time=self.origintime,
                                 start_pad=self.config.start_pad,
                                 end_pad=self.config.end_pad
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
            if len(st) > 0:  # is this necessary?
                st.merge()
                st.trim(starttime=self.origintime-self.config.start_pad,
                        endtime=self.origintime+self.config.end_pad
                        )
                return st
        else:
            logger.info(
                "No waveforms found for {} for given directories".format(
                    station_code)
            )
            raise FileNotFoundError()

    def _fetch_by_event(self, station_code):
        """
        Synthetic data is saved by event id and model number,
        give a hardcoded search path to get to the data.
        If data hasn't been converted to miniseed format,
        call on a conversion function to get from SPECFEM3D-native ascii to
        an obspy stream object in miniseed format

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type paths_to_waveforms: list of str
        :param paths_to_waveforms: absolute pathways for mseed file locations
        :rtype stream: obspy.core.stream.Stream
        :return stream: stream object containing relevant waveforms
        """
        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")

        # Specfem denotes the units of synthetic traces by the naming of the
        # ASCII files, corresponding to 'd' for displacement, 'v' for velocity,
        # and 'a' for acceleration. Take this value from the user config.
        specfem_id = self.config.synthetic_unit[0].lower()

        # Generate information necessary to search for data
        net, sta, _, cha = station_code.split('.')
        specfem_fid_template = '{net}.{sta}.*{cmp}.sem{dva}'

        # Check through paths given in Config
        for path_ in self.config.paths['synthetics']:
            if not os.path.exists(path_):
                continue
            full_path = os.path.join(path_, specfem_fid_template)
            st = Stream()
            for filepath in glob.glob(full_path.format(
                    net=net, sta=sta, cmp=cha[2:], dva=specfem_id)):
                try:
                    # Convert the ASCII file to a miniseed
                    st += ascii_to_mseed(filepath, self.origintime)
                except UnicodeDecodeError:
                    # If the data file is for some reason already in miniseed
                    st += read(filepath)
                logger.info("stream fetched by event {}".format(
                    os.path.basename(filepath))
                )

            if len(st) > 0:
                st.merge()
                st.trim(starttime=self.origintime - self.config.start_pad,
                        endtime=self.origintime + self.config.end_pad
                        )
                return st
        else:
            logger.info(
                "No synthetic waveforms for {} found for given event".format(
                    station_code)
            )
            raise FileNotFoundError()

    def obs_waveform_fetch(self, station_code):
        """
        Main waveform fetching function for observation data.
        Will return a FileNotFoundError if no internal data is found.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type tag: str
        :param tag: The tag to save the waveform by, if an ASDF dataset is given
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        if self.ds:
            try:
                return self._asdf_waveform_fetch(station_code, tag=tag)
            except KeyError:
                st_obs = self._fetch_by_directory(station_code)
                self.ds.add_waveforms(waveform=st_obs,
                                      tag=self.config.observed_tag)
                return st_obs
        else:
            return self._fetch_by_directory(station_code)

    def syn_waveform_fetch(self, station_code):
        """
        If synthetics are freshly generated from specfem, they should
        be placed into a folders separated by event id and model iteration.
        Check if they are already saved into a pyasdf dataset first

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        if self.ds:
            try:
                return self._asdf_waveform_fetch(station_code,
                                                 tag=self.config.synthetic_tag)
            except KeyError:
                st_syn = self._fetch_by_event(station_code)
                self.ds.add_waveforms(waveform=st_syn,
                                      tag=self.config.synthetic_tag)
                return st_syn
        else:
            return self._fetch_by_event(station_code)

    def station_fetch(self, station_code):
        """
        Main station fetching function for observation data.
        Will raise a FileNotFoundError if no internal station data is found

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        if self.ds:
            try:
                return self._asdf_station_fetch(station_code)
            except (KeyError, AttributeError):
                inv = self._fetch_response(station_code)
                self.ds.add_stationxml(inv)
                return inv
        else:
            return self._fetch_response(station_code)










