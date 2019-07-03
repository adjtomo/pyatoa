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

from obspy import Stream, read, read_inventory

from pyatoa import logger
from pyatoa.utils.operations.calculations import overlapping_days
from pyatoa.utils.operations.conversions import ascii_to_mseed
from pyatoa.utils.operations.source_receiver import merge_inventories


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
        """
        self.config = config
        self.ds = ds
        self.origintime = origintime

    def asdf_event_fetch(self):
        """
        Return event information from pyasdf.

        Assumes that the ASDF Dataset will only contain one event, which is
        dictated by the structure of Pyatoa.

        :rtype event: obspy.core.event.Event
        :return event: event object
        """
        event = self.ds.events[0]
        self.origintime = event.preferred_origin().time
        return event

    def asdf_station_fetch(self, station_code):
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
        net, sta, loc, cha = station_code.split('.')
        return self.ds.waveforms[
            '{n}_{s}'.format(n=net, s=sta)].StationXML.select(channel=cha)

    def asdf_waveform_fetch(self, station_code, tag):
        """
        Return stream based on tags from pyasdf. Allows for wildcard selection
        of component. Does not select by channel because synthetic channels
        will differ from observation channels. Select only by component.

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
        net, sta, loc, cha = station_code.split('.')
        return self.ds.waveforms[
            '{n}_{s}'.format(n=net, s=sta)][tag].select(component=cha[2:])

    def fetch_resp_by_dir(self, station_code, paths_to_responses=None,
                          dir_structure='{sta}.{net}',
                          file_template='RESP.{net}.{sta}.{loc}.{cha}'):
        """
        Fetch station xml from given internal pathing. Search through all
        paths given until corresponding inventories are found or until nothing
        is found.

        Note: Obspy does not have any type of inventory merge property, and
            because SEED format saves response by channel, we must read
            individual stationxml files in as individual inventory objects, and
            add them together. This means that the output inventory has multiple
            networks that are essentially the same. It works, but be careful
            when looping through an inventory object.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type paths_to_responses: list of str
        :param paths_to_responses: absolute pathways for response file locations
        :type dir_structure: str
        :param dir_structure: a hardcoded directory structure to search for
            response files. Follows the SEED convention
        :type file_template: str
        :param file_template: a hardcoded file naming template to search for
            response files. Follows the SEED convention
        :rtype inv: obspy.core.inventory.Inventory
        :return inv: inventory containing relevant network and stations
        """
        if not paths_to_responses:
            paths_to_responses = self.config.cfgpaths['responses']

        net, sta, loc, cha = station_code.split('.')

        inv = None
        for path_ in paths_to_responses:
            if not os.path.exists(path_):
                continue
            # Inventory() requires some positional arguements, use None to skip
            fid = os.path.join(path_, dir_structure, file_template).format(
                net=net, sta=sta, cha=cha, loc=loc)
            for filepath in glob.glob(fid):
                # The first inventory becomes the main inv to return
                if inv is None:
                    inv = read_inventory(filepath)
                # All other inventories are appended to the original
                else:
                    inv_append = read_inventory(filepath)
                    inv = merge_inventories(inv, inv_append)

                logger.debug("response found at {}".format(filepath))

        # Merge inventory objects
        if inv is None:
            logger.debug("no response found for given paths")
            raise FileNotFoundError()

        return inv

    def fetch_obs_by_dir(
            self, station_code, paths_to_waveforms=None,
            dir_structure='{year}/{net}/{sta}/{cha}*',
            file_template='{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'):
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
        :type dir_structure: str
        :param dir_structure: a hardcoded directory structure to search for
            observation data. Follows the SEED convention
        :type file_template: str
        :param file_template: a hardcoded file naming template to search for
            observation data. Follows the SEED convention
        :rtype stream: obspy.core.stream.Stream
        :return stream: stream object containing relevant waveforms
        """
        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")
        if paths_to_waveforms is None:
            paths_to_waveforms = self.config.cfgpaths['waveforms']

        net, sta, loc, cha = station_code.split('.')
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
                    logger.debug("stream fetched from directory {}".format(
                        filepath))
            if len(st) > 0:  # is this necessary?
                st.merge()
                st.trim(starttime=self.origintime-self.config.start_pad,
                        endtime=self.origintime+self.config.end_pad
                        )
                return st
        else:
            logger.debug(
                "no waveforms found for {} for given directories".format(
                    station_code)
            )
            raise FileNotFoundError()

    def fetch_syn_by_dir(self, station_code, event_id='',
                         specfem_fid_template='{net}.{sta}.*{cmp}.sem{dva}'):
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
        :type specfem_fid_template: str
        :param specfem_fid_template: The naming template of Specfem ascii files
        :type event_id: str
        :param event_id: The event id, if given, this function will search
            the directory for the event id first, rather than just searching the
            given directory. Useful for when data is stored by event id
            e.g. /path/to/synthetics/{event_id}/*
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
        net, sta, loc, cha = station_code.split('.')

        # Check through paths given in Config
        for path_ in self.config.cfgpaths['synthetics']:
            if not os.path.exists(path_):
                continue

            # Here the path is determined for search. If event_id is given,
            # the function will search for an event_id directory.
            full_path = os.path.join(path_, event_id, specfem_fid_template)

            st = Stream()
            for filepath in glob.glob(full_path.format(
                    net=net, sta=sta, cmp=cha[2:], dva=specfem_id)):
                try:
                    # Convert the ASCII file to a miniseed
                    st += ascii_to_mseed(filepath, self.origintime)
                except UnicodeDecodeError:
                    # If the data file is for some reason already in miniseed
                    st += read(filepath)
                logger.debug("stream fetched by event {}".format(
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
                "no synthetic waveforms for {} found for given event".format(
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
                # Search the given asdf dataset first
                logger.debug("fetching obs internal asdf")
                return self.asdf_waveform_fetch(station_code,
                                                tag=self.config.observed_tag)
            except KeyError:
                logger.debug("obs internal asdf failed, fetching by dir")
                # If asdf dataset does not contain data, search internal dir.
                st_obs = self.fetch_obs_by_dir(station_code)
                self.ds.add_waveforms(waveform=st_obs,
                                      tag=self.config.observed_tag)
                return st_obs
        else:
            return self.fetch_obs_by_dir(station_code)

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
                logger.debug("fetching syn internal asdf")
                return self.asdf_waveform_fetch(station_code,
                                                tag=self.config.synthetic_tag)
            except KeyError:
                logger.debug("syn internal asdf failed, fetching by dir")
                st_syn = self.fetch_syn_by_dir(station_code)
                self.ds.add_waveforms(waveform=st_syn,
                                      tag=self.config.synthetic_tag)
                return st_syn
        else:
            return self.fetch_syn_by_dir(station_code)

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
                logger.debug("searching station internal asdf")
                return self.asdf_station_fetch(station_code)
            except (KeyError, AttributeError):
                logger.debug("internal asdf station ")
                inv = self.fetch_resp_by_dir(station_code)
                self.ds.add_stationxml(inv)
                return inv
        else:
            return self.fetch_resp_by_dir(station_code)










