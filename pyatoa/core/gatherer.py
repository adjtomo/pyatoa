#!/usr/bin/env python
"""
Mid and Low level data gathering classes to retrieve data from local filesystems
or to query data from FDSN webservices via ObsPy.

.. note::
    * Fetch: used to denote searching local filesystems for data (internal)
    * Get: used to denote querying FDSN webservices via ObsPy (external), naming
      convention from the obspy.fdsn.client.Client function names


Gatherer directly called by the Manager class and shouldn't need to be called
by the User unless for bespoke data gathering functionality.
"""
import os
import glob
import warnings
import traceback

from pyasdf import ASDFWarning
from obspy.core.event import Event
from obspy.clients.fdsn import Client
from obspy import Stream, read, read_inventory, read_events
from obspy.clients.fdsn.header import FDSNException

from pyatoa import logger
from pyatoa.utils.read import (read_sem, read_specfem2d_source, 
                               read_forcesolution)
from pyatoa.utils.form import format_event_name
from pyatoa.utils.calculate import overlapping_days
from pyatoa.utils.srcrcv import merge_inventories


class GathererNoDataException(Exception):
    """
    Custom exception to be thrown generally to Manager class in the case that
    gathering of any data fails.
    """
    pass


class ExternalGetter:
    """
    Low-level gathering classs to retrieve data via FDSN webservices.
    Calls made through ObsPy. Functionality is inhereted by the Gatherer class.

    .. warning::
        Not to be used standalone, class attributes are inherted by the
        mid-level Gatherer class.
    """
    def event_get(self, event_id=None):
        """
        Return event information parameters pertaining to a given event id
        if an event id is given, else by origin time. Catches FDSN exceptions.

        :rtype event: obspy.core.event.Event or None
        :return event: event object if found, else None.
        """
        if not self.Client:
            return None
        if event_id is None:
            event_id = self.config.event_id

        event, origintime = None, None
        if event_id is not None:
            try:
                # Get events via event id, only available from certain clients
                logger.debug(f"event ID: {event_id}, querying "
                             f"client {self.config.client}")
                event = self.Client.get_events(eventid=event_id)[0]
            except FDSNException:
                pass
        if self.origintime and event is None:
            try:
                # If getting by event id doesn't work, try based on origintime
                logger.debug(f"origintime: {self.origintime}, querying"
                             f"client {self.config.client}")
                event = self.Client.get_events(starttime=self.origintime,
                                               endtime=self.origintime
                                               )
                if len(event) > 1:
                    # Getting by origin time may result in multiple events 
                    # found in the catalog, this is hard to control and will
                    # probably need to be addressed manually.
                    logger.warning(f"{len(event)} events found, expected 1."
                                   f"Returning first entry, manual revision "
                                   f"may be required."
                                   )
                event = event[0]
            except FDSNException:
                pass
        return event

    def station_get(self, code, **kwargs):
        """
        Call for ObsPy FDSN client to download station dataless information.
        Defaults to retrieving response information.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations

        Keyword Arguments
        ::
            str station_level:
                The level of the station metadata if retrieved using the ObsPy
                Client. Defaults to 'response'
        """
        level = kwargs.get("station_level", "response")

        if not self.Client:
            return None

        logger.debug(f"querying client {self.config.client}")
        net, sta, loc, cha = code.split('.')
        try:
            inv = self.Client.get_stations(
                network=net, station=sta, location=loc, channel=cha,
                starttime=self.origintime - self.config.start_pad,
                endtime=self.origintime + self.config.end_pad, level=level
            )
            return inv
        except FDSNException:
            return None

    def obs_waveform_get(self, code):
        """
        Call for ObsPy FDSN webservice client to download waveform data.

        .. Note:
            ObsPy sometimes returns traces with varying sample lengths,
            so we use a 10 second cushion on start and end time and trim after
            retrieval to make sure traces are the same length.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stream
        :raises FDSNException: if no data found using Client
        """
        if not self.Client or self.config.synthetics_only:
            return None

        logger.debug(f"querying client {self.config.client}")
        net, sta, loc, cha = code.split('.')
        try:
            st = self.Client.get_waveforms(
                network=net, station=sta, location=loc, channel=cha,
                starttime=self.origintime - (self.config.start_pad + 10),
                endtime=self.origintime + (self.config.end_pad + 10)
            )
            # Sometimes FDSN queries return improperly cut start and end times,
            # so we retrieve +/-10 seconds and then cut down
            st.trim(starttime=self.origintime - self.config.start_pad,
                    endtime=self.origintime + self.config.end_pad)
            return st
        except FDSNException:
            return None

    def _obs_get_multithread(self, code, **kwargs):
        """
        A small function to gather StationXMLs and observed waveforms together.
        Used for multithreading where IO queries are sped up.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type return_count: int
        :param return_count: if not None, determines how many data items must be
            collected for the station to be saved into the ASDFDataSet. e.g.
            StationXML and 3 component waveforms would equal 4 pieces of data,
            so a return_count == 4 means stations that do not return all
            components and metadata will not be saved to the dataset.
        :rtype status: int
        :return status: a simple status check that lets the user know how many
            items were collected.
        """
        level = kwargs.get("station_level", "response")
        return_count = kwargs.get("return_count", None)

        data_count, inv, st = 0, None, None
        net, sta, loc, cha = code.split(".")
        try:
            inv = self.Client.get_stations(
                network=net, station=sta, location=loc, channel=cha,
                starttime=self.origintime - self.config.start_pad,
                endtime=self.origintime + self.config.end_pad, level=level
            )
            data_count += 1
        except FDSNException:
            pass

        try:
            st = self.Client.get_waveforms(
                network=net, station=sta, location=loc, channel=cha,
                starttime=self.origintime - (self.config.start_pad + 10),
                endtime=self.origintime + (self.config.end_pad + 10),
                attach_response=True
            )
            # Sometimes FDSN queries return improperly cut start and end times,
            # so we retrieve +/-10 seconds and then cut down
            st.trim(starttime=self.origintime - self.config.start_pad,
                    endtime=self.origintime + self.config.end_pad
                    )
            data_count += len(st)
        except FDSNException:
            pass

        # Additional check for saving data if not all requested data found
        if (return_count is not None) and (data_count < return_count):
            _save = False
        else:
            _save = True

        # Save data to ASDFDataSet if save criteria are met and dats available
        if _save and inv is not None:
            self.ds.add_stationxml(inv)
        if _save and st is not None:
            self.ds.add_waveforms(waveform=st, tag=self.config.observed_tag)

        return data_count


class InternalFetcher:
    """
    Low-level data gatherer to search for data on disk.
    
    * Initially looks through User provided ASDFDataSets, if not available or no
      data locatable, search through the local directories.
    * Defaults directory structure follows SEED convention but can be
      overwritten.

    .. warning::
        Not to be used standalone, class attributes are inherted by the
        mid-level Gatherer class.
    """
    def asdf_event_fetch(self):
        """
        Return Event information from ASDFDataSet.

        .. note::
            Assumes that the ASDF Dataset will only contain one event, which is
            dictated by the structure of Pyatoa.
        .. note::
            TO DO:
            * Remove the logger statement, and write to corresponding functions
            * event_fetch: calls this function and (1) in succession

        :rtype event: obspy.core.event.Event
        :return event: event object
        :raises AttributeError: if no event attribute found in ASDFDataSet
        :raises IndexError: if event attribute found but no events
        """
        event = self.ds.events[0]
        logger.debug(f"matching event found: {format_event_name(event)}")
        return event

    def asdf_station_fetch(self, code):
        """
        Return StationXML from ASDFDataSet based on station code.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.network.Network
        :return: network containing relevant station information
        :raises KeyError: if no matching StationXML found
        """
        net, sta, loc, cha = code.split(".")
        return self.ds.waveforms[f"{net}_{sta}"].StationXML.select(channel=cha)

    def asdf_waveform_fetch(self, code, tag):
        """
        Return waveforms as Stream objects based from ASDFDataSet.

        .. note:
            * Allows for wildcard selection of component (? or *)
            * Selects by component because synthetic channel naming may differ
            from observation channels.
            * Component is assumed to be the last index in the channel,
            following SEED convention.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type tag: str
        :param tag: internal asdf tag labelling waveforms
        :rtype: obspy.core.stream.Stream
        :return: waveform contained in a stream
        :raises KeyError: if no matching waveforms found.
        """
        net, sta, loc, cha = code.split(".")
        return self.ds.waveforms[f"{net}_{sta}"][tag].select(component=cha[-1])

    def fetch_event_by_dir(self, event_id, prefix="", suffix="", format_=None, 
                           **kwargs):
        """
        Fetch event information via directory structure on disk. Developed to
        parse CMTSOLUTION and QUAKEML files, but theoretically accepts any 
        format that the ObsPy read_events() function will accept.

        Will search through all paths given until a matching source file found.

        .. note::
            This function will search for the following path
            /path/to/event_dir/{prefix}{event_id}{suffix}
            
            so, if e.g., searching for a CMTSOLUTION file in the current dir:
            ./CMTSOLUTION_{event_id}

            Wildcards are okay but the function will return the first match

        :type event_id: str
        :param event_id: Unique event identifier to search source file by.
            e.g., a New Zealand earthquake ID '2018p130600'. A prefix or suffix
            will be tacked onto this 
        :rtype event: obspy.core.event.Event or None
        :return event: event object if found, else None.
        :type prefix: str
        :param prefix Prefix to prepend to event id for file name searching.
            Wildcards are okay.
        :type suffix: str
        :param suffix: Suffix to append to event id for file name searching.
            Wildcards are okay.
        :type format_: str or NoneType
        :param format_: Expected format of the file to read, e.g., 'QUAKEML', 
            passed to ObsPy read_events. NoneType means read_events() will guess
        """
        # Ensure that the paths are a list so that iterating doesnt accidentally
        # try to iterate through a string.
        paths = self.config.paths["events"]
        if not isinstance(paths, list):
            paths = [paths]

        event = None
        for path_ in paths:
            if not os.path.exists(path_):
                continue
            # Search for available event files
            fid = os.path.join(path_, f"{prefix}{event_id}{suffix}")
            for filepath in glob.glob(fid):
                logger.debug(f"searching for event data: {filepath}")
                if os.path.exists(filepath):
                    try:
                        # Allow input of various types of source files
                        if "SOURCE" in prefix:
                            logger.info(f"reading SPECFEM2D SOURCE: {filepath}")
                            cat = [read_specfem2d_source(filepath)]
                        elif "FORCESOLUTION" in prefix:
                            logger.info(f"reading FORCESOLUTION: {filepath}")
                            cat = [read_forcesolution(filepath)]
                        else:
                            logger.info(
                                    f"reading source using ObsPy: {filepath}")
                            cat = read_events(filepath, format=format_)

                        if len(cat) != 1:
                            logger.warning(
                                f"{filepath} event file contains more than one "
                                 "event, returning 1st entry")
                        event = cat[0]
                        break
                    except Exception as e:
                        logger.warning(f"{filepath} event file read error {e}")

        if event is not None:
            logger.info(f"retrieved local file:\n{filepath}")
        else:
            logger.info(f"no local event file found")

        return event

    def fetch_resp_by_dir(self, code, **kwargs):
        """
        Fetch station dataless via directory structure on disk.
        Will search through all paths given until StationXML found.

        .. note::
            Default path naming follows SEED convention, that is:
            path/to/dataless/{NET}.{STA}/RESP.{NET}.{STA}.{LOC}.{CHA}
            e.g. path/to/dataless/NZ.BFZ/RESP.NZ.BFZ.10.HHZ

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype inv: obspy.core.inventory.Inventory or None
        :return inv: inventory containing relevant network and stations

        Keyword Arguments
        ::
            str resp_dir_template:
                Directory structure template to search for response files.
                By default follows the SEED convention,
                'path/to/RESPONSE/{sta}.{net}/'
            str resp_fid_template:
                Response file naming template to search for station dataless.
                By default, follows the SEED convention
                'RESP.{net}.{sta}.{loc}.{cha}'
        """
        resp_dir_template = kwargs.get("resp_dir_template", "{sta}.{net}")
        resp_fid_template = kwargs.get("resp_fid_template",
                                       "RESP.{net}.{sta}.{loc}.{cha}")

        inv = None
        net, sta, loc, cha = code.split('.')

        # Ensure that the paths are a list so that iterating doesnt accidentally
        # try to iterate through a string.
        paths = self.config.paths["responses"]
        if not isinstance(paths, list):
            paths = [paths]

        for path_ in paths:
            if not os.path.exists(path_):
                continue
            # Attempting to instantiate an empty Inventory requires some 
            # positional arguements we dont have, so don't do that
            fid = os.path.join(path_, resp_dir_template, resp_fid_template)
            fid = fid.format(net=net, sta=sta, cha=cha, loc=loc)
            logger.debug(f"searching for responses: {fid}")
            for filepath in glob.glob(fid):
                if inv is None:
                    # The first inventory becomes the main inv to return
                    inv = read_inventory(filepath)
                else:
                    # All other inventories are appended to the original
                    inv_append = read_inventory(filepath)
                    # Merge inventories to remove repeated networks
                    inv = merge_inventories(inv, inv_append)
                logger.info(f"retrieved response locally:\n{filepath}")

        return inv

    def fetch_obs_by_dir(self, code, **kwargs):
        """
        Fetch observation waveforms via directory structure on disk.

        .. note::
            Default waveform directory structure assumed to follow SEED
            convention. That is:
            path/to/data/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/{FID}
            e.g. path/to/data/2017/NZ/OPRZ/HHZ.D/NZ.OPRZ.10.HHZ.D

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype stream: obspy.core.stream.Stream or None
        :return stream: stream object containing relevant waveforms, else None

        Keyword Arguments
        ::
            str obs_dir_template:
                directory structure to search for observation data.
                Follows the SEED convention:
                'path/to/obs_data/{year}/{net}/{sta}/{cha}'
            str obs_fid_template:
                File naming template to search for observation data.
                Follows the SEED convention:
                '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
        """
        obs_dir_template = kwargs.get("obs_dir_template",
                                   "{year}/{net}/{sta}/{cha}*")
        obs_fid_template = kwargs.get("obs_fid_template",
                                   "{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}")

        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")

        net, sta, loc, cha = code.split('.')
        # If waveforms contain midnight, multiple files need to be read
        jdays = overlapping_days(origin_time=self.origintime,
                                 start_pad=self.config.start_pad,
                                 end_pad=self.config.end_pad
                                 )

        # Ensure that the paths are a list so that iterating doesnt accidentally
        # try to iterate through a string.
        paths = self.config.paths["waveforms"]
        if not isinstance(paths, list):
            paths = [paths]

        for path_ in paths:
            if not os.path.exists(path_):
                continue
            full_path = os.path.join(path_, obs_dir_template, obs_fid_template)
            pathlist = []
            for jday in jdays:
                pathlist.append(full_path.format(net=net, sta=sta, cha=cha,
                                                 loc=loc, jday=jday,
                                                 year=self.origintime.year)
                                )
            st = Stream()
            for fid in pathlist:
                logger.debug(f"searching for observations: {fid}")
                for filepath in glob.glob(fid):
                    st += read(filepath)
                    logger.info(f"retrieved observations locally:\n{filepath}")
            if len(st) > 0:
                # Take care of gaps in data by converting to masked data
                st.merge()
                st.trim(starttime=self.origintime-self.config.start_pad,
                        endtime=self.origintime+self.config.end_pad
                        )
                # Check if trimming retains data
                if len(st) > 0:
                    return st   
                else:
                    logger.warning("data does not fit origin time +/- pad time")
                    return None
        else:
            return None

    def fetch_syn_by_dir(self, code, **kwargs):
        """
        Fetch synthetic waveforms from Specfem3D via directory structure on
        disk, if necessary convert native ASCII format to Stream object.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype stream: obspy.core.stream.Stream or None
        :return stream: stream object containing relevant waveforms

        Keyword Arguments
        ::
            str syn_pathname:
                Config.paths key to search for synthetic data. Defaults to
                'synthetics', but for the may need to be set to 'waveforms'
                in certain use-cases.
            str syn_unit:
                Optional argument to specify the letter used to identify the
                units of the synthetic data: For Specfem3D: ["d", "v", "a", "?"]
                'd' for displacement, 'v' for velocity,  'a' for acceleration.
                Wildcards okay. Defaults to '?'
            str syn_dir_template:
                Directory structure template to search for synthetic waveforms.
                Defaults to empty string
            str syn_fid_template:
                The naming template of synthetic waveforms defaults to
                "{net}.{sta}.*{cmp}.sem{syn_unit}"
        """
        syn_cfgpath = kwargs.get("syn_cfgpath", "synthetics")
        syn_unit = kwargs.get("syn_unit", "?")
        syn_dir_template = kwargs.get("syn_dir_template", "")
        syn_fid_template = kwargs.get("syn_fid_template",
                                      "{net}.{sta}.*{cmp}.sem{dva}"
                                      )

        if self.origintime is None:
            raise AttributeError("'origintime' must be specified")

        # Generate information necessary to search for data
        net, sta, loc, cha = code.split('.')

        # Ensure that the paths are a list so that iterating doesnt accidentally
        # try to iterate through a string.
        paths = self.config.paths[syn_cfgpath]
        if not isinstance(paths, list):
            paths = [paths]

        for path_ in paths:
            if not os.path.exists(path_):
                continue

            # Here the path is determined for search. If event_id is given,
            # the function will search for an event_id directory.
            full_path = os.path.join(path_, syn_dir_template, syn_fid_template)
            logger.debug(f"searching for synthetics: {full_path}")
            st = Stream()
            for filepath in glob.glob(full_path.format(
                    net=net, sta=sta, cmp=cha[2:], dva=syn_unit.lower())):
                try:
                    # Convert the ASCII file to a miniseed
                    st += read_sem(filepath, self.origintime)
                except UnicodeDecodeError:
                    # If the data file is for some reason already in miniseed
                    st += read(filepath)
                logger.info(f"retrieved synthetics locally:\n{filepath}")
            if len(st) > 0:
                st.merge()
                st.trim(starttime=self.origintime - self.config.start_pad,
                        endtime=self.origintime + self.config.end_pad
                        )
                return st
        else:
            return None

    def obs_waveform_fetch(self, code, **kwargs):
        """
        Mid-level internal fetching function for observation waveform data.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream or None
        :return: stream object containing relevant waveforms, or None
        """
        if self.ds:
            try:
                # Search the given ASDFDataSet first
                logger.info("searching ASDFDataSet for observations")
                return self.asdf_waveform_fetch(code,
                                                tag=self.config.observed_tag)
            except KeyError:
                pass
        logger.info("searching local filesystem for observations")
        if self.config.synthetics_only:
            return self.fetch_syn_by_dir(code, syn_cfgpath="waveforms",
                                         **kwargs)
        else:
            return self.fetch_obs_by_dir(code, **kwargs)

    def syn_waveform_fetch(self, code, **kwargs):
        """
        Mid-level internal fetching function for synthetic waveform data.

        .. note::
            Checks if synthetics are already saved into a ASDFDataSet first.
            If synthetics are freshly generated from Specfem3D, they should have
            been placed into folders separated by event id and model iteration.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream or None
        :return: stream object containing relevant waveforms, or None if no data
            is found
        """
        if self.ds:
            try:
                logger.info("searching ASDFDataSet for synthetics")
                return self.asdf_waveform_fetch(code,
                                                tag=self.config.synthetic_tag)
            except KeyError:
                pass
        logger.info("searching local filesystem for synthetics")
        return self.fetch_syn_by_dir(code, **kwargs)

    def event_fetch(self, event_id, **kwargs):
        """
        Mid-level internal fetching function for event information.
        Search ASDFDataSet for corresponding evemt, else look on disk.

        :type event_id: str
        :param event_id: Unique event identifier to search source file by.
            e.g., a New Zealand earthquake ID '2018p130600'. A prefix or suffix
            will be tacked onto this 
        :rtype event: obspy.core.event.Event or None
        :return event: event object if found, else None.
        """
        if self.ds:
            try:
                logger.info("searching ASDFDataSet for event info")
                return self.asdf_event_fetch()
            except (IndexError, AttributeError):
                pass
        logger.info("searching local filesystem for event info")
        return self.fetch_event_by_dir(event_id, **kwargs)

    def station_fetch(self, code, **kwargs):
        """
        Mid-level internal fetching function for station dataless information.
        Search ASDFDataSet for corresponding dataless, else look on disk.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.Inventory or None
        :return: inventory containing relevant network and stations, or None
            if no data is found
        """
        if self.ds:
            try:
                logger.info("searching ASDFDataSet for station info")
                return self.asdf_station_fetch(code)
            except (KeyError, AttributeError):
                pass
        logger.info("searching local filesystem for station info")
        return self.fetch_resp_by_dir(code, **kwargs)


class Gatherer(InternalFetcher, ExternalGetter):
    """
    A mid-level wrapper class used to fetch data internally and externally.
    Inherets internal and external data gathering from lower level
    InternalFetcher and ExternalGetter classes. 

    All saving to ASDFDataSet taken care of by the Gatherer class.
    """
    def __init__(self, config, ds=None, origintime=None):
        """
        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset for internal data searching and saving
        """
        self.ds = ds
        self.config = config
        self.origintime = origintime
        if self.config.client is not None:
            self.Client = Client(self.config.client)
        else:
            self.Client = None

    def gather_event(self, event_id=None, try_fm=True, **kwargs):
        """
        Gather an ObsPy Event object by searching disk then querying webservices

        .. note::
            Event info need only be retrieved once per Pyatoa workflow.

        :type try_fm: bool
        :param try_fm: try to find correspondig focal mechanism.
        :rtype: obspy.core.event.Event 
        :return: event retrieved either via internal or external methods
        :raises GathererNoDataException: if no event information is found.
        """
        logger.debug("gathering event")
        event = None

        # Attempt to gather event information internally
        event = self.event_fetch(event_id, **kwargs)
        # If no data internally, query FDSN
        if event is None and self.Client:
            event = self.event_get(event_id)
            # Append focal mechanism or moment tensor information, which is 
            # likely stored in a separate catalog
            if try_fm:
                event = append_focal_mechanism(event, 
                                               client=self.config.client)
        # If no event after internal/external checks, throw error
        if event is None:
            raise GathererNoDataException(f"no Event information found for "
                                          f"{self.config.event_id}")
        # Otherwise state success and grab important origin information
        else:        
            self.origintime = event.preferred_origin().time
        # Save event information to dataset if necessary
        if self.ds and self.config.save_to_ds:
            try:
                self.ds.add_quakeml(event)
                logger.debug(f"event QuakeML added to ASDFDataSet")
            # Trying to re-add an event to the ASDFDataSet throws ValueError
            except ValueError:
                pass
        return event

    def gather_station(self, code, **kwargs):
        """
        Gather station dataless information. Check disk then query webservices.
        Save station information to ASDFDataSet if requested.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        logger.info("gathering StationXML")
        inv = self.station_fetch(code, **kwargs)
        if inv is None:
            inv = self.station_get(code, **kwargs)
            if inv is None:
                raise GathererNoDataException(
                    f"no StationXML for {code} found"
                    )
        logger.info("matching StationXML found")
        if (self.ds is not None) and self.config.save_to_ds:
            # !!! This is a temp fix for PyASDF 0.6.1 where re-adding StationXML 
            # !!! that contains comments throws a TypeError. Issue #59
            try: 
                self.ds.add_stationxml(inv)
                logger.info("saved to ASDFDataSet")
            except TypeError:
                pass

        return inv

    def gather_observed(self, code, **kwargs):
        """
        Gather observed waveforms as ObsPy streams.
        Check disk, else query webservice. Save to ASDFDataSet if requested.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        logger.info("gathering observed waveforms")
        st_obs = self.obs_waveform_fetch(code, **kwargs)
        if st_obs is None:
            st_obs = self.obs_waveform_get(code)
            if st_obs is None:
                raise GathererNoDataException(
                    f"no observed waveforms for {code} found"
                    )
        logger.info("matching observed waveforms found")
        self._save_waveforms_to_dataset(st_obs, self.config.observed_tag)

        return st_obs

    def gather_synthetic(self, code, **kwargs):
        """
        Gather synthetic waveforms as ObsPy streams.
        Only possible to check ASDFDataSet and local filesystem.

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        :raises GathererNoDataException: if no synthetic data is found
        """
        logger.info("gathering synthetic waveforms")
        st_syn = self.syn_waveform_fetch(code, **kwargs)
        if st_syn is None:
            raise GathererNoDataException(f"no synthetic waveforms found "
                                          f"for {code}"
                                          )
        logger.info("matching synthetic waveforms found")
        self._save_waveforms_to_dataset(st_syn, self.config.synthetic_tag)

        return st_syn

    def gather_obs_multithread(self, codes, max_workers=None,
                               print_exception=False, **kwargs):
        """
        A multithreaded function that fetches all observed data (waveforms and
        StationXMLs) for a given event and store it to an ASDFDataSet.
        Multithreading is used to provide significant speed up for these request
        based tasks.

        :type codes: list of str
        :param codes: A list of station codes where station codes must be in the
            form NN.SSSS.LL.CCC (N=network, S=station, L=location, C=channel)
        :type max_workers: int
        :param max_workers: number of concurrent threads to use, passed to the
            ThreadPoolExecutor. If left as None, conurrent futures will
            automatically choose the system's number of cores.

        Keyword Arguments
        ::
            int return_count:
                if not None, determines how many data items must be collected
                for the station to be saved into the ASDFDataSet.
                e.g. StationXML and 3 component waveforms would equal 4 pieces
                of data, so a return_count == 4 means stations that do not
                return all components and metadata will not be saved to the
                dataset.
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed

        logger.info("mass gathering observation data")

        assert(self.ds is not None), \
            "Mass gathering requires a dataset `ds` for data storage"
        assert(self.Client is not None), \
            "Mass gathering requires a Client for data queries"
        assert(self.origintime is not None), \
            "Mass gathering requires an origintime for data queries"

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(self._obs_get_multithread, code,
                                       **kwargs):
                           code for code in codes
                       }
            for future in as_completed(futures):
                code = futures[future]
                try:
                    status = future.result()
                except Exception as e:
                    print(f"{code} exception: {e}\n")
                    if print_exception:
                        traceback.print_exc()
                else:
                    print(f"{code} data count: {status}")

        #  GCMT

    def _save_waveforms_to_dataset(self, st, tag):
        """
        Save waveformsm to the ASDFDataSet with a simple check for existence
        of dataset and save parameter. Passes if waveforms already exist while
        ignoring the PyASDF warning that gets thrown if waveforms exist.

        :type st: obspy.core.stream.Stream
        :param st: Stream object to be saved into the dataset
        :type tag: str
        :param tag: unique identifier to save the waveforms under
        """
        if (self.ds is not None) and self.config.save_to_ds:
            # Catch ASDFWarning that occurs when data already exists
            with warnings.catch_warnings():
                warnings.filterwarnings("error")
                try:
                    self.ds.add_waveforms(waveform=st, tag=tag)
                    logger.info(f"saved to ASDFDataSet with tag '{tag}'")
                except ASDFWarning:
                    pass


def append_focal_mechanism(event, client=None, overwrite=False):
    """
    Attempt to find focal mechanism information with a given Event object.

    .. note::
        FDSN fetched events are devoid of a few bits of information that are
        useful for our applications, e.g. moment tensor, focal mechanisms.
        This function will perform the conversions and append the necessary
        information to the event located in the dataset.

    :type event: obspy.core.event.Event
    :param event: Event object to append a focal mechanism to.
    :type overwrite: bool
    :param overwrite: If the event already has a focal mechanism, this will
        overwrite that focal mechanism
    :raises TypeError: if event is not provided as an obspy.core.event.Event
    """
    if isinstance(event, Event):
        event_id = format_event_name(event)

        # If the event already has a focal mechanism attribute, don't gather
        if hasattr(event, 'focal_mechanisms') and \
                event.focal_mechanisms and not overwrite:
            return event
        if client and client.upper() == "GEONET":
            # Query GeoNet moment tensor catalog if using GeoNet catalog
            from pyatoa.plugins.new_zealand.gather import  geonet_mt
            event, _ = geonet_mt(event_id=event_id, event=event, units="nm")
            logger.info("GeoNet moment tensor appended to Event")
        else:
            try:
                # Try to query GCMT web-based catalog for matching event
                event = get_gcmt_moment_tensor(
                    origintime=event.preferred_origin().time,
                    magnitude=event.preferred_magnitude().mag
                )
            except FileNotFoundError:
                logger.info("no GCMT moment tensor for event found")
    else:
        raise TypeError("'event' must be an ObsPy Event object")

    return event


def get_gcmt_moment_tensor(origintime, magnitude, time_wiggle_sec=120,
                           magnitude_wiggle=0.5):
    """
    Query GCMT moment tensor catalog for moment tensor components

    :type origintime: UTCDateTime or str
    :param origintime: event origin time
    :type magnitude: float
    :param magnitude: centroid moment magnitude for event lookup
    :type time_wiggle_sec: int
    :param time_wiggle_sec: padding on catalog filtering criteria realted to
        event origin time
    :type magnitude_wiggle: float
    :param magnitude_wiggle: padding on catalog filter for magnitude
    :rtype: obspy.core.event.Event
    :return: event object for given earthquake
    """
    from urllib.error import HTTPError
    from obspy import UTCDateTime, read_events

    if not isinstance(origintime, UTCDateTime):
        #  GCMT
        origintime = UTCDateTime(origintime)

    # Determine filename using datetime properties
    month = origintime.strftime('%b').lower()  # e.g. 'jul'
    year_short = origintime.strftime('%y')  # e.g. '19'
    year_long = origintime.strftime('%Y')  # e.g. '2019'

    fid = f"{month}{year_short}.ndk"
    logger.info("querying GCMT database for moment tensor")
    try:
        cat = read_events(
            "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
            f"catalog/NEW_MONTHLY/{year_long}/{fid}"
        )
    except HTTPError:
        cat = read_events(
            "http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
            "catalog/NEW_QUICK/qcmt.ndk"
        )

    # GCMT catalogs contain all events for a span of time
    # filter catalogs using ObsPy to find events with our specifications.
    # Magnitudes and origintimes are not always in agreement between agencies
    # So allow for some wiggle room
    cat_filt = cat.filter(f"time > {str(origintime - time_wiggle_sec)}",
                          f"time < {str(origintime + time_wiggle_sec)}",
                          f"magnitude >= {magnitude - magnitude_wiggle}",
                          f"magnitude <= {magnitude + magnitude_wiggle}",
                          )
    # Filtering may remove all events from catalog, return multiple events, or
    # may return the event of choice
    if not len(cat_filt):
        logger.info(f"no GCMT event found for {origintime} and M{magnitude}")
        raise FileNotFoundError("No events found")
    elif len(cat_filt) > 1:
        logger.info(f"multiple events found for {origintime} and M{magnitude}")
        print(f"{len(cat_filt)} events found, choosing first")
        return cat_filt[0]
    else:
        logger.info("GCMT event found matching criteria")
        return cat_filt[0]


def get_usgs_moment_tensor(origintime, magnitude, time_wiggle_sec=120,
                           **kwargs):
    """
    Query FDSN webservices USGS client for moment tensors

    Kwargs passed to Client.get_events() for additional event constraint pars

    :type origintime: UTCDateTime or str
    :param origintime: event origin time
    :type magnitude: float
    :param magnitude: centroid moment magnitude for event lookup
    :type time_wiggle_sec: int
    :param time_wiggle_sec: padding on catalog filtering criteria realted to
        event origin time
    :type magnitude_wiggle: float
    :param magnitude_wiggle: padding on catalog filter for magnitude
    :rtype: obspy.core.event.Event
    :return: event object for given earthquake
    """
    c = Client("USGS")
    cat = c.get_events(starttime=origintime - time_wiggle_sec, 
                       endtime=origintime + time_wiggle_sec, 
                       includeallorigins=True, **kwargs
                       )
    return cat
