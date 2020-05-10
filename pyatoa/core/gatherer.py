#!/usr/bin/env python
"""
Mid and Low level data gathering classes to retrieve data from local filesystems
or to query data from FDSN webservices via ObsPy.

Fetch: used to denote searching local filesystems for data (internal)
Get: used to denote querying FDSN webservices via ObsPy (external), naming
     *convention from the obspy.fdsn.client.Client function names

Gatherer directly called by the Manager class and shouldn't need to be called
by the User unless for bespoke data gathering functionality.
"""
import os
import glob
import warnings

from pyasdf import ASDFWarning
from obspy.core.event import Event
from obspy.clients.fdsn import Client
from obspy import Stream, read, read_inventory
from obspy.clients.fdsn.header import FDSNNoDataException, FDSNException

from pyatoa import logger
from pyatoa.utils.read import read_ascii
from pyatoa.utils.calculate import overlapping_days
from pyatoa.utils.srcrcv import merge_inventories


class ExternalGetter:
    """
    Low-level gathering classs to retrieve data via FDSN webservices.
    Calls made through ObsPy. Functionality is inhereted by the Gatherer class.
    """
    def event_get(self):
        """
        Return event information parameters pertaining to a given event id
        if an event id is given

        :rtype event: obspy.core.event.Event
        :return event: event object
        """
        event = None
        logger.debug(f"fetching event from {self.config.client}")
        if self.config.event_id is not None:
            # Get events via event id, only available through certain clients
            try:
                event = self.Client.get_events(eventid=self.config.event_id)[0]
                self.origintime = event.origins[0].time
            except FDSNNoDataException:
                logger.warning(f"no event found for {self.config.event_id} "
                               f"from {self.config.client}")
                event = None
        if self.origintime and event is None:
            # If getting by event id doesn't work, try based on origintime
            try:
                event = self.Client.get_events(starttime=self.origintime,
                                               endtime=self.origintime)
            except FDSNNoDataException:
                logger.warning(
                    f"no event found for origin time {self.origintime}"
                    f"from {self.config.client}"
                )
            if event is not None and len(event) > 1:
                logger.warning(f"{len(event)} events found, expected only 1,"
                               f"returning first entry but manual revision "
                               f"may be required."
                               )
                event = event[0]
        return event

    def station_get(self, station_code, level="response"):
        """
        Call for ObsPy FDSN client to download station dataless information.
        Defaults to retrieving response information.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type level: str
        :param level: level argument to be passed to obspy
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        logger.debug(f"fetching station from {self.config.client}")
        net, sta, loc, cha = station_code.split('.')
        return self.Client.get_stations(
            network=net, station=sta, location=loc, channel=cha,
            starttime=self.origintime - self.config.start_pad,
            endtime=self.origintime + self.config.end_pad, level=level
        )

    def waveform_get(self, station_code):
        """
        Call for ObsPy FDSN client to download waveform data.

        Note:
            ObsPy sometimes returns traces with varying sample lengths,
            so we use a 10 second cushion on start and end time and trim after
            retrieval to make sure traces are the same length.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stream
        """
        logger.debug(f"fetching observations from {self.config.client}")

        net, sta, loc, cha = station_code.split('.')
        st = self.Client.get_waveforms(
            network=net, station=sta, location=loc, channel=cha,
            starttime=self.origintime - (self.config.start_pad + 10),
            endtime=self.origintime + (self.config.end_pad + 10)
        )
        # Sometimes FDSN queries return improperly cut start and end times, so
        # we retrieve +/-10 seconds and then cut down
        st.trim(starttime=self.origintime - self.config.start_pad,
                endtime=self.origintime + self.config.end_pad)

        logger.debug(f"stream got external {station_code}")
        return st


class InternalFetcher:
    """
    Low-level data gatherer to search for data on disk.
    
    Initially looks through User provided ASDFDataSets, if not available or no
    data locatable, search through the local directories.

    Defaults directory structure follows SEED convention but can be overwritten.

    Class attributes inhereted by the Gatherer class.
    """
    def asdf_event_fetch(self):
        """
        Return Event information from ASDFDataSet.

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
        Return station dataless from ASDFDataSet based on station code.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.network.Network
        :return: network containing relevant station information
        """
        net, sta, loc, cha = station_code.split(".")
        return self.ds.waveforms[f"{net}_{sta}"].StationXML.select(channel=cha)

    def asdf_waveform_fetch(self, station_code, tag):
        """
        Return waveforms as Stream objects based from ASDFDataSet.

        Note:
            -Allows for wildcard selection of component (? or *)
            -Selects by component because synthetic channel naming may differ
            from observation channels.
            -Component is assumed to be the last index in the channel, following
            SEED convention.

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
        return self.ds.waveforms[f"{net}_{sta}"][tag].select(component=cha[-1])

    def fetch_resp_by_dir(self, station_code, paths_to_responses=None,
                          dir_structure="{sta}.{net}",
                          file_template="RESP.{net}.{sta}.{loc}.{cha}"):
        """
        Fetch station dataless via directory structure on disk.

        Will search through all paths given until corresponding inventories are
        found else raises a FileNotFoundError.

        Default path naming follows SEED convention, that is:
            path/to/dataless/{NET}.{STA}/RESP.{NET}.{STA}.{LOC}.{CHA}
            e.g. path/to/dataless/NZ.BFZ/RESP.NZ.BFZ.10.HHZ

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
            paths_to_responses = self.config.cfgpaths["responses"]

        inv = None
        net, sta, loc, cha = station_code.split('.')
        for path_ in paths_to_responses:
            if not os.path.exists(path_):
                continue
            # Inventory() requires some positional arguements, use None to skip
            fid = os.path.join(path_, dir_structure, file_template).format(
                net=net, sta=sta, cha=cha, loc=loc)
            for filepath in glob.glob(fid):
                if inv is None:
                    # The first inventory becomes the main inv to return
                    inv = read_inventory(filepath)
                else:
                    # All other inventories are appended to the original
                    inv_append = read_inventory(filepath)
                    # Merge inventories to remove repeated networks
                    inv = merge_inventories(inv, inv_append)

                logger.debug(f"response found at {filepath}")
        if inv is None:
            logger.debug("no response found for given paths")
            raise FileNotFoundError()

        return inv

    def fetch_obs_by_dir(
            self, station_code, paths_to_waveforms=None,
            dir_structure="{year}/{net}/{sta}/{cha}*",
            file_template="{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}"):
        """
        Fetch observation waveforms via directory structure on disk.

        Default waveform directory structure assumed to follow SEED convention.
        That is:
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
        # If waveforms contain midnight, multiple files need to be read
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
                    logger.debug(f"stream fetched from directory {filepath}")
            if len(st) > 0:
                # Take care of gaps in data by converting to masked data
                st.merge()
                st.trim(starttime=self.origintime-self.config.start_pad,
                        endtime=self.origintime+self.config.end_pad
                        )
                return st
        else:
            logger.debug(
                f"no waveforms found for {station_code} for given directories"
            )
            raise FileNotFoundError()

    def fetch_syn_by_dir(self, station_code, event_id="", pathname="synthetics",
                         specfem_fid_template="{net}.{sta}.*{cmp}.sem{dva}"):
        """
        Fetch synthetic waveforms from Specfem3D via directory structure on
        disk, if necessary convert native ASCII format to Stream object.

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
        :type pathname: str
        :param pathname: the key in cfgpaths dictionary to search for data
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
        for path_ in self.config.cfgpaths[pathname]:
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
                    st += read_ascii(filepath, self.origintime)
                except UnicodeDecodeError:
                    # If the data file is for some reason already in miniseed
                    st += read(filepath)
                logger.debug(
                    f"stream fetched by event {os.path.basename(filepath)}"
                )
            if len(st) > 0:
                st.merge()
                st.trim(starttime=self.origintime - self.config.start_pad,
                        endtime=self.origintime + self.config.end_pad
                        )
                return st
        else:
            # This needs to be a warn because there are no more checks for
            # synthetic data after this
            logger.warn(
                f"no synthetic waveforms for {station_code} found for event"
            )
            raise FileNotFoundError()

    def obs_waveform_fetch(self, station_code):
        """
        Low-level internal fetching function for observation waveform data.
        Will return a FileNotFoundError if no internal data is found.

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
                # Search the given asdf dataset first
                logger.debug("fetching obs internal asdf")
                return self.asdf_waveform_fetch(station_code,
                                                tag=self.config.observed_tag)
            except KeyError:
                logger.debug("obs internal asdf failed, fetching by dir")
                # If ADSFDataSet does not contain data, search internal
                # If this is a synthetic-synthetic example, search the waveforms
                # directory, using the synthetic search function
                if self.config.synthetics_only:
                    logger.debug("synthetics only case, searching for syn data")
                    st_obs = self.fetch_syn_by_dir(station_code,
                                                   pathname="waveforms"
                                                   )
                # Else look for observations using SEED convention
                else:
                    st_obs = self.fetch_obs_by_dir(station_code)
                return st_obs
        else:
            return self.fetch_obs_by_dir(station_code)

    def syn_waveform_fetch(self, station_code):
        """
        Low-level internal fetching function for synthetic waveform data.

        Checks if synthetics are already saved into a ASDFDataSet first.
        If synthetics are freshly generated from Specfem3D, they should have
        been placed into folders separated by event id and model iteration.

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
                return st_syn
        else:
            return self.fetch_syn_by_dir(station_code)

    def station_fetch(self, station_code):
        """
        Low-level internal fetching function for station dataless information.

        Search ASDFDataSet for corresponding dataless, else look on disk.
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
                return inv
        else:
            return self.fetch_resp_by_dir(station_code)


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
        self.Client = Client(self.config.client)

    def gather_event(self, append_focal_mechanism=True):
        """
        Gather an ObsPy Event object by searching disk then querying webservices

        Event need only be retrieved once per Pyatoa workflow.

        :type append_focal_mechanism: bool
        :param append_focal_mechanism: try to find correspondig focal mechanism.
        :rtype: obspy.core.event.Event
        :return: event retrieved either via internal or external methods
        """
        # If dataset is given, search for event
        if self.ds:
            try:
                event = self.asdf_event_fetch()
            except (AttributeError, IndexError):
                event = self.event_get()
                if append_focal_mechanism:
                    event = self.append_focal_mechanism(event)
                if self.config.save_to_ds:
                    self.ds.add_quakeml(event)
                    logger.debug(
                        f"event retrieved from client {self.config.client} "
                        f"added to pyasdf dataset"
                    )
                else:
                    logger.debug("event not saved")
        # Else, query FDSN for event information
        else:
            event = self.event_get()
            if append_focal_mechanism:
                event = self.append_focal_mechanism(event)
            logger.debug(f"event retrieved from client {self.config.client}")

        return event

    def append_focal_mechanism(self, event, overwrite=False):
        """
        Attempt to find focal mechanism information with a given Event object.

        Focal mechanisms are unforunately not something that is standard to
        store, as calculations differ between regions and agencies.

        That means this function is not always bulletproof.

        FDSN fetched events are devoid of a few bits of information that are
        useful for our applications, e.g. moment tensor, focal mechanisms.
        This function will perform the conversions and append the necessary
        information to the event located in the dataset.

        :type event: obspy.core.event.Event
        :param event: Event object to append a focal mechanism to.
        :type overwrite: bool
        :param overwrite: If the event already has a focal mechanism, this will
            overwrite that focal mechanism
        """
        if isinstance(event, Event):
            # If the event already has a focal mechanism attribute, don't gather
            if hasattr(event, 'focal_mechanisms') and \
                    event.focal_mechanisms and not overwrite:
                return event
            if self.config.client.upper() == "GEONET":
                # Query GeoNet moment tensor catalog
                from pyatoa.plugins.new_zealand.gather import \
                                                    geonet_focal_mechanism
                event, _ = geonet_focal_mechanism(event_id=self.config.event_id,
                                                  event=event
                                                  )
                logger.info("appending GeoNet moment tensor to event")
            else:
                try:
                    # Try to query GCMT web-based catalog for matching event
                    event = get_gcmt_moment_tensor(
                        origintime=event.preferred_origin().time,
                        magnitude=event.preferred_magnitude().mag
                    )
                except FileNotFoundError:
                    logger.info("No GCMT event found")
        else:
            raise NotImplementedError("`event` must be an obspy Event object")

        return event

    def gather_station(self, station_code):
        """
        Gather station dataless information. Check disk then query webservices.
        Save station information to ASDFDataSet if requested.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        try:
            return self.station_fetch(station_code)
        except FileNotFoundError:
            logger.debug(
                "internal station information not found, searching ext.")
            inv = self.station_get(station_code)
            if self.ds is not None and self.config.save_to_ds:
                self.ds.add_stationxml(inv)
            else:
                logger.debug("station information is not being saved")
            return inv

    def gather_observed(self, station_code):
        """
        Gather observed waveforms as ObsPy streams.
        Check disk, else query webservice. Save to ASDFDataSet if requested.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        try:
            st_obs = self.obs_waveform_fetch(station_code)
        except FileNotFoundError:
            logger.debug("no internal obs data, searching external")
            try:
                st_obs = self.waveform_get(station_code)
            # Catch all FDSN Exceptions
            except FDSNException:
                logger.warning("no obs stream can be returned")
                st_obs = None
        # Save waveforms to ASDFDataSet
        if (self.ds is not None) and self.config.save_to_ds and (
                st_obs is not None):
            # Ignore the ASDFWarning that occurs if data exists in dataset
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ASDFWarning)
                self.ds.add_waveforms(waveform=st_obs,
                                      tag=self.config.observed_tag)
        else:
            logger.debug("observed waveforms not being saved")

        return st_obs

    def gather_synthetic(self, station_code):
        """
        Gather synthetic waveforms as ObsPy streams.
        Check disk, else query webservice. Save to ASDFDataSet if requested.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        try:
            st_syn = self.syn_waveform_fetch(station_code)
            if (self.ds is not None) and self.config.save_to_ds and (
                st_syn is not None):
                # Ignore the ASDFWarning that occurs if data exists in dataset
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", ASDFWarning)
                    self.ds.add_waveforms(waveform=st_syn,
                                          tag=self.config.synthetic_tag)
            return st_syn
        except FileNotFoundError:
            return None


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
    :rtype event: obspy.core.event.Event
    :return event: event object for given earthquake
    """
    from urllib.error import HTTPError
    from obspy import UTCDateTime, read_events

    if not isinstance(origintime, UTCDateTime):
        datetime = UTCDateTime(origintime)

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
    # So allow fro some wiggle room
    cat_filt = cat.filter(f"time > {str(origintime - time_wiggle_sec)}",
                          f"time < {str(origintime + time_wiggle_sec)}",
                          f"magnitude >= {magnitude - magnitude_wiggle}",
                          f"magnitude <= {magnitude + magnitude_wiggle}",
                          )
    # Filtering may remove all events from catalog, return multiple events, or
    # may return the event of choice
    if not len(cat_filt):
        logger.info(f"no gcmt event found for {datetime} and M{magnitude}")
        raise FileNotFoundError("No events found")
    elif len(cat_filt) > 1:
        logger.info(f"multiple events found for {datetime} and M{magnitude}")
        print(f"{len(cat_filt)} events found, choosing first")
        return cat_filt[0]
    else:
        logger.info("gcmt event found matching criteria")
        return cat_filt[0]

