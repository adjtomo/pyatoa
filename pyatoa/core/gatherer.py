#!/usr/bin/env python
"""
Mid and Low level data gathering classes to retrieve data from local filesystems
either via an ASDFDataSet or through a pre-defined data directory structure.

Gatherer directly called by the Manager class and shouldn't need to be called
by the User unless for bespoke data gathering functionality.
"""
import os
import glob
import warnings

from pyasdf import ASDFWarning
from pysep.utils.io import read_sem, read_specfem2d_source, read_forcesolution
from obspy import Stream, read, read_inventory, read_events

from pyatoa import logger
from pyatoa.utils.form import format_event_name
from pyatoa.utils.calculate import overlapping_days
from pyatoa.utils.srcrcv import merge_inventories


class GathererNoDataException(Exception):
    """
    Custom exception to be thrown generally to Manager class in the case that
    gathering of any data fails.
    """
    pass


class Gatherer:
    """
    A mid-level data gathering class used to get data internally and externally.
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

    def gather_event(self, event_id=None, **kwargs):
        """
        Gather an ObsPy Event object by searching disk
        Event info need only be retrieved once per Pyatoa workflow.

        :type event_id: str
        :param event_id: a unique event idenfitier to search and tag event info
        :rtype: obspy.core.event.Event
        :return: event retrieved either via internal or external methods
        :raises GathererNoDataException: if no event information is found.
        """
        logger.info("gathering event QuakeML")
        # Try 1: fetch from ASDFDataSet
        if self.ds:
            logger.debug("searching ASDFDataSet for event info")
            event = self.fetch_event_from_dataset()
        else:
            event = None
        # Try 2: look at local filesystem
        if event is None:
            logger.debug("searching local filesystem for QuakeML")
            event = self.fetch_event_by_dir(event_id, **kwargs)
        # Abort if none of the three attempts returned successfully
        if event is None:
            raise GathererNoDataException(f"no QuakeML found for "
                                          f"{self.config.event_id}")

        # Overwrite origintime with the catalog value to be consistent
        self.origintime = event.preferred_origin().time
        logger.debug(f"event origin time is set to: {self.origintime}")

        # Save event information to dataset if necessary
        if self.ds and self.config.save_to_ds:
            try:
                self.ds.add_quakeml(event)
                logger.debug(f"event QuakeML added to ASDFDataSet")
            # Trying to re-add an event to the ASDFDataSet throws ValueError
            except ValueError as e:
                logger.debug(f"event QuakeML was not added to ASDFDataSet: {e}")
                pass

        logger.info(f"matching event found: '{format_event_name(event)}'")
        return event

    def gather_station(self, code, **kwargs):
        """
        Gather StationXML information from disk

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        logger.info(f"gathering StationXML for: {code}")

        # Try 1: fetch from ASDFDataSet
        if self.ds:
            logger.debug("searching ASDFDataSet for StationXML")
            inv = self.fetch_inv_from_dataset(code)
        else:
            inv = None
        # Try 2: fetch from local filesystem
        if inv is None:
            logger.debug("searching local filesystem for StationXML")
            inv = self.fetch_inv_by_dir(code, **kwargs)
        # Abort if none of the three attempts returned successfully
        if inv:
            logger.info(f"matching StationXML found: {code}")
            if (self.ds is not None) and self.config.save_to_ds:
                # !!! This is a temp fix for PyASDF 0.6.1 where re-adding
                # !!! StationXML that contains comments throws a TypeError.
                # !!! Issue #59
                try:
                    self.ds.add_stationxml(inv)
                    logger.debug("saved inventory to ASDFDataSet")
                except TypeError:
                    pass
        else:
            raise GathererNoDataException(f"no StationXML for {code} found")

        return inv

    def gather_observed(self, code, **kwargs):
        """
        Gather observed waveforms from disk as ObsPy streams

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype: obspy.core.stream.Stream
        :return: stream object containing relevant waveforms
        """
        logger.info(f"gathering observed waveforms: {code}")

        # Try 1: fetch from ASDFDataSet
        if self.ds:
            logger.debug("searching ASDFDataSet for observations")
            st_obs = self.fetch_waveform_from_dataset(
                code=code, tag=self.config.observed_tag
            )
        else:
            st_obs = None
        # Try 2: fetch from local file system, different approaches if we're
        # looking for data or synthetic "data"
        if st_obs is None:
            logger.debug("searching local filesystem for observations")
            if self.config.synthetics_only:
                st_obs = self.fetch_synthetic_by_dir(code,
                                                     syn_cfgpath="waveforms",
                                                     **kwargs)
            else:
                st_obs = self.fetch_observed_by_dir(code, **kwargs)
        # Abort if none of the three attempts returned successfully
        if st_obs:
            logger.info(f"matching observed waveforms found: {code}")
            self.save_waveforms_to_dataset(st_obs, self.config.observed_tag)
        else:
            raise GathererNoDataException(f"no observed waveforms found "
                                          f"for: {code}")
        return st_obs

    def gather_synthetic(self, code, **kwargs):
        """
        Gather synthetic waveforms as ObsPy streams.

        Only possible to check ASDFDataSet and local filesystem, cannot gather
        synthetics from webservice.

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
        logger.info(f"gathering synthetic waveforms: {code}")

        # Try 1: fetch from ASDFDataSet
        if self.ds:
            logger.debug("searching ASDFDataSet for synthetics")
            st_syn = self.fetch_waveform_from_dataset(
                code=code, tag=self.config.synthetic_tag
            )
        else:
            st_syn = None
        # Try 2: fetch from local directories
        if st_syn is None:
            logger.debug("searching local filesystem for synthetics")
            st_syn = self.fetch_synthetic_by_dir(code, **kwargs)
        # Abort if none of the attempts returned successfully
        if st_syn:
            logger.debug(f"matching synthetic waveforms found: {code}")
            self.save_waveforms_to_dataset(st_syn, self.config.synthetic_tag)
        else:
            raise GathererNoDataException(f"no synthetic waveforms found "
                                          f"for: {code}"
                                          )
        return st_syn

    def fetch_event_from_dataset(self):
        """
        Return Event information from ASDFDataSet.

        .. note::
            Assumes that the ASDF Dataset will only contain one event, which is
            dictated by the structure of Pyatoa.

        :rtype event: obspy.core.event.Event
        :return event: event object
        :raises AttributeError: if no event attribute found in ASDFDataSet
        :raises IndexError: if event attribute found but no events
        """
        try:
            event = self.ds.events[0]
        except (IndexError, AttributeError):
            logger.debug(f"no matching event found in dataset")
            event = None
        return event

    def fetch_inv_from_dataset(self, code):
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
        try:
            inv = self.ds.waveforms[f"{net}_{sta}"].StationXML
            inv = inv.select(channel=cha)
        except (KeyError, AttributeError):
            logger.debug(f"no matching inventory in dataset: {net}_{sta}")
            inv = None
        return inv

    def fetch_waveform_from_dataset(self, code, tag):
        """
        Return waveforms as Stream objects from ASDFDataSet.

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
        :return: waveform contained in a stream, or None if no matching value
        """
        net, sta, loc, cha = code.split(".")
        comp = cha[-1]
        try:
            st = self.ds.waveforms[f"{net}_{sta}"][tag].select(component=comp)
        except KeyError:
            logger.debug(f"no matching waveform in dataset: {net}_{sta}_{tag}")
            st = None
        return st

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
                logger.debug(f"event search path does not exist: {path_}")
                continue
            # Search for available event files
            fid = os.path.join(path_, f"{prefix}{event_id}{suffix}")
            for filepath in glob.glob(fid):
                logger.debug(f"searching for event data: {filepath}")
                if os.path.exists(filepath):
                    try:
                        # Allow input of various types of source files
                        if "SOURCE" in prefix.upper():
                            logger.debug(f"found SPECFEM2D SOURCE: {filepath}")
                            cat = [read_specfem2d_source(filepath)]
                        elif "FORCESOLUTION" in prefix.upper():
                            logger.debug(f"found FORCESOLUTION: {filepath}")
                            cat = [read_forcesolution(filepath)]
                        # ObsPy can handle QuakeML and CMTSOLUTION
                        else:
                            logger.debug(f"found event file: {filepath}")
                            cat = read_events(filepath, format=format_)

                        if len(cat) != 1:
                            logger.warning(
                                f"{filepath} event file contains {len(cat)} "
                                f"events, returning 0th index")
                        event = cat[0]
                        break
                    except Exception as e:
                        logger.warning(f"{filepath} event file read error: {e}")

        if event is not None:
            logger.info(f"retrieved local file: {filepath}")
        else:
            logger.debug(f"no local event file found")

        return event

    def fetch_inv_by_dir(self, code, resp_dir_template="{sta}.{net}",
                         resp_fid_template="RESP.{net}.{sta}.{loc}.{cha}",
                         **kwargs):
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
        :type resp_dir_template: str
        :param resp_dir_template: Directory structure template to search for
            response files. By default follows the SEED convention:
            'path/to/RESPONSE/{sta}.{net}/'
        :type resp_fid_template: str
        :param resp_fid_template: Response file naming template to search for
            station dataless. By default, follows the SEED convention:
            'RESP.{net}.{sta}.{loc}.{cha}'
        :rtype inv: obspy.core.inventory.Inventory or None
        :return inv: inventory containing relevant network and stations
        """
        inv = None
        net, sta, loc, cha = code.split(".")

        # Ensure that the paths are a list so that iterating doesnt accidentally
        # try to iterate through a string.
        paths = self.config.paths["responses"]
        if not isinstance(paths, list):
            paths = [paths]

        for path_ in paths:
            if not os.path.exists(path_):
                logger.debug(f"StationXML search path does not exist: {path_}")
                continue
            # Attempting to instantiate an empty Inventory requires some
            # positional arguements we dont have, so don't do that
            fid = os.path.join(path_, resp_dir_template, resp_fid_template)
            fid = fid.format(net=net, sta=sta, cha=cha, loc=loc)
            logger.debug(f"searching for StationXML: {fid}")

            for filepath in glob.glob(fid):
                if inv is None:
                    # The first inventory becomes the main inv to return
                    inv = read_inventory(filepath)
                else:
                    # All other inventories are appended to the original
                    inv_append = read_inventory(filepath)
                    # Merge inventories to remove repeated networks
                    inv = merge_inventories(inv, inv_append)
                logger.info(f"retrieved StationXML locally: {filepath}")

        return inv

    def fetch_observed_by_dir(
            self, code,  obs_dir_template="{year}/{net}/{sta}/{cha}",
            obs_fid_template="{net}.{sta}.{loc}.{cha}.{year}.{jday:0>3}",
            **kwargs):
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
        :type obs_dir_template: str
        :param obs_dir_template: directory structure to search for observation
            data. Follows the SEED convention:
            'path/to/obs_data/{year}/{net}/{sta}/{cha}'
        :type obs_fid_template: str
        :param obs_fid_template: File naming template to search for observation
            data. Follows the SEED convention:
            '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
        :rtype stream: obspy.core.stream.Stream or None
        :return stream: stream object containing relevant waveforms, else None
        """
        if self.origintime is None:
            raise AttributeError("`origintime` must be specified")

        net, sta, loc, cha = code.split('.')
        # If waveforms contain midnight, multiple files need to be read.
        # Checks one hour before and after origintime
        jdays = overlapping_days(origin_time=self.origintime, start_pad=3600,
                                 end_pad=3600)

        # Ensure that the paths are a list so that iterating doesnt accidentally
        # try to iterate through a string.
        paths = self.config.paths["waveforms"]
        if not isinstance(paths, list):
            paths = [paths]

        st = Stream()
        for path_ in paths:
            if not os.path.exists(path_):
                logger.debug(f"waveform search path does not exist: {path_}")
                continue
            full_path = os.path.join(path_, obs_dir_template, obs_fid_template)
            pathlist = []
            for jday in jdays:
                pathlist.append(full_path.format(net=net, sta=sta, cha=cha,
                                                 loc=loc, jday=jday,
                                                 year=self.origintime.year)
                                )
            for fid in pathlist:
                logger.debug(f"searching for observations: {fid}")
                for filepath in glob.glob(fid):
                    st += read(filepath)
                    logger.info(f"retrieved observations locally: {filepath}")
            break
        # Take care of gaps in data by converting to masked data
        if len(st) > 0:
            st.merge()

        # If empty stream either due to no data or trimming removes all data,
        # we will return None
        if len(st) == 0:
            logger.warning(f"no matching observed waveforms found: {code}")
            st = None

        return st

    def fetch_synthetic_by_dir(self, code, syn_cfgpath="synthetics",
                               syn_unit="?", syn_dir_template="",
                               syn_fid_template="{net}.{sta}.*{cmp}.sem{dva}*",
                               **kwargs):
        """
        Fetch synthetic waveforms from Specfem3D via directory structure on
        disk, if necessary convert native ASCII format to Stream object.

        .. note::
            By default, synthetics will be searched for with the following path
            config.paths[syn_cfgpath]/syn_dir_template/syn_fid_template.format()

        :type code: str
        :param code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type syn_cfgpath: str
        :param syn_cfgpath: Config.paths key to search for synthetic data.
            Defaults to 'synthetics', but for the may need to be set to
            'waveforms' in certain use-cases.
        :type syn_unit: str
        :param syn_unit: Optional argument to specify the letter used to
            identify the units of the synthetic data: For Specfem3D:
            ["d", "v", "a", "?"] 'd' for displacement, 'v' for velocity,
            'a' for acceleration. Wildcards okay. Defaults to '?'
        :type syn_dir_template: str
        :param syn_dir_template: Directory structure template to search for
            synthetic waveforms. Defaults to empty string
        :type syn_fid_template: str
        :param syn_fid_template: The naming template of synthetic waveforms
            defaults to "{net}.{sta}.*{cmp}.sem{syn_unit}"
        :rtype stream: obspy.core.stream.Stream or None
        :return stream: stream object containing relevant waveforms
        """
        if self.origintime is None:
            raise AttributeError("`origintime` must be specified")

        # Generate information necessary to search for data
        net, sta, loc, cha = code.split(".")

        # Ensure that the paths are a list so that iterating doesnt accidentally
        # try to iterate through a string.
        paths = self.config.paths[syn_cfgpath]
        if not isinstance(paths, list):
            paths = [paths]

        st = Stream()
        for path_ in paths:
            if not os.path.exists(path_):
                logger.debug(f"synthetic search path does not exist: {path_}")
                continue
            # Expand the full path for searching for synthetics
            full_path = os.path.join(path_, syn_dir_template, syn_fid_template)
            filenames = glob.glob(full_path.format(net=net, sta=sta,
                                                   cmp=cha[2:],
                                                   dva=syn_unit.lower())
                                       )
            logger.debug(f"found {len(filenames)} synthetics: {full_path}")
            if filenames:
                for filename in filenames:
                    try:
                        # Convert the ASCII file to a miniseed
                        st += read_sem(filename, self.origintime)
                    except UnicodeDecodeError:
                        # If the data file is already in miniseed
                        st += read(filename)
                    logger.info(f"retrieved synthetics locally: {filename}")
            else:
                continue
        # Take care of gaps in data by converting to masked data
        if len(st) > 0:
            st.merge()
        # If empty stream either due to no data or trimming removes all data,
        # we will return None
        if len(st) == 0:
            logger.warning(f"no matching synthetic data found: {code}")
            st = None

        return st

    def save_waveforms_to_dataset(self, st, tag):
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
                    logger.debug(f"saved waveform to ASDFDataSet as: '{tag}'")
                except ASDFWarning:
                    pass
