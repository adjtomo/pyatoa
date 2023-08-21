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


def read_events_plus(fid, fmt, **kwargs):
    """
    Given a path `fid`, read an event/source file in as an ObsPy Event object.
    Wrapper for ObsPy's `read_events` that provides additional support for
    SPECFEM-specific source files.

    :type fid: str
    :param fid: full path to the event file to be read
    :type fmt: str
    :param fmt: Expected format of the file to read, available are 'SOURCE'
        (SPECFEM2D SOURCE file), 'FORCESOLUTION' (SPECFEM3D/3D_GLOBE) or any
        acceptable values of `format` in ObsPy's `read_events` function.
    :rtype event: obspy.core.event.Event or None
    :return event: event object if found, else None.
    """
    filename = os.path.basename(fid)
    fmt = fmt.upper()

    # Allow input of various types of source files not allowed in ObsPy
    if fmt == "SOURCE":
        cat = [read_specfem2d_source(fid)]
    elif fmt == "FORCESOLUTION":
        cat = [read_forcesolution(fid)]
    # ObsPy can handle QuakeML and CMTSOLUTION
    else:
        cat = read_events(fid, format=fmt)
        if len(cat) != 1:
            logger.warning(f"Catalog contains {len(cat)} events when only 1 "
                           f"was expected. Returning zeroth index")

    return cat[0]


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
