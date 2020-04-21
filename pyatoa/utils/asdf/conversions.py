"""
Functions to convert ASDF dataset data into directory structures and human
readable formats. This makes it easy to get back into the 'standard'
inversion format where data is stored in directory structures and auxiliary
data is stored in human readable JSON files or XML, csv etc.

kwargs available for dump:
    station_format
    station_directory
    waveform_directory
"""
import os
from pyatoa.utils.form import event_name


def dump(ds, **kwargs):
    """
    Convenience function to dump everything inside a dataset
    :return:
    """
    raise NotImplementedError


def dir_check(directory):
    """
    Convenience function to check and make a directory if it doesn't exist

    :type directory: str
    :param directory: directory to check
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def stations(ds, path="./", **kwargs):
    """
    Convert stations into XML files
    :param ds:
    :return:
    """
    path_tag = kwargs.get("station_directory", "stations")
    fmt = kwargs.get("station_format", "STATIONXML")

    # Set up the directory structure
    sta_dir = os.path.join(path, path_tag)
    dir_check(sta_dir)

    for sta in ds.waveforms.list():
        stationxml = ds.waveforms[sta].StationXML
        stationxml.write(os.path.join(sta_dir, sta), fmt)


def waveforms(ds, path="./", **kwargs):
    """
    Convert waveforms into individual files
    :param ds:
    :return:
    """
    path_tag = kwargs.get("waveform_directory", "waveforms")
    fmt = kwargs.get("waveform_format", "MSEED")

    # Set up the directory structure
    for sta in ds.waveforms.list():
        tags = ds.waveforms[sta].get_waveform_tags()
        for tag in tags:
            st = ds.waveforms[sta][tag]
            # Set up the directory structure for this tag
            wav_dir = os.path.join(path, path_tag, event_name(ds), tag)
            dir_check(wav_dir)
            # Make sure separators are underscores
            st.write(os.path.join(wav_dir, sta.replace(".", "_")), format=fmt)


def windows(ds, path="./", **kwargs):
    """

    :param ds:
    :param path:
    :param kwargs:
    :return:
    """
    pass




