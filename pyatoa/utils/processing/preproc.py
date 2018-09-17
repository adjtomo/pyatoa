"""
Pre processing functionality to put raw seismic waveoforms into the proper
format for use in analysis
"""
import warnings
import numpy as np

from pyatoa import logger


def _zero_pad_stream(st, pad_length_in_seconds):
    """
    zero pad the data of a stream, change the starttime to reflect the change

    :type st: obspy.stream.Stream
    :param st: stream to be zero padded
    :type pad_length_in_seconds: int
    :param pad_length_in_seconds: length of padding front and back
    :rtype st: obspy.stream.Stream
    :return st: stream with zero padded data object
    """
    for tr in st:
        array = tr.data
        pad_width = int(pad_length_in_seconds * tr.stats.sampling_rate)
        tr.data = np.pad(array, pad_width, mode='constant')
        tr.stats.starttime -= pad_length_in_seconds
    return st


def trimstreams(st_a, st_b):
    """
    Trim two streams to common start and end times, do some basic preprocessing
    before trimming.

    :type st_?: obspy.stream.Stream
    :param st_?: streams to be trimmed
    :rtype st_?: obspy.stream.Stream
    :return st_?: trimmed streams
    """
    st_trimmed = st_a.copy() + st_b.copy()
    start_set, end_set = 0, 1E10
    for tr in st_trimmed:
        start_hold = tr.stats.starttime
        end_hold = tr.stats.endtime
        if start_hold > start_set:
            start_set = start_hold
        if end_hold < end_set:
            end_set = end_hold
    for st in [st_a, st_b]:
        st.trim(start_set, end_set)
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=0.05)

    return st_a, st_b


def preproc(st, inv=None, resample=5, pad_length_in_seconds=20,
            output="VEL", back_azimuth=None, filterbounds=None, water_level=60):
    """
    Preprocess waveform data

    :type st: obspy.stream.Stream
    :param st: stream object to process
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing relevant network and stations
    :type resample: int
    :param resample: sampling rate to resample to in Hz
    :type pad_length_in_seconds: int
    :param pad_length_in_seconds: length of padding front and back
    :type output: str
    :param output: output of response removal, available: 'DISP', 'VEL', 'ACC'
    :type back_azimuth: float
    :param back_azimuth: back azimuth in degrees
    :type filterbounds: list of float
    :param filterbounds: (min period, max_period)
    :type water_level: int
    :param water_level: water level for response removal
    :rtype st: obspy.stream.Stream
    :return st: preprocessed stream object
    """
    warnings.filterwarnings("ignore", category=FutureWarning)

    st.resample(resample)
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05)
    if inv:
        st.attach_response(inv)
        st.remove_response(output=output, water_level=water_level, plot=False)
        import matplotlib.pyplot as plt; plt.show()
        logger.info("removing response with water level {}".format(water_level))
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=0.05)
    # no inventory means synthetic data
    elif not inv:
        if output == "DISP":
            st.integrate(method="cumtrapz")
        elif output == "ACC":
            st.differentiate(method="gradient")
        st.taper(max_percentage=0.05)
    if back_azimuth is not None:
        st.rotate(method="NE->RT", back_azimuth=back_azimuth)
        logger.info("rotating NE->RT by {} degrees".format(back_azimuth))
    st = _zero_pad_stream(st, pad_length_in_seconds)
    logger.info("zero padding front and back by {}s".format(
        pad_length_in_seconds))
    if filterbounds is not None:
        st.filter('bandpass', freqmin=1/filterbounds[1],
                  freqmax=1/filterbounds[0], corners=4, zerophase=True)
        logger.info("filtering stream at {}s to {}s".format(filterbounds[0],
                                                            filterbounds[1])
                    )
    return st

