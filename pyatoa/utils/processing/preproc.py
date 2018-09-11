"""
Pre processing functionality to put raw seismic waveoforms into the proper
format for use in analysis
"""
import warnings
import numpy as np

from pyatoa import logger


def _zero_pad_stream(st, pad_length_in_seconds):
    """zero pad the data of a stream, change the starttime to reflect the change
    """
    for tr in st:
        array = tr.data
        pad_width = int(pad_length_in_seconds * tr.stats.sampling_rate)
        tr.data = np.pad(array, pad_width, mode='constant')
        tr.stats.starttime -= pad_length_in_seconds
    return st


def trimstreams(st_a, st_b):
    """trim streams to common start and end times
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
            output="VEL", filter=False):
    """
    preprocess waveform data
    """
    logger.info("ignoring FutureWarnings due to obspy trace warnings")
    warnings.filterwarnings("ignore", category=FutureWarning)

    st.resample(resample)
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05)
    if inv:
        st.attach_response(inv)
        st.remove_response(output=output, water_level=60, plot=False)
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
    st = _zero_pad_stream(st, pad_length_in_seconds)
    logger.info("zero padding by {}s".format(pad_length_in_seconds))
    if filter:
        st.filter('bandpass', freqmin=1/filter[1], freqmax=1/filter[0],
                  corners=4, zerophase=True)
        logger.info("filtering stream at {}s to {}s".format(
            filter[0],filter[1])
        )
    return st

