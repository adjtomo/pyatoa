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
            output="VEL", back_azimuth=None, filterbounds=None, water_level=60):
    """
    preprocess waveform data
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

