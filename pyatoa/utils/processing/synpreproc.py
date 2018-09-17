"""
a suite of functions used to preprocess synthetic seismograms so that they
are comparable to realworld observations
"""
import numpy as np
from scipy import signal

from pyatoa import logger


def stf_convolve(st, half_duration, window="bartlett", time_shift=False):
    """
    Convolve source time function with a stream, time shift if needed

    :type st: obspy.stream.Stream
    :param st: stream object containing traces of data
    :type half_duration: float
    :param half_duration: half duration of stf in seconds
    :type window: str
    :param window: window type to return
    :type time_shift: float
    :param time_shift: change the starttime
    :rtype st_out: obspy.stream.Stream
    :return st_out: stream object with data convolved with stf
    """
    sampling_rate = st[0].stats.sampling_rate
    half_duration_in_samples = round(half_duration * sampling_rate)
    stf = signal.get_window(window=window,
                            Nx=(half_duration_in_samples * 2) - 1)
    logger.info("convolving synthetic data with "
                "{0} window of half duration {1:.2f}s".format(window,
                                                              half_duration)
                )

    # make sure window touches 0 at the end
    if stf[-1] != 0:
        stf = np.append(stf, 0)
    stf *= (2/len(stf))
    st_out = st.copy()
    for tr in st_out:
        if time_shift:
            tr.stats.starttime += time_shift
        data_out = np.convolve(tr.data, stf, mode="same")
        tr.data = data_out
    return st_out


def half_duration_from_m0(moment):
    """
    Empirical formula for half duration used by Harvard CMT, stated in
    Daheln and Tromp (1998, p.178).

    :type moment: float
    :param moment: seismic moment in N*m
    :rtype: float
    :return: empirically scaled half duration
    """
    return 2.4E-6 * moment**(1/3)

