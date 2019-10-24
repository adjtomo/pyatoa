"""
a suite of functions used to preprocess synthetic seismograms so that they
are comparable to realworld observations
"""
import numpy as np
from scipy import signal

from pyatoa import logger


def stf_convolve(st, half_duration, window="bartlett", time_shift=None):
    """
    Generic convolve source time function with a stream, time shift if needed

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
    logger.debug("convolve syn  w/ {0} of T_half={1:.2f}s".format(
                                                          window, half_duration)
                 )

    # make sure window touches 0 at the end
    if stf[-1] != 0:
        stf = np.append(stf, 0)
    # stf *= (2/len(stf))  # not sure if this normalization is necessary?
    st_out = st.copy()
    for tr in st_out:
        if time_shift:
            tr.stats.starttime += time_shift
        data_out = np.convolve(tr.data, stf, mode="same")
        tr.data = data_out

        # Add a processing Tag
        tr.stats.processing += ["Pyatoa: convolve(window={0}::time_shift={1}::"
                                "half_duration_in_samples={2}".format(
            window, time_shift, half_duration_in_samples)
        ]
    return st_out


def stf_convolve_gaussian(st, half_duration, time_shift=None):
    """
    Convolve function with a Gaussian window. 
    Following taken from specfem "comp_source_time_function.f90"

    hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
    with SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68

    This gaussian uses a strong decay rate to avoid non-zero onset times, while
    still miicking a triangle source time function

    :param st:
    :param half_duration:
    :param time_shift:
    :return:
    """
    logger.debug("convolving synthetic data with gaussian "
                 "window of half duration {:.2f}s".format(half_duration)
                 )
    sampling_rate = st[0].stats.sampling_rate
    half_duration_in_samples = round(half_duration * sampling_rate)

    # generate gaussian function
    source_decay = 4
    # source_decay_mimic_triangle = 1.68
    decay_rate = half_duration_in_samples / source_decay
    a = 1 / (decay_rate ** 2)
    t = np.arange(-half_duration_in_samples, half_duration_in_samples, 1)
    gaussian_stf = np.exp(-a * t**2) / (np.sqrt(np.pi) * decay_rate)

    st_out = st.copy()
    for tr in st_out:
        if time_shift:
            tr.stats.starttime += time_shift
        data_out = np.convolve(tr.data, gaussian_stf, mode="same")
        tr.data = data_out

    return st_out




