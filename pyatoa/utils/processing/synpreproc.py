"""
a suite of functions used to preprocess synthetic seismograms so that they
are comparable to realworld observations
"""
import numpy as np
from scipy import signal

def stf_convolve(st,half_duration,window="bartlett",time_shift=False):
    """convolve source time function with a stream, time shift if needed
    :type st: obspy.stream
    :param st: stream object containing traces of data
    :type half_duration: float
    :param half_duration: half duration of stf in seconds
    :type window: str
    :param window: window type to return
    :type time_shift: float
    :param time_shift: change the starttime
    ========================================
    boxcar, triang, blackman, hamming, hann,
    bartlett, flattop, parzen, bohman,
    blackmanharris, nuttall, barthann,
    kaiser (needs beta),
    gaussian (needs standard deviation),
    general_gaussian (needs power, width),
    slepian (needs width),
    chebwin (needs attenuation),
    exponential (needs decay scale),
    tukey (needs taper fraction)
    NOTE: bartlett window is a triangle that touches 0
    ========================================
    :return new_st:
    """
    npts = st[0].stats.npts
    sampling_rate = st[0].stats.sampling_rate
    half_duration_in_samples = round(half_duration * sampling_rate)
    stf = signal.get_window(window=window,
                            Nx=(half_duration_in_samples * 2) -1)

    # make sure window touches 0 at the end
    if stf[-1] != 0:
        stf = np.append(stf,0)
    stf *= (2/len(stf))
    st_out = st.copy()
    for tr in st_out:
        if time_shift:
            tr.stats.starttime = tr.stats.starttime + time_shift
        data_out = np.convolve(tr.data,stf,mode="same")
        tr.data = data_out
    return new_st


def timeshift_halfduration(event_id):
    """
    get the absolute time shift between centroid time and hypocenter time to
    shift the synthetic seismogram into absolute time using the equation:

    t_abs = t_pde + time shift + t_syn

    also get the half duration from the GCMT solution
    """
    MT = getdata.get_GCMT_solution(event_id)
    CMTSOLUTIONPATH = (pathnames()['kupedata'] +
                            'CMTSOLUTIONS/{}CMTSOLUTION'.format(event_id))
    CMTSOLUTION = read_events(CMTSOLUTIONPATH)

    CMTSOLUTION_time = CMTSOLUTION[0].origins[0].time
    CENTROID_time = [i.time for i in MT.origins
                                        if i.origin_type == "centroid"][0]

    time_shift = abs(CMTSOLUTION_time - CENTROID_time)

    # half duration
    moment_tensor = MT.focal_mechanisms[0].moment_tensor
    half_duration = (moment_tensor.source_time_function['duration'])/2

    return time_shift, half_duration

