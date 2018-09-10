"""
Pre processing functionality to put raw seismic waveoforms into the proper
format for use in analysis
"""
import numpy as np
from obspy import Stream

def _zero_pad_stream(st,pad_length_in_seconds):
    """zero pad the data of a stream, change the starttime to reflect the change
    """
    for tr in st:
        array = tr.data
        pad_width = pad_length_in_seconds * tr.stats.sampling_rate
        tr.data = np.pad(array, pad_width, mode='constant')
        tr.stats.starttime -= pad_length_in_seconds
    return st


def preproc(st, inv=None, resample=50, pad_length_in_seconds=20,
               output="VEL", filter=False):
    """
    preprocess waveform data
    """
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
    if filter:
        st.filter('bandpass', freqmin=filter[0], freqmax=filter[1], corners=4,
                  zerophase=True)
    return st

