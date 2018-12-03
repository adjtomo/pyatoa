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


def trimstreams(st_a, st_b, force=None):
    """
    Trim two streams to common start and end times, do some basic preprocessing
    before trimming. Allows user to force one stream to conform to another

    :type st_?: obspy.stream.Stream
    :param st_?: streams to be trimmed
    :rtype st_?: obspy.stream.Stream
    :return st_?: trimmed streams
    """
    if force:
        force = force.lower()
        if force == "a":
            start_set = st_a[0].stats.starttime
            end_set = st_a[0].stats.endtime
        elif force == "b":
            start_set = st_b[0].stats.starttime
            end_set = st_b[0].stats.endtime
    else:
        st_trimmed = st_a.copy() + st_b.copy()
        start_set, end_set = 0, 1E10
        for st in st_trimmed:
            start_hold = st.stats.starttime
            end_hold = st.stats.endtime
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


def preproc(st, inv=None, **kwargs):
    """
    Preprocess waveform data. Assumes synthetics are in units of displacement.

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
    :type corners: int
    :param corners: value of the filter corners, i.e. steepness of filter edge
    :rtype st: obspy.stream.Stream
    :return st: preprocessed stream object
    """
    warnings.filterwarnings("ignore", category=FutureWarning)
    resample = kwargs.get("resample", None)
    pad_length_in_seconds = kwargs.get("pad_length_in_seconds", 20)
    output = kwargs.get("output", "VEL").upper()
    back_azimuth = kwargs.get("back_azimuth", None)
    filter_bounds = kwargs.get("filter_bounds", (10,30))
    water_level = kwargs.get("water_level", 60)
    corners = kwargs.get("corners", 4)

    if resample:
        st.resample(resample)
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05)
    if inv:
        st.attach_response(inv)
        st.remove_response(output=output, water_level=water_level, plot=False)
        logger.info("removing response with water level {}".format(water_level))
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=0.05)
        st.rotate(method="->ZNE", inventory=inv)
    # no inventory means synthetic data
    elif not inv:
        # TO DO: dynamically check what the raw specfem output units are
        # for now we assume they are velocity
        st.integrate(method="cumtrapz")
        # if output != self.raw_syn_comp:
        #     diff_dict = {"DISP":1, "VEL":2, "ACC":3}
        #     desired = diff_dict[output]
        #     given = diff_dict[self.raw_syn_comp]
        #     separation = desired - given
        #     if separation == 1:
        #         st.integrate(method="cumtrapz")
        #     elif separation == 2:
        #         st.integrate(method="cumtrapz").integrate(method="cumtrapz")
        #     elif separation == -1:
        #         st.differentiate(method="gradient")
        #     elif separation == -2:
        #         st.differentiate(method="gradient").differentiate(
        #             method="gradient")

        st.taper(max_percentage=0.05)
    if back_azimuth is not None:
        st.rotate(method="NE->RT", back_azimuth=back_azimuth)
        logger.info("rotating NE->RT by {} degrees".format(back_azimuth))
    # st = _zero_pad_stream(st, pad_length_in_seconds)
    # logger.info("zero padding front and back by {}s".format(
    #     pad_length_in_seconds))
    if filter_bounds is not None:
        st.filter('bandpass', freqmin=1/filter_bounds[1],
                  freqmax=1/filter_bounds[0], corners=corners, zerophase=True)
        msg = "filtering streams at {t0}s to {t1}s with a {c} corner {f} filter"
        logger.info(msg.format(t0=filter_bounds[0], t1=filter_bounds[1],
                               c=corners, f="Butterworth"))

    return st

