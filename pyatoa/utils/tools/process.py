"""
Tools for processing obspy.Stream or obspy.Trace objects
Used for preprocessing data through filtering and tapering, zero padding etc.
Also contains tools for synthetic traces such as source time function
convolutions
"""
import warnings
import numpy as np

from pyatoa import logger


def zero_pad_stream(st, pad_length_in_seconds):
    """
    Zero pad the data of a stream, change the starttime to reflect the change

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
    Trim two streams to common start and end times,
    Do some basic preprocessing before trimming.
    Allows user to force one stream to conform to another.
    Assumes all traces in a stream have the same time.
    Prechecks make sure that the streams are actually different

    :type st_a: obspy.stream.Stream
    :param st_a: streams to be trimmed
    :type st_b: obspy.stream.Stream
    :param st_b: streams to be trimmed
    :type force: str
    :param force: "a" or "b"; force trim to the length of "st_a" or to "st_b",
        if not given, trims to the common time
    :rtype st_?: obspy.stream.Stream
    :return st_?: trimmed stream
    """
    # Check if the times are already the same
    diff = 1E-2
    if st_a[0].stats.starttime - st_b[0].stats.starttime < diff and \
            st_a[0].stats.endtime - st_b[0].stats.endtime < diff:
        return st_a, st_b

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


def is_preprocessed(st):
    """
    Small check to make sure a stream object has not yet been run through
    preprocessing. Simple, as it assumes that a fresh stream will have no
    processing attribute in their stats, or if they do, will not have been
    filtered (getting cut waveforms from FDSN appends a 'trim' stat).
    :type st: obspy.stream.Stream
    :param st: stream to check processing on
    :rtype: bool
    :return: if preprocessing has occurred
    """
    for tr in st:
        if hasattr(tr.stats, 'processing'):
            for processing in tr.stats.processing:
                # A little hacky, but processing flag will have the str
                # ..': filter(options'... to signify that a filter is applied
                if 'filter(' in processing:
                    warnings.warn("stream already preprocessed", UserWarning)
                    return True

    # If nothing found, return False
    return False


def preproc(st_original, inv=None, unit_output="VEL", synthetic_unit=None,
            back_azimuth=None, filter_bounds=(10, 30), water_level=60,
            corners=4, taper_percentage=0.05):
    """
    Preprocess waveform data. Because synthetics can be in any unit,
    requires user to specify what units they are.

    :type st_original: obspy.stream.Stream
    :param st_original: stream object to process
    :type inv: obspy.core.inventory.Inventory or None
    :param inv: inventory containing relevant network and stations
    :type unit_output: str
    :param unit_output: output of response removal, available:
        'DISP', 'VEL', 'ACC'
    :type synthetic_unit: str
    :param synthetic_unit: units of synthetic traces, same available as unit
    :type back_azimuth: float
    :param back_azimuth: back azimuth in degrees
    :type filter_bounds: list or tuple of float
    :param filter_bounds: (min period, max_period), units of Seconds
    :type water_level: int
    :param water_level: water level for response removal
    :type corners: int
    :param corners: value of the filter corners, i.e. steepness of filter edge
    :type taper_percentage: float
    :param taper_percentage: amount to taper ends of waveform
    :rtype st: obspy.stream.Stream
    :return st: preprocessed stream object
    """
    # Copy the stream to avoid editing in place
    st = st_original.copy()
    if is_preprocessed(st):
        return st

    # Standard preprocessing
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=taper_percentage)

    # If inventory is given, assume working with observation data
    if inv:
        # Occasionally, inventory issues arise, as ValueErrors due to 
        # station availability, e.g. NZ.COVZ. Try/except to catch these.
        try:
            st.attach_response(inv)
            st.remove_response(output=unit_output,
                               water_level=water_level,
                               plot=False)
        except ValueError:
            logger.debug(f"Error removing response from {st[0].get_id()}")
            return None
    
        logger.debug("remove response, units of {}".format(unit_output))

        # Clean up streams after response removal
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=taper_percentage)

        # Rotate streams if they are not in the ZNE coordinate system
        st.rotate(method="->ZNE", inventory=inv)

    # No inventory means synthetic data
    elif not inv:
        if unit_output != synthetic_unit:
            logger.debug("unit output and synthetic output do not match, "
                         "adjusting")
            st = change_syn_units(st, current=unit_output,
                                  desired=synthetic_unit)
            st.detrend("linear")
            st.detrend("demean")

        st.taper(max_percentage=taper_percentage)
    
    # Rotate the given stream from standard NEZ to RTZ
    if back_azimuth:
        st.rotate(method="NE->RT", back_azimuth=back_azimuth)
        logger.debug(f"rotating NE->RT by {back_azimuth} degrees")

    # Filter data using ObsPy Butterworth filters. Zerophase avoids phase shift
    if filter_bounds is not None:
        st.filter('bandpass', freqmin=1/filter_bounds[1],
                  freqmax=1/filter_bounds[0], corners=corners, zerophase=True
                  )
        logger.debug(f"filter {filter_bounds[0]}s to {filter_bounds[1]}s")

    return st


def change_syn_units(st, current, desired):
    """
    Change the synthetic units based on the desired unit output and the current
    unit output

    :type st: obspy.stream.Stream
    :param st: obspy stream with synthetic data
    :type current: str
    :param current: current units, 'DISP', 'VEL', or 'ACC'
    :type desired: str
    :param desired: desired unit output key, same keys as `current`
    :rtype: obspy.stream.Stream
    :return: stream with desired units
    """
    # Copy to avoid editing in place
    st_diff = st.copy()

    # Determine the difference between synthetic unit and observed unit
    diff_dict = {"DISP": 1, "VEL": 2, "ACC": 3}
    difference = diff_dict[current] - diff_dict[desired]

    # Integrate or differentiate stream to retrieve correct units
    if difference == 1:
        st_diff.integrate(method="cumtrapz")
    elif difference == 2:
        st_diff.integrate(method="cumtrapz").integrate(method="cumtrapz")
    elif difference == -1:
        st_diff.differentiate(method="gradient")
    elif difference == -2:
        st_diff.differentiate(method="gradient").differentiate(
            method="gradient")

    return st_diff


def stf_convolve(st, half_duration, source_decay=4., time_shift=None,
                 time_offset=None):
    """
    Convolve function with a Gaussian window.
    Following taken from specfem "comp_source_time_function.f90"

    hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
    with SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68

    This gaussian uses a strong decay rate to avoid non-zero onset times, while
    still miicking a triangle source time function

    :type st: obspy.stream.Stream
    :param st: stream object to convolve with source time function
    :type half_duration: float
    :param half_duration: the half duration of the source time function,
        usually provided in moment tensor catalogs
    :type source_decay: float
    :param source_decay: the decay strength of the source time function, the
        default value of 4 gives a Gaussian. A value of 1.68 mimics a triangle.
    :type time_shift: float
    :param time_shift: Time shift of the source time function in seconds
    :type time_offset: If simulations have a value t0 that is negative, i.e. a
        starttime before the event origin time. This value will make sure the
        source time function doesn't start convolving before origin time to
        avoid non-zero onset times
    :rtype: obspy.stream.Stream
    :return: stream object which has been convolved with a source time function
    """
    logger.debug(f"convolve w/ gaussian half-dur={half_duration:.2f}s")

    sampling_rate = st[0].stats.sampling_rate
    half_duration_in_samples = round(half_duration * sampling_rate)

    # generate gaussian function
    decay_rate = half_duration_in_samples / source_decay
    a = 1 / (decay_rate ** 2)
    t = np.arange(-half_duration_in_samples, half_duration_in_samples, 1)
    gaussian_stf = np.exp(-a * t**2) / (np.sqrt(np.pi) * decay_rate)

    # prepare time offset machinery
    if time_offset:
        time_offset_in_samp = int(time_offset * sampling_rate)

    # convolve each trace with the soure time function and time shift if needed
    st_out = st.copy()
    for tr in st_out:
        if time_shift:
            tr.stats.starttime += time_shift
        # if time offset is given, split the trace into before and after origin
        # only convolve after zero time

        # This doesn't actually work? It just adds more samples to the trace
        # which is not what we want
        # if time_offset:
        #     after_origin = np.convolve(tr.data[time_offset_in_samp:],
        #                                gaussian_stf, mode="same"
        #                                )
        #     data_out = np.concatenate([tr.data[:time_offset_in_samp],
        #                                after_origin])
        #     tr.data = data_out
        # else:
        #     data_out = np.convolve(tr.data, gaussian_stf, mode="same")
        #     tr.data = data_out
        data_out = np.convolve(tr.data, gaussian_stf, mode="same")
        tr.data = data_out

    return st_out
