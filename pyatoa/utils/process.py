"""
Tools for processing obspy.Stream or obspy.Trace objects
Used for preprocessing data through filtering and tapering, zero padding etc.
Also contains tools for synthetic traces such as source time function
convolutions
"""
import numpy as np
from pyatoa import logger


def default_process(mgmt, choice, **kwargs):
    """
    Default preprocessing function to process  waveform data from a Manager
    Preprocessing is slightly different for obs and syn waveforms. Each
    processing function is split into a separate function so that they can
    be called by custom preprocessing functions.

    :type mgmt: pyatoa.core.manager.Manager
    :param mgmt: Manager class that should contain a Config object as well as
        waveform data and inventory
    :type choice: str
    :param choice: option to preprocess observed, synthetic or both
        available: 'obs', 'syn'
    :rtype: obspy.core.stream.Stream
    :return: preprocessed stream object pertaining to `choice`

    Keyword Arguments
    ::
        int water_level:
            water level for response removal
        float taper_percentage:
            amount to taper ends of waveform
        bool remove_response:
            remove instrument response using the Manager's inventory object.
            Defaults to True
        bool apply_filter:
            filter the waveforms using the Config's min_period and max_period
            parameters. Defaults to True
        bool convolve_with_stf:
            Convolve synthetic data with a Gaussian source time function if a
            half duration is provided.
    """
    assert choice in ["obs", "syn"], "choice must be 'obs' or 'syn"

    water_level = kwargs.get("water_level", 60)
    taper_percentage = kwargs.get("taper_percentage", 0.05)
    zerophase = kwargs.get("zerophase", True)
    remove_response = kwargs.get("remove_response", True)
    apply_filter = kwargs.get("apply_filter", True)
    convolve_with_stf = kwargs.get("convolve_with_stf", True)

    # Copy the stream to avoid editing in place. Synthetic variable used to
    # denote if the waveforms are synthetic or not, these require special
    # processing steps.
    if choice == "syn":
        st = mgmt.st_syn.copy()
        is_synthetic_data = True
    elif choice == "obs":
        st = mgmt.st_obs.copy()
        is_synthetic_data = mgmt.config.synthetics_only

    if is_preprocessed(st):
        return st

    # Get rid of any long period trends that may affect that data
    st.detrend("simple").detrend("demean").taper(taper_percentage)
    st = taper_time_offset(st, taper_percentage, mgmt.stats.time_offset_sec)

    # Observed specific data preprocessing includes response and rotating to ZNE
    if remove_response and not is_synthetic_data:
        logger.debug(f"removing response, units to {mgmt.config.unit_output}")
        st.remove_response(inventory=mgmt.inv, output=mgmt.config.unit_output,
                           water_level=water_level, plot=False)

        # Rotate streams if not in ZNE, e.g. Z12. Only necessary for observed
        logger.debug("rotating from generic coordinate system to ZNE")
        st.rotate(method="->ZNE", inventory=mgmt.inv)
        st.detrend("simple").detrend("demean").taper(taper_percentage)
    else:
        logger.debug("no response removal, synthetic data or requested not to")

    # Rotate the given stream from standard NEZ to RTZ if BAz given
    if mgmt.baz:
        logger.debug(f"rotating NE->RT by {mgmt.baz} degrees")
        st.rotate(method="NE->RT", back_azimuth=mgmt.baz)

    # Filter data based on the given period bounds
    if apply_filter:
        st = filters(st, min_period=mgmt.config.min_period,
                     max_period=mgmt.config.max_period, 
                     corners=mgmt.config.filter_corners,
                     zerophase=zerophase
                     )
        st.detrend("simple").detrend("demean").taper(taper_percentage)
    else:
        logger.debug(f"no filter applied to data")

    # Convolve synthetic data with a Gaussian source time function
    if convolve_with_stf and is_synthetic_data and mgmt.stats.half_dur:
        st = stf_convolve(st=st, half_duration=mgmt.stats.half_dur)

    return st


def filters(st, min_period=None, max_period=None, min_freq=None, max_freq=None,
            corners=2, zerophase=True, **kwargs):
    """
    Choose the appropriate filter depending on the ranges given.
    Either periods or frequencies can be given. Periods will be prioritized.
    Uses Butterworth filters by default.

    Filters the stream in place. Kwargs passed to filter functions.

    :type st: obspy.core.stream.Stream
    :param st: stream object to be filtered
    :type min_period: float
    :param min_period: minimum filter bound in units of seconds
    :type max_period: float
    :param max_period: maximum filter bound in units of seconds
    :type min_freq: float
    :param min_freq: optional minimum filter bound in units of Hz, will be
        overwritten by `max_period` if given
    :type max_freq: float
    :param max_freq: optional maximum filter bound in units of Hz, will be
        overwritten by `min_period` if given
    :type corners: int
    :param corners: number of filter corners to be passed to ObsPy
        filter functions
    :type zerophase: bool
    :param zerophase: if True, run filter backwards and forwards to avoid
        any phase shifting
    :rtype: obspy.core.stream.Stream
    :return: Filtered stream object
    """
    # Ensure that the frequency and period bounds are the same
    if not min_period and max_freq:
        min_period = 1 / max_freq
    if not max_period and min_freq:
        max_period = 1 / min_freq
    if not max_freq:
        max_freq = 1/ min_period
    if not min_freq:
        min_freq = 1 / max_period

    # Bandpass if both bounds given
    if min_period and max_period:
        st.filter("bandpass", corners=corners, zerophase=zerophase,
                  freqmin=min_freq, freqmax=max_freq, **kwargs)
        logger.debug(f"bandpass filter: {min_period} - {max_period}s w/ "
                     f"{corners} corners")

    # Highpass if only minimum period given
    elif min_period:
        st.filter("highpass", freq=max_freq, corners=corners,
                  zerophase=zerophase, **kwargs)
        logger.debug(f"highpass filter: {min_period}s w/ {corners} corners")

    # Lowpass if only minimum period given
    elif max_period:
        st.filter("lowpass", freq=min_freq, corners=corners, zerophase=True,
                  **kwargs)
        logger.debug(f"lowpass filter: {max_period}s w/ {corners} corners")

    return st


def taper_time_offset(st, taper_percentage=0.05, time_offset_sec=0):
    """
    Taper the leading edge of the waveform. If a time offset is given,
    e.g. 20s before the event origin time (T_0), taper all the way up from
    T=0 to T=T_0, to ensure that there are no impulse-like signals prior to the
    event origin.

    :type st: obspy.core.stream.Stream
    :param st: Stream object to be tapered
    :type taper_percentage: float
    :param taper_percentage: default taper percentage
    :type time_offset_sec: float
    :param time_offset_sec: Any time offset between the start of the stream to
        the event origin time. All time between these two points will be tapered
        to reduce any signals prior to the event origin.
    :rtype: obspy.core.stream.Stream
    :return: tapered Stream object
    """
    taper_amount = st[0].stats.npts * taper_percentage * st[0].stats.delta

    if taper_amount > abs(time_offset_sec):
        logger.warning("taper amount exceeds time offset, taper may affect "
                       "data if source receiver distance is short")
    elif taper_amount < abs(time_offset_sec):
        logger.info(f"adjusting taper to cover time offset {time_offset_sec}")
        taper_percentage = (abs(time_offset_sec) /
                            st[0].stats.npts * st[0].stats.delta)

    # Get rid of extra long period signals which may adversely affect processing
    st.detrend("simple").taper(taper_percentage, side="left")

    return st


def zero_pad(st, pad_length_in_seconds, before=True, after=True):
    """
    Zero pad the data of a stream, change the starttime to reflect the change.
    Useful for if e.g. observed data starttime comes in later than synthetic.

    :type st: obspy.stream.Stream
    :param st: stream to be zero padded
    :type pad_length_in_seconds: int
    :param pad_length_in_seconds: length of padding front and back
    :type before: bool
    :param before: pad the stream before the origin time
    :type after: bool
    :param after: pad the stream after the last sample
    :rtype st: obspy.stream.Stream
    :return st: stream with zero padded data object
    """
    pad_before, pad_after = 0, 0
    st_pad = st.copy()
    for tr in st_pad:
        array = tr.data
        pad_width = int(pad_length_in_seconds * tr.stats.sampling_rate)
        # Determine if we should pad before or after
        if before:
            pad_before = pad_width
        if after:
            pad_after = pad_width
        logger.debug(f"zero pad {tr.id} ({pad_before}, {pad_after}) samples")
        # Constant value is default 0
        tr.data = np.pad(array, (pad_before, pad_after), mode='constant')
        tr.stats.starttime -= pad_length_in_seconds
        logger.debug(f"new starttime {tr.id}: {tr.stats.starttime}")

    return st_pad


def trim_streams(st_a, st_b, precision=1E-3, force=None):
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
    :type precision: float
    :param precision: precision to check UTCDateTime differences
    :type force: str
    :param force: "a" or "b"; force trim to the length of "st_a" or to "st_b",
        if not given, trims to the common time
    :rtype: tuple of obspy.stream.Stream
    :return: trimmed stream objects in the same order as input
    """
    # Check if the times are already the same
    if st_a[0].stats.starttime - st_b[0].stats.starttime < precision and \
            st_a[0].stats.endtime - st_b[0].stats.endtime < precision:
        logger.debug(f"start and endtimes already match to {precision}")
        return st_a, st_b

    # Force the trim to the start and end times of one of the streams
    if force:
        if force.lower() == "a":
            start_set = st_a[0].stats.starttime
            end_set = st_a[0].stats.endtime
        elif force.lower() == "b":
            start_set = st_b[0].stats.starttime
            end_set = st_b[0].stats.endtime
    # Get starttime and endtime base on min values
    else:
        st_trimmed = st_a + st_b
        start_set, end_set = 0, 1E10
        for st in st_trimmed:
            start_hold = st.stats.starttime
            end_hold = st.stats.endtime
            if start_hold > start_set:
                start_set = start_hold
            if end_hold < end_set:
                end_set = end_hold
    
    # Trim to common start and end times    
    st_a_out = st_a.copy()
    st_b_out = st_b.copy()
    for st in [st_a_out, st_b_out]:
        st.trim(start_set, end_set)

    # Trimming doesn't always make the starttimes exactly equal if the precision
    # of the UTCDateTime object is set too high.
    # Artificially shift the starttime of the streams iff the amount shifted
    # is less than the sampling rate
    for st in [st_a_out, st_b_out]:
        for tr in st:
            dt = start_set - tr.stats.starttime
            if 0 < dt < tr.stats.sampling_rate:
                logger.debug(f"shifting {tr.id} starttime by {dt}s")
                tr.stats.starttime = start_set
            elif dt >= tr.stats.delta:
                logger.warning(f"{tr.id} starttime is {dt}s greater than delta")

    return st_a_out, st_b_out


def match_npts(st_a, st_b, force=None):
    """
    Resampling can cause sample number differences which will lead to failure
    of some preprocessing or processing steps. This function ensures that `npts` 
    matches between traces by extending one of the traces with zeros. 
    A small taper is applied to ensure the new values do not cause 
    discontinuities.

    Note:
        its assumed that all traces within a single stream have the same `npts`

    :type st_a: obspy.stream.Stream
    :param st_a: one stream to match samples with
    :type st_b: obspy.stream.Stream
    :param st_b: one stream to match samples with
    :type force: str
    :param force: choose which stream to use as the default npts,
        defaults to 'a', options: 'a', 'b'
    :rtype: tuple (obspy.stream.Stream, obspy.stream.Stream)
    :return: streams that may or may not have adjusted npts, returned in the 
        same order as provided
    """
    # Assign the number of points, copy to avoid editing in place
    if not force or force == "a":
        npts = st_a[0].stats.npts
        st_const = st_a.copy()
        st_change = st_b.copy()
    else:
        npts = st_b[0].stats.npts
        st_const = st_b.copy()
        st_change = st_a.copy()

    for tr in st_change:
        diff = abs(tr.stats.npts - npts)
        if diff:
            logger.debug(f"appending {diff} zeros to {tr.get_id()}")
            tr.data = np.append(tr.data, np.zeros(diff))

    # Ensure streams are returned in the correct order
    if not force or force == "a":
        return st_const, st_change
    else:
        return st_change, st_const


def is_preprocessed(st, filter_only=True):
    """
    Check to make sure a stream object has not yet been run through
    preprocessing.
    Assumes that a fresh stream will have no processing attribute in their
    stats, or if they do, will not have been filtered
    (getting cut waveforms from FDSN appends a 'trim' stat).

    :type st: obspy.stream.Stream
    :param st: stream to check processing on
    :type filter_only: bool
    :param filter_only: only check if the stream has been filtered, as other
        processing steps (e.g., demeaning, rotating) will also lead to a
        'processing' stat. Usually this is what you want to check as filtering
        is one of the last steps in the processing chain.
    :rtype: bool
    :return: if preprocessing has occurred
    """
    for tr in st:
        if hasattr(tr.stats, "processing"):
            if filter_only:
                for processing in tr.stats.processing:
                    # A little hacky, but processing flag will have the str
                    # ..': filter(options'... to signify that a filter is applied
                    if "filter(options" in processing:
                        return True
            else:
                return bool(tr.stats.processing)
    # If nothing found, return False
    return False


def stf_convolve(st, half_duration, source_decay=4., time_shift=None,
                 time_offset=None):
    """
    Convolve function with a Gaussian window source time function.
    Design follows Specfem3D Cartesian "comp_source_time_function.f90"

    `hdur` given is `hdur`_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
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
    logger.debug(f"convolving data w/ Gaussian (t/2={half_duration:.2f}s)")

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
        data_out = np.convolve(tr.data, gaussian_stf, mode="same")
        tr.data = data_out

    return st_out



