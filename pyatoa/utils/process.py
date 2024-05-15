"""
Tools for processing obspy.Stream or obspy.Trace objects
Used for preprocessing data through filtering and tapering, zero padding etc.
Also contains tools for synthetic traces such as source time function
convolutions
"""
import numpy as np
from pyatoa import logger


def apply_filter(st, min_period=None, max_period=None, min_freq=None,
                 max_freq=None, corners=2, zerophase=True, **kwargs):
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
    if min_period is None and max_period is None:
        logger.info(f"no filter bounds given, no filtering will be applied")
        return st

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
        logger.debug(f"applying 'bandpass' filter: "
                     f"{min_period} - {max_period}s w/ {corners} corners")

    # Minimum period only == lowpass filter
    elif min_period:
        st.filter("lowpass", freq=max_freq, corners=corners,
                  zerophase=zerophase, **kwargs)
        logger.debug(f"applying 'lowpass' filter: "
                     f"{min_period}s w/ {corners} corners")

    # Maximum period only == highpass filter
    elif max_period:
        st.filter("highpass", freq=min_freq, corners=corners, 
                  zerophase=zerophase, **kwargs)
        logger.debug(f"applying 'highpass' filter: "
                     f"{max_period}s w/ {corners} corners")

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


def trim_streams(st_a, st_b, force=None):
    """
    Trim two streams to common start and end times, allowing user to `force` 
    one stream to conform to another. Assumes all traces in a stream have the 
    same time.
    
    :type st_a: obspy.stream.Stream
    :param st_a: streams to be trimmed
    :type st_b: obspy.stream.Stream
    :param st_b: streams to be trimmed
    :type force: str
    :param force: "a" or "b"; force trim to the length of "st_a" or to "st_b",
        if not given, trims to the common time
    :rtype: tuple of obspy.stream.Stream
    :return: trimmed stream objects in the same order as input
    :raises AssertionError: if the streams cannot be trimmed successfully
    """
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
                logger.info(f"shifting {tr.id} starttime by {dt}s")
                tr.stats.starttime = start_set
            elif dt >= tr.stats.delta:
                logger.warning(f"{tr.id} starttime is {dt}s greater than delta")

    for tr_a, tr_b in zip(st_a_out, st_b_out):
        assert(tr_a.stats.starttime == tr_b.stats.starttime), \
            "unable to trim streams and match starttimes"
        assert(tr_a.stats.endtime == tr_b.stats.endtime), \
            "unable to trim streams and match endtimes"

    return st_a_out, st_b_out


def match_npts(st_a, st_b, force=None):
    """
    Resampling can cause sample number differences which will lead to failure
    of some preprocessing steps. This function ensures that `npts` 
    matches between traces by extending one of the traces with zeros, or 
    removing data from the end of the trace. 

    This is used for the case where the 'obs' data is too short w.r.t 'syn' 
    data, so we pad out the end of the 'obs' with 0s. 

    .. warning::

        It is assumed you have resampled and trimmed the streams and that this
        function is only used to make up sub-second differences in the number of
        samples. Also assumed you will taper the end of the trace otherwise the
        appending of zeros will cause issues.

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
    :raises AssertionError: if the number of points cannot be matched
    """
    # Quick check to make sure all traces have the same length within one stream
    for st in [st_a, st_b]:
        for tr in st:
            assert(tr.stats.npts == st[0].stats.npts), \
                "all traces in stream must have the same number of samples"

    # Assign the number of points, copy to avoid editing in place
    if force is None or force == "a":
        npts = st_a[0].stats.npts
        st_const = st_a.copy()
        st_change = st_b.copy()
    else:
        npts = st_b[0].stats.npts
        st_const = st_b.copy()
        st_change = st_a.copy()

    for tr in st_change:
        diff = npts - tr.stats.npts 
        if diff > 0:
            logger.info(f"appending {diff} zeros ({diff * tr.stats.delta}s) to "
                        f"{tr.get_id()} to match npts")
            tr.data = np.append(tr.data, np.zeros(diff))
        elif diff < 0:
            logger.info(f"removing {diff} zeros ({diff * tr.stats.delta}s) "
                        f"from {tr.get_id()} to match npts")
            tr.data = tr.data[:diff]
        elif diff == 0:
            continue

    for tr_a, tr_b in zip(st_const, st_change):
        assert(tr_a.stats.npts == tr_b.stats.npts), \
            "unable to match npts between streams"

    # Ensure streams are returned in the correct order
    if not force or force == "a":
        return st_const, st_change
    else:
        return st_change, st_const


def normalize(st_a, st_b, choice):
    """
    Normalize amplitudes in traces to the absolute maximum amplitude of the
    other, or normalize all traces so that their maximum value is a given value.

    :type st_a: obspy.stream.Stream
    :param st_a: one stream to match samples with
    :type st_b: obspy.stream.Stream
    :param st_b: one stream to match samples with
    :type choice: str or int or float
    :param choice: choose which stream to use as the default npts,
        defaults to 'a', options: 'a', 'b' or a value that you want to set the
        max amplitude to be equal to, e.g., 1
    :rtype: tuple (obspy.stream.Stream, obspy.stream.Stream)
    :return: streams that have been normalized based on `choice`
    :raises NotImplementedError: If incorrect value of `choice` is provided
    """
    if isinstance(choice, str):
        if choice == "a":
            st_const = st_a.copy()
            st_change = st_b.copy()
        elif choice == "b":
            st_change = st_a.copy()
            st_const = st_b.copy()
        else:
            raise NotImplementedError("normalize `choice` must be 'a' or 'b'")

        for tr_const, tr_change in zip(st_const, st_change):
            tr_change.data *= (abs(tr_const.max()) / abs(tr_change.max()) )

        if choice == "a":
            st_a_out = st_const
            st_b_out = st_change
        else:
            st_a_out = st_change
            st_b_out = st_const
    elif isinstance(choice, (int, float)):
        st_a_out = st_a.copy()
        st_b_out = st_b.copy()

        for tr_a, tr_b in zip(st_a, st_b):
            tr_a.data *= (choice / abs(tr_a.max()))
            tr_b.data *= (choice / abs(tr_b.max()))
    else:
        raise NotImplementedError(f"invalid choice {choice} for normalization")

    return st_a_out, st_b_out


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
    logger.info(f"convolving data w/ Gaussian (t/2={half_duration:.2f}s)")

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
    for tr in st:
        if time_shift:
            tr.stats.starttime += time_shift
        data_out = np.convolve(tr.data, gaussian_stf, mode="same")
        tr.data = data_out

    return st



