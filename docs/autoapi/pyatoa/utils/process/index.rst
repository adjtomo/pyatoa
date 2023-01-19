:py:mod:`pyatoa.utils.process`
==============================

.. py:module:: pyatoa.utils.process

.. autoapi-nested-parse::

   Tools for processing obspy.Stream or obspy.Trace objects
   Used for preprocessing data through filtering and tapering, zero padding etc.
   Also contains tools for synthetic traces such as source time function
   convolutions



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.process.default_process
   pyatoa.utils.process.filters
   pyatoa.utils.process.taper_time_offset
   pyatoa.utils.process.zero_pad
   pyatoa.utils.process.trim_streams
   pyatoa.utils.process.match_npts
   pyatoa.utils.process.is_preprocessed
   pyatoa.utils.process.stf_convolve



.. py:function:: default_process(st, choice, inv=None, baz=None, min_period=None, max_period=None, unit_output=None, half_dur=None, rotate_to_rtz=False, apply_filter=True, remove_response=True, convolve_with_stf=True, **kwargs)

   Generalized, default preprocessing function to process waveform data.
   Preprocessing is slightly different for obs and syn waveforms. Each
   processing function is split into a separate function so that they can
   be called by custom preprocessing functions.

   :type st: obspy.core.stream.Stream
   :param st: Stream object to be preprocessed
   :type inv: obspy.core.inventory.Inventory
   :choice inv: Inventory containing response information for real waveform
       data. If not provided, reseponse will not be removed. Also used for
       stream rotation from ZNE -> RTZ
   :type choice: str
   :param choice: 'obs' or 'syn' for observed or synthetic waveform. 'obs' will
       attempt to remove response information with `inv`. 'syn' will attempt
       to convolve with a half duration.
   :type remove_response: bool
   :param remove_response: flag, remove instrument response from choice=='obs'
       data using the provided `inv`. Defaults to True
   :type apply_filter: bool
   :param apply_filter: flag, filter the waveforms using the
       `min_period` and `max_period` parameters. Defaults to True
   :type convolve_with_stf: bool
   :param convolve_with_stf: flag, convolve choice=='syn' data with a Gaussian
       source time function if a `half_dur` (half duration) is provided.
       Defaults to true
   :type rotate_to_rtz: bool
   :param rotate_to_rtz: flag, use the `rotate_baz` variable to rotate streams
       from ZNE components to RTZ
   :rtype: obspy.core.stream.Stream
   :return: preprocessed stream object pertaining to `choice`

   :keyword int water_level: water level for response removal
   :keyword float taper_percentage: amount to taper ends of waveform


.. py:function:: filters(st, min_period=None, max_period=None, min_freq=None, max_freq=None, corners=2, zerophase=True, **kwargs)

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


.. py:function:: taper_time_offset(st, taper_percentage=0.05, time_offset_sec=0)

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


.. py:function:: zero_pad(st, pad_length_in_seconds, before=True, after=True)

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


.. py:function:: trim_streams(st_a, st_b, precision=0.001, force=None)

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


.. py:function:: match_npts(st_a, st_b, force=None)

   Resampling can cause sample number differences which will lead to failure
   of some preprocessing or processing steps. This function ensures that `npts`
   matches between traces by extending one of the traces with zeros.
   A small taper is applied to ensure the new values do not cause
   discontinuities.

   .. note:: its assumed that all traces within a single stream have the same `npts`

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


.. py:function:: is_preprocessed(st, filter_only=True)

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


.. py:function:: stf_convolve(st, half_duration, source_decay=4.0, time_shift=None, time_offset=None)

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


