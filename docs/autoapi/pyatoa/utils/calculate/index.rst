:py:mod:`pyatoa.utils.calculate`
================================

.. py:module:: pyatoa.utils.calculate

.. autoapi-nested-parse::

   Custom math functions for faster calculations in other parts of Pyatoa



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.calculate.abs_max
   pyatoa.utils.calculate.myround
   pyatoa.utils.calculate.enforce_angle_pi
   pyatoa.utils.calculate.overlapping_days
   pyatoa.utils.calculate.normalize_a_to_b
   pyatoa.utils.calculate.amplitude_anomaly
   pyatoa.utils.calculate.vrl
   pyatoa.utils.calculate.time_offset



.. py:function:: abs_max(array)

   Find the absolute maximum of an array

   :type array: np.array
   :param array: array to find abs max of
   :rtype: float
   :return: absolute maximum of array


.. py:function:: myround(x, base=5, choice='near')

   Round value x to nearest base, round 'up','down' or to 'near'est base

   :type x: float
   :param x: value to be rounded
   :type base: int
   :param base: nearest integer to be rounded to
   :type choice: str
   :param choice: method of rounding, 'up', 'down' or 'near'
   :rtype roundout: int
   :return: rounded value


.. py:function:: enforce_angle_pi(angle, bound=180)

   Used for mapping mostly, sometimes when adding values to map bounds, you
   can cross over the accepted threshold, e.g., angle > 180deg.
   This function just keeps things within bounds

   https://stackoverflow.com/questions/2320986/
                         easy-way-to-keeping-angles-between-179-and-180-degrees


.. py:function:: overlapping_days(origin_time, start_pad=20, end_pad=200)

   Helper function to return a list of julian days based on a given
   origin time with a specific padding on either side. used to catch if an
   origin time sits too close to midnight and two days need to be fetched

   :type origin_time: obspy.core.UTCDateTime
   :param origin_time: event origin time
   :param start_pad: padding in seconds before the origin time of an event
       for waveform fetching, to be fed into lower level functions.
   :type end_pad: int
   :param end_pad: padding in seconds after the origin time of an event
       for wavefomr fetching.
   :rtype: list of int
   :return: list of available julian days


.. py:function:: normalize_a_to_b(array, a=0, b=1)

   normalize an array from a to b for e.g. plotting, maths

   :type array: list
   :param array: values to be normalized
   :type a: int
   :param a: lower bound of normalization
   :type b: int
   :param b: upper bound of normalization
   :rtype z: numpy.array
   :return z: normalized array


.. py:function:: amplitude_anomaly(a, b, dt)

   Calculate the amplitude differences between two waveforms, a la.
   Equation A2 from Tape et al. 2010, which states that:

   DlnA = ln(a/b) = 0.5 * ln[integral(a(t)**2 dt)/integral(b(t)**2 dt)]
   where a and b represent data and synthetics, respectively

   .. note::
       It is expected that a and b have the same value of dt, if they do not,
       they should be resampled before being passed to this function.

   :type a: np.array
   :param a: waveform data to act as numerator of misfit definition
   :type b: np.array
   :param b: waveform data to act as denominator of misfit definition
   :type dt: float
   :param dt: sampling rate for the integration
   :rtype: float
   :return: the value of DlnA, the amplitude anomaly


.. py:function:: vrl(d, s1, s2)

   Logarithmic variance reduction. Approximately Gaussian distributed
   reduction in variance based on full waveform difference. Based on
   Equation 8 in Tape et al. 2010.

   :type d: np.array
   :param d: data array to compare against synthetic arrays
   :type s1: np.array
   :param s1: synthetic array for model 1
   :type s2: np.array
   :param s2: synthetic array for model 2
   :rtype: float
   :return: logarithmic variance reduction


.. py:function:: time_offset(st, origintime)

   Oft repeated calculation of finding the time offset between the start of
   a stream and the origin time of an event. Important when dealing with
   SPECFEM seismograms, which explicitely start at the time offset value,
   w.r.t. streams which set starttime at 0 but contain their own
   origin time information

   .. note::
       convention is that if the waveform starts BEFORE the event origin time,
       then `time_offset` will be NEGATIVE.

   :type st: obspy.core.stream.Stream
   :param st: stream to check time offset of
   :type origintime: obspy.UTCDateTime
   :param origintime: the origin time of the event
   :rtype: float
   :return: the time offset, or: origintime - stream_origin_time


