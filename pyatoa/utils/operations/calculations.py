"""
Functions for calculating values for use in other functions
"""
import numpy as np


def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


def overlapping_days(origin_time, start_pad=20, end_pad=200):
    """
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
    """
    if (origin_time - start_pad).julday != origin_time.julday:
        return [(origin_time-start_pad).julday, origin_time.julday]
    elif (origin_time + end_pad*2).julday != origin_time.julday:
        return [origin_time.julday, (origin_time+end_pad*2).julday]
    else:
        return [origin_time.julday]


def normalize_a_to_b(array, a=0, b=1):
    """
    normalize an array from a to b for easy plotting

    :type array: list
    :param array: values to be normalized
    :type a: int
    :param a: lower bound of normalization
    :type b: int
    :param b: upper bound of normalization
    :rtype z: numpy.array
    :return z: normalized array
    """
    array = np.array(array)
    z = ((b-a) * (array-array.min()) / (array.max()-array.min())) + a

    return z
