"""
functions for calculating values for use in other functions
"""
import numpy as np


def myround(x, base=5, choice='near'):
    """round value x to nearest base, round 'up','down' or to 'near'est base
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


def overlapping_days(origin_time, startpad=20, endpad=200):
    """helper function to return a list of julian days based on a given
    origin time with a specific padding on either side. used to catch if an
    origin time sits too close to midnight and two days need to be fetched
    """
    if (origin_time - startpad).julday != origin_time.julday:
        return [(origin_time-startpad).julday, origin_time.julday]
    elif (origin_time + endpad*2).julday != origin_time.julday:
        return [origin_time.julday, (origin_time+endpad*2).julday]
    else:
        return [origin_time.julday]