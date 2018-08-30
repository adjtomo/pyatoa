"""file format conversion functions
"""
import os
import numpy as np
from obspy import Stream, Trace


def ascii_to_mseed(path, origintime, location=''):
    """convert native output of specfem into mseed format
    :type event: obspy.event
    :param event: specific event information
    """
    time = np.loadtxt(fname=path,usecols=0)
    data = np.loadtxt(fname=path,usecols=1)
    delta = round(time[1]-time[0],3) #assume dt is constant after 3 dec. points

    net, sta, cha, fmt = os.path.basename(path).split('.')

    stats = {"network": net, "station": sta, "location": location,
             "channel": cha, "starttime": origintime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'}
             }
    st = Stream([Trace(data=data,header=stats)])

    return st

