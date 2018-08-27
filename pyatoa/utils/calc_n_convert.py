"""file format conversion functions and calculation functinos
"""
import os
import numpy as np
from obspy import Stream, Trace

def myround(x,base=5,choice='near'):
    """round value x to nearest base, round 'up','down' or to 'near'est base
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout

def ascii_to_mseed(event,ascii_path,mseed_path=None):
    """convert native output of specfem into mseed format
    :type event: obspy.event
    :param event: specific event information
    """
    starttime = event.origins[0].time
    time = np.loadtxt(fname=ascii_path,usecols=0)
    data = np.loadtxt(fname=ascii_path,usecols=1)
    delta = round(time[1]-time[0],3) #assume dt is constant after 3 dec. points

    net,sta,cha,cmp = os.path.basename(ascii_path).split('.')
    stats = {"network":network,"station":station,"location":location,
        "channel":channel,"starttime":starttime,"npts":len(data),"delta":delta,
        "mseed":{"dataquality":'D'}
             }
    st = Stream([Trace(data=data,header=stats)])
    if mseed_path:
        st.write(os.path.join(mseed_path,os.path.basename(ascii_path)+'.mseed'),
            format="MSEED")

    return st

def caculate_julian_day(origin_time,startpad=20,endpad=200):
    """helper function to return a list of julian days based on a given
    origin time with a specific padding on either side. used to catch if an
    origin time sits too close to midnight and two days need to be fetched
    """
    julday = origin_time.julday
    # check if beginning of waveform too close to previous day
    if (origin_time - pad) != julday:
        return [(origin_time-pad).julday,julday]
    elif (origin_time + pad*2) != julday:
        return [julday,(origin_time+pad*2).julday]
    else:
        return [julday]
