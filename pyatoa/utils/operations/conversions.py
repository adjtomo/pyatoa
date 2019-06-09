"""
File format conversion functions
"""
import os
import numpy as np
from obspy import Stream, Trace


def read_fortran_binary(path):
    """
    convert a fortran .bin file into a NumPy array,
    stolen from fortran_binary.py _read() in SeisFlows
    
    :type path: str
    :param path: path to fortran .bin file
    :rtype: np.array
    :return: fortran binary data as a numpy array
    """
    nbytes = os.path.getsize(path)
    with open(path, "rb") as f:
        f.seek(0)
        n = np.fromfile(f, dtype="int32", count=1)[0]
        if n == nbytes - 8:
            f.seek(4)
            data = np.fromfile(f, dtype="float32")
            return data[:-1]
        else:
            f.seek(0)
            data = np.fromfile(f, dtype="float32")
            return data


def ascii_to_mseed(path, origintime, location=''):
    """
    Specfem3D outputs seismograms to ASCII (.sem) files
    Pyatoa expects seismograms as obspy Stream objects. 
    This convenience function converts the .sem files into Stream objects
    with the correct header information.    

    :type path: str
    :param path: path of the given ascii file
    :type origintime: obspy.UTCDateTime
    :param origintime: UTCDatetime object for the origintime of the event
    :type location: str
    :param location: location value for a given station/component
    :rtype st: obspy.Stream.stream
    :return st: stream containing header and data info taken from ascii file
    """
    time = np.loadtxt(fname=path, usecols=0)
    data = np.loadtxt(fname=path, usecols=1)
    delta = round(time[1]-time[0], 3)  # assume dt constant after 3 dec. points

    origintime += time[0]  # specfem doesn't start exactly on 0, honor that
    net, sta, cha, fmt = os.path.basename(path).split('.')
    stats = {"network": net, "station": sta, "location": location,
             "channel": cha, "starttime": origintime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'},
             "time_offset": time[0], "format":fmt
             }
    st = Stream([Trace(data=data, header=stats)])

    return st


def mt_transform(mt, method):
    """
    Transform moment tensor between XYZ and RTP coordinates

    Acceptable formats for the parameter mt:
        1) [m11,m22,m33,m12,m13,m23]
        2) [mxx,myy,mzz,mxy,mxz,myz]
        3) [mrr,mtt,mpp,mrt,mrp,mtp]
    
    Based on equation ?? from Aki and Richards Quantitative Seismology
    TO DO: find the correct equation number

    :type mt: dict
    :param mt: moment tensor in format above
    :type method: str
    :param method: type of conversion, "rtp2xyz" or "xyz2rtp"
    :rtype: dict
    :return: converted moment tensor dictionary
    """
    if method == "xyz2rtp":
        if "m_xx" not in mt.keys():
            print("for xyz2rtp, dict must have keys in xyz")
        m_rr = mt["m_zz"]
        m_tt = mt["m_xx"]
        m_pp = mt["m_yy"]
        m_rt = mt["m_xz"]
        m_rp = -1 * mt["m_yz"]
        m_tp = -1 * mt["m_xy"]
        return {"m_rr": m_rr, "m_tt": m_tt, "m_pp": m_pp, "m_rt": m_rt,
                "m_rp": m_rp, "m_tp": m_tp}

    if method == "rtp2xyz":
        if "m_tt" not in mt.keys():
            print("for rtp2xyz, dict must have keys in rtp")
        m_xx = mt["m_tt"]
        m_yy = mt["m_pp"]
        m_zz = mt["m_rr"]
        m_xy = -1 * mt["m_tp"]
        m_xz = mt["m_rt"]
        m_yz = -1 * mt["m_rp"]
        return {"m_xx": m_xx, "m_yy": m_yy, "m_zz": m_zz, "m_xy": m_xy,
                "m_xz": m_xz, "m_yz": m_yz}
    else:
        print("Invalid transformation method, xyz2rtp or rtp2xyz")
        return None

