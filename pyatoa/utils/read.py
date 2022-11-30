"""
Utilities for reading various file types, mostly from Specfem3D to ObsPy classes
These are meant to be standalone functions so they may repeat some functionality
found elsewhere in the package.
"""
import os
import numpy as np
from obspy import Stream, Trace, UTCDateTime, Inventory
from obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station
from obspy.core.event import Event, Origin, Magnitude
from pyatoa import logger


def read_fortran_binary(path):
    """
    Convert a Specfem3D fortran .bin file into a NumPy array,
    Copied verbatim from Seisflows/plugins/solver_io/fortran_binary.py/_read()

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



def read_station_codes(path_to_stations, loc="??", cha="*",
                       seed_template="{net}.{sta}.{loc}.{cha}"):
    """
    Read the SPECFEM3D STATIONS file and return a list of codes (Pyatoa format)
    that are accepted by the Manager and Pyaflowa classes. Since the STATIONS
    file only provides NET and STA information, the user must provide the
    location and channel information, which can be wildcards.

    :type path_to_stations: str
    :param path_to_stations: full path to the STATIONS file
    :type loc: str
    :param loc: formatting of the location section of the code, defaults to
        '??' two-digit wildcard
    :type cha: str
    :param cha: formatting of the channel section fo the code, defaults to
        'HH?' for wildcard component of a high-gain seismometer. Follows SEED
        convention (see IRIS).
    :type seed_template: str
    :param seed_template: string template to be formatted with some combination
        of 'net', 'sta', 'loc' and 'cha', used for generating station codes
    :rtype: list of str
    :return: list of codes to be used by the Manager or Pyaflowa classes for
        data gathering and processing
    """
    codes = []
    stations = np.loadtxt(path_to_stations, dtype="str")
    if stations.size == 0:
        return codes

    # Deal with the special case where the stations file is only 1 station long
    # otherwise we end up iterating over a string and not an ndarray
    if not isinstance(stations[0], np.ndarray):
        stations = [stations]

    for station in stations:
        sta = station[0]
        net = station[1]
        codes.append(seed_template.format(net=net, sta=sta, loc=loc, cha=cha))

    return codes


