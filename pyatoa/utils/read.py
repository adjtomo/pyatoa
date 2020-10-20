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


def read_sem(path, origintime=None, location='', precision=4):
    """
    Specfem3D outputs seismograms to ASCII (.sem?) files

    Pyatoa expects seismograms as Obspy Stream objects.
    This convenience function converts the .sem? files into Stream objects
    with the correct header information.

    Works with Specfem3D git version 6895e2f7

    Note: if origintime is None, the default start time is 1970-01-01T00:00:00,
          this is fine for quick-looking the data but it is recommended that an
          actual starttime is given.

    :type path: str
    :param path: path of the given ascii file
    :type origintime: obspy.UTCDateTime
    :param origintime: UTCDatetime object for the origintime of the event
    :type location: str
    :param location: location value for a given station/component
    :rtype st: obspy.Stream.stream
    :return st: stream containing header and data info taken from ascii file
    """
    # This was tested up to SPECFEM3D Cartesian git version 6895e2f7
    try:
        times = np.loadtxt(fname=path, usecols=0)
        data = np.loadtxt(fname=path, usecols=1)

    # At some point in 2018, the Specfem developers changed how the ascii files
    # were formatted from two columns to comma separated values, and repeat
    # values represented as 2*value_float where value_float represents the data
    # value as a float
    except ValueError:
        times, data = [], []
        with open(path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            try:
                time_, data_ = line.strip().split(',')
            except ValueError:
                if "*" in line:
                    time_ = data_ = line.split('*')[-1]
                else:
                    raise ValueError
            times.append(float(time_))
            data.append(float(data_))

        times = np.array(times)
        data = np.array(data)

    if origintime is None:
        print("No origintime given, setting to default 1970-01-01T00:00:00")
        origintime = UTCDateTime("1970-01-01T00:00:00")

    # We assume that dt is constant after 'precision' decimal points
    delta = round(times[1] - times[0], precision)

    # Honor that Specfem doesn't start exactly on 0
    origintime += times[0]

    # Write out the header information
    net, sta, cha, fmt = os.path.basename(path).split('.')
    stats = {"network": net, "station": sta, "location": location,
             "channel": cha, "starttime": origintime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'},
             "time_offset": times[0], "format": fmt
             }
    st = Stream([Trace(data=data, header=stats)])

    return st


def read_stations(path_to_stations):
    """
    Convert a Specfem3D STATIONS file into an ObsPy Inventory object.

    Specfem3D STATION files contain no channel or location information, so
    the inventory can only go down to the station level.

    Note:
        This assumes a row structure for the station file is
        STA, NET, LAT [deg], LON [deg], ELEVATION [m], BURIAL [m]

    :type path_to_stations: str
    :param path_to_stations: the path to the STATIONS file that is associated
        with the Specfem3D DATA directory
    :rtype: obspy.core.inventory.Inventory
    :return: a station-level Inventory object
    """
    stations = np.loadtxt(path_to_stations, dtype="str")

    # Get all the unique network names, try-except to catch when there is only
    # one station in the file
    try:
        networks = {_: [] for _ in np.unique(stations[:, 1])}
    except IndexError:
        networks = {stations[1]: []}
        stations = [stations]

    for sta in stations:
        # Parse the station information
        station_ = sta[0]
        network_ = sta[1]
        latitude_ = float(sta[2])
        longitude_ = float(sta[3])
        elevation_ = float(sta[4])
        burial_ = float(sta[5])  # burial isnt an option in ObsPy

        # Create the station object, temp store in a network
        station = Station(code=station_, latitude=latitude_,
                          longitude=longitude_, elevation=elevation_,
                          creation_date=UTCDateTime()
                          )
        networks[network_].append(station)

    # Create the network objects
    list_of_networks = []
    for network, stations in networks.items():
        list_of_networks.append(Network(code=network, stations=stations))

    return Inventory(networks=list_of_networks, source="PYATOA")


def read_station_codes(path_to_stations, loc="??", cha="HH?"):
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
    :rtype: list of str
    :return: list of codes to be used by the Manager or Pyaflowa classes for
        data gathering and processing
    """
    SEED_TEMPLATE = "{net}.{sta}.{loc}.{cha}"

    codes = []
    stations = np.loadtxt(path_to_stations, dtype="str")
    for station in stations:
        sta = station[0]
        net = station[1] 
        codes.append(SEED_TEMPLATE.format(net=net, sta=sta, loc=loc, cha=cha)

    return codes


def read_specfem_vtk(path_to_vtk):
    """
    Read the unstructured grid VTK files that are output by Specfem for model,
    gradient, and kernel visualizations. Returns a header as a dictionary, and
    the lines of the data file. Useful for manipulating VTK files in place.

    :type path_to_vtk: str
    :param path_to_vtk: full path to the .vtk file to read
    :rtype lines: list
    :return lines: data line by line from readlines()
    :rtype header_dict: dic
    :return header_dict: dictionary with all relevant header information
        from an unstructured_grid vtk file
    """
    with open(path_to_vtk, "r") as f:
        lines = f.readlines()

    # determine important line numbers, headers are specific to specfem3d output
    for i, line in enumerate(lines):
        if "POINTS" in line:
            points_n = int(line.split()[1])
            points_line = i + 1
        elif "CELL" in line and "CELL_TYPE" not in line:
            cells_n = int(line.strip().split()[1])
            cells_size = int(line.strip().split()[2])
            cells_line = i
        elif "CELL_TYPES" in line:
            cell_types_n = int(line.strip().split()[1])
            cell_types_line = i
        elif "POINT_DATA" in line:
            point_data_n = int(line.split()[1])
            point_data_line = i
        elif "SCALARS" in line:
            scalars = line.strip().split()[1]
            data_line = i + 2

    # easier returns in a dictionary
    header_dict = {"points_n": points_n, "points_line": points_line,
                   "cells_n": cells_n, "cells_size": cells_size,
                   "cells_line": cells_line, "cell_types_n": cell_types_n,
                   "cell_types_line": cell_types_line,
                   "points_data_n": point_data_n,
                   "points_data_line": point_data_line, "scalars": scalars,
                   "data_line": data_line
                   }

    return lines, header_dict

