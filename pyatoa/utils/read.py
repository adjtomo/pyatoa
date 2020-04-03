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
from obspy.core.event import CreationInfo, Event, Catalog, EventDescription
from obspy.core.event.origin import Origin
from obspy.core.event.magnitude import Magnitude
from obspy.core.event.source import Tensor, MomentTensor, FocalMechanism


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


def read_ascii(path, origintime=None, location=''):
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
    # This was tested up to version 6895e2f7
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

    # We assume that dt is constant after 3 decimal points
    delta = round(times[1] - times[0], 3)

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

    # Get all the unique network names
    networks = {_: [] for _ in np.unique(stations[:, 1])}
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


def read_cmtsolution(path_to_cmtsolution, rtype="event"):
    """
    Convert a Specfem3D CMTSOLUTION file into an ObsPy Event object

    Note: ResourceID's are not handled, they will be auto-set by ObsPy.

    The values in the CMTSOLUTION are expected to be (in order):
     event_name, time_shift, half_duration, latitude, longitude, depth,
     Mrr, Mtt, Mpp, Mrt, Mrp, Mtp

    The header should have the format:
     pde, year, month, day, hour, minute, second, lat, lon, depth, mb, ms, name

     To Do:
        Convert moment tensor to strike dip rake object

    :type path_to_cmtsolution: str
    :param path_to_cmtsolution: path to the CMTSOLUTION file associated with
        a Specfem3D DATA directory
    :type rtype: str
    :param rtype: return type, choice between 'Event' and 'Catalog' objects
        defaults to returning an Event object
    :rtype: obspy.core.event.Event or obspy.core.event.Catalog
    :return: converted CMTSOLUTION into an Event, or a Catalog object.
    """
    assert(rtype in ["event", "catalog"]), \
        "Return type must be 'event' or 'catalog'"

    def seismic_moment(mt):
        """
        Return the seismic moment based on a moment tensor. Can take a list of
        tensor components, or a Tensor object from ObsPy.

        Same as pyatoa.utils.srcrcv.seismic_moment()

        :type mt: list of floats or obspy.core.event.source.Tensor
        :param mt: the components of the moment tensor M_ij
        :rtype: float
        :return: the seismic moment, in units of N*m
        """
        if isinstance(mt, Tensor):
            # Little one liner to spit out moment tensor components into a list
            mt_temp = [getattr(mt, key) for key in mt.keys()
                                            if not key.endswith("errors")]
            assert(len(mt_temp) == 6), "Moment tensor should have 6 components"
            mt = mt_temp
        return 1/np.sqrt(2) * np.sqrt(sum([_**2 for _ in mt]))

    def moment_magnitude(moment):
        """
        Return the moment magitude based on a seismic moment, from
        Hanks & Kanamori (1979)

        Same as pyatoa.utils.srcrcv.moment_magnitude()

        :type moment: float
        :param moment: the seismic moment, in units of N*m
        :rtype: float
        :return: moment magnitude, M_w
        """
        return 2/3 * np.log10(moment) - 10.7

    # Read in header and body
    header = np.genfromtxt(path_to_cmtsolution, dtype="str", max_rows=1)
    cmtsolution = np.genfromtxt(path_to_cmtsolution, dtype="str", skip_header=1,
                                delimiter=":")

    # Parse header, get origin time and name
    year = header[1]
    month = header[2]
    day = header[3]
    hour = header[4]
    minute = header[5]
    seconds = header[6]
    origintime = UTCDateTime(f"{year}-{month}-{day}T{hour}:{minute}:{seconds}")
    name = header[-1]  # a text description of where the event occurred

    # Parse the body of the CMTSOLUTION, event ID gets separate treatment
    # Replace spaces in tags with underscores
    cmt = dict()
    cmt[cmtsolution[0][0].replace(" ", "_")] = cmtsolution[0][1].strip()
    for value in cmtsolution[1:]:
        cmt[value[0].replace(" ", "_")] = float(value[1].strip())

    # Tensor to be put in the MomentTensor
    tensor = Tensor(m_rr=cmt["Mrr"], m_tt=cmt["Mtt"], m_pp=cmt["Mpp"],
                    m_rt=cmt["Mrt"], m_rp=cmt["Mrp"], m_tp=cmt["Mtp"]
                    )

    # Use the tensor to get M0 and Mw
    scalar_moment = seismic_moment(tensor)
    mw = moment_magnitude(scalar_moment)

    # Create various objects used to fill the Event to be more info rich
    origin = Origin(time=origintime, longitude=cmt["longitude"],
                    latitude=cmt["latitude"], depth=cmt["depth"]*1E3,
                    )
    moment_tensor = MomentTensor(tensor=tensor, scalar_moment=scalar_moment)
    focal_mechanism = FocalMechanism(moment_tensor=moment_tensor)
    magnitude = Magnitude(mag=float(f"{mw:.2f}"), magnitude_type="Mw")
    info = CreationInfo(author="Pyatoa from CMTSOLUTION",
                        creation_time=UTCDateTime())
    description = EventDescription(text=name)

    # Create the Event using all the pieces created, set preferences
    event = Event(origins=[origin], focal_mechanisms=[focal_mechanism],
                  magnitudes=[magnitude], event_descriptions=[description],
                  creation_info=info)

    event.preferred_origin_id = origin.resource_id
    event.preferred_focal_mechanism_id = focal_mechanism.resource_id
    event.preferred_magnitude_id = magnitude.resource_id

    # Determine if an Event or Catalog is returned
    if rtype == "event":
        return event
    elif rtype == "catalog":
        return Catalog(events=[event])

