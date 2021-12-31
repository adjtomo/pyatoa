"""
Formatting functionality

Pyatoa relies on data structure being ordered and consistent throughout all the
various bits of data required, as well as a few standardized string formatters
to keep everything playing nice. Functions here will aid in reshaping data
into the correct formats.
"""
import os
from pyasdf import ASDFDataSet
from obspy.core.event import Event
from pyatoa.utils.srcrcv import Source


def format_iter(iteration):
    """
    Format the iteration to be used in internal path naming and labelling.
    Standard is to format with a leading 'i' character followed by two digits. 
    Inputs can be strings or integers. Assumes iterations won't go above 99.

    :type iteration: str or int
    :param iteration: input model number, e.g. 0, '0', 'i0', 'i00'
    :rtype: str
    :return: formatted model number, e.g. 'i00', or None if no matching format
    """
    if isinstance(iteration, str):
        # If e.g. iteration = "0"
        if not iteration[0] == "i":
            iternum = f"i{iteration:0>2}"
        # If e.g. iteration = "m00"
        else:
            iternum = iteration
    # If e.g. iteration = 0
    elif isinstance(iteration, int):
        iternum = f"i{iteration:0>2}"
    else:
        iternum = None

    return iternum


def format_step(count):
    """
    Same as for iteration but step count is formatted with a leading 's'

    :type count: str or int
    :param count: input model number, e.g. 0, '0', 's0', 's00'
    :rtype: str
    :return: formatted model number, e.g. 's00', or None if no matching format
    """
    if isinstance(count, str):
        if not count[0] == "s":
            stpcnt = f"s{count:0>2}"
        else:
            stpcnt = count
    elif isinstance(count, int):
        stpcnt = f"s{count:0>2}"
    else:
        stpcnt = None

    return stpcnt


def format_event_name(ds_or_event):
    """
    Formalize the definition of Event ID in Pyatoa

    :type ds_or_event: pyasdf.ASDFDataSet or obspy.core.event.Event or str
    :param ds_or_event: get dataset event name from the event resource_id
    :rtype: str
    :return: the event name to be used for naming schema in the workflow
    """
    # Extract the resource id dependent on the input file type
    if isinstance(ds_or_event, (Event, Source)):
        rid = ds_or_event.resource_id.id
    elif isinstance(ds_or_event, str):
        rid = ds_or_event
    elif isinstance(ds_or_event, ASDFDataSet):
        rid = ds_or_event.events[0].resource_id.id
    else:
        raise TypeError("format_event_name() only accepts pyasdf.ASDFDataSet "
                        "or obspy.core.event.Event objects")

    rid_up = rid.upper()
    # GeoNet Client: smi:nz.org.geonet/2018p130600
    if "GEONET" in rid_up:
        return rid.split("/")[-1]
    # IRIS Client: smi:service.iris.edu/fdsnws/event/1/query?eventid=5197722
    elif "IRIS" in rid_up:
        return rid.split("eventid=")[-1]
    # SPUD, GCMT: smi:local/ndk/C202005010101A/event
    # or CMTSOLUTION: smi:local/cmtsolution/2013p617227/event
    elif ("NDK" in rid_up) or ("CMTSOLUTION" in rid_up):
        return rid.split("/")[-2]
    # USGS Client: quakeml:earthquake.usgs.gov/fdsnws/event/1/...
    #                          ...query?eventid=ak0198gdhuwa&format=quakeml
    elif "USGS" in rid_up:
        rid_ = rid.split("eventid=")[-1]
        return rid_.split("&")[0]
    # USGS ANSS ComCat: quakeml:us.anss.org/event/20005ysu
    elif "ANSS" in rid_up:
        return rid.split("event/")[-1]
    # PYATOA pyatoa.read.read_force_solution: pyatoa:source/
    elif "PYATOA" in rid_up:
        return rid.split("source/")[-1]
    else:
        raise NotImplementedError(f"Unknown resource id format {rid}, "
                                  "Please raise a GitHub issue and the"
                                  "developers will address this")


def channel_code(dt):
    """
    Specfem outputs seismograms with channel band codes set by IRIS. Instrument
    codes are always X for synthetics, but band code will vary with the sampling
    rate of the data, return the correct code given a sampling rate.
    Taken from Appenix B of the Specfem3D cartesian manual (June 15, 2018)

    :type dt: float
    :param dt: sampling rate of the data in seconds
    :rtype: str
    :return: band code as specified by Iris
    """
    if dt >= 1:
        return "L"  # long period
    elif 0.1 < dt < 1:
        return "M"  # mid period
    elif 0.0125 < dt <= 0.1:
        return "B"  # broad band
    elif 0.001 <= dt <= 0.0125:
        return "H"  # high broad band
    elif 0.004 <= dt < 0.001:
        return "C"
    elif 0.001 <= dt < 0.0002:
        return "F"
    else:
        print("Channel code does not exist for this value of 'dt'")
        return None

