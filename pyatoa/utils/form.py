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
    Formalize the defition of Event ID in Pyatoa

    :type ds_or_event: pyasdf.ASDFDataSet or obspy.core.event.Event
    :param ds_or_event: get dataset event name from the filename
    :type event: obspy.core.event.Event
    :param event: get event name from the resource id
    :rtype: str
    :return: the event name to be used for naming schema in the workflow
    """
    if isinstance(ds_or_event, Event):
        return ds_or_event.resource_id.id.split('/')[-1]
    elif isinstance(ds_or_event, ASDFDataSet):
        return os.path.basename(ds_or_event.filename).split(".")[0]
    else:
        raise TypeError("format_event_name() only accepts pyasdf.ASDFDataSet "
                        "or obspy.core.event.Event objects")


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

