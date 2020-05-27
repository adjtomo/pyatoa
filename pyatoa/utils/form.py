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


def format_model_number(iteration):
    """
    The model number is based on the current iteration and should be formatted
    with a leading 'm' character and two digits. Inputs can be strings or
    integers. Assuming model numbers won't go above 99.

    :type iteration: str or int
    :param iteration: input model number, e.g. 0, '0', 'm0', 'm00'
    :rtype: str
    :return: formatted model number, e.g. 'm00', or None if no matching format
    """
    if isinstance(iteration, str):
        # If e.g. model_number = "0"
        if not iteration[0] == "m":
            mdlnmbr = f"m{iteration:0>2}"
        # If e.g. model_number = "m00"
        else:
            mdlnmbr = iteration
    # If e.g. model_number = 0
    elif isinstance(iteration, int):
        mdlnmbr = f"m{iteration:0>2}"
    else:
        mdlnmbr = None

    return mdlnmbr


def format_step_count(count):
    """
    Same as for model number but step count is formatted with a leading 's'

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


def distribute_dataless(path_to_response, inv):
    """
    Response files written through Obspy come out as a single object, but Pyatoa
    will look for response information from individual components and individual
    stations. Distrubute this dataless information into the necessary components

    Note: The template here is hardcoded with SEED convention

    :type path_to_response: str
    :param path_to_response: path to save the new response files to
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory with response to be delinieated into separate objects
    """
    inner_folder = "{STA}.{NET}"
    fid_template = "RESP.{NET}.{STA}.{LOC}.{CHA}"
    full_template = os.path.join(path_to_response, inner_folder, fid_template)
    for net in inv:
        for sta in net:
            try:
                # Create the container directory unless it exists
                os.mkdir(os.path.join(path_to_response, inner_folder.format(
                    STA=sta.code, NET=net.code))
                    )
            except FileExistsError:
                pass
            for cha in sta:
                # Write the individual channel as a STATIONXML file
                inv_temp = inv.select(network=net.code, station=sta.code,
                                      location=cha.location_code,
                                      channel=cha.code)
                inv_temp.write(full_template.format(
                    STA=sta.code, NET=net.code, LOC=cha.location_code,
                    CHA=cha.code),
                    format='STATIONXML'
                    )


def create_window_dictionary(window):
    """
    HDF5 doesnt play nice with nonstandard objects in dictionaries, e.g.
    nested dictionaries, UTCDateTime objects. So remake the Pyflex window
    JSON dictionary into something that will sit well in a Pyasdf object

    :type window: pyflex.Window
    :param window: misfit window calcualted by pyflex
    :rtype win_dict: dict
    :return win_dict: formatted dictionary of values for pyasdf auxiliary data
    """
    win_dict = window._get_json_content()

    # change UTCDateTime objects into strings
    win_dict["absolute_endtime"] = str(win_dict["absolute_endtime"])
    win_dict["absolute_starttime"] = str(win_dict["absolute_starttime"])
    win_dict["time_of_first_sample"] = str(win_dict["time_of_first_sample"])

    phase_arrivals = win_dict["phase_arrivals"]
    for phase in phase_arrivals:
        win_dict["phase_arrival_{}".format(phase["name"])] = phase["time"]
    win_dict.pop("phase_arrivals")

    return win_dict


