"""
pyatoa relies on data structure being ordered and consistent throughout all the
various bits of data required. functions here will aid in reshaping data
into the correct formats
"""
import os


def distribute_dataless(path_to_response,inventory):
    """
    Response files written through obspy come out as a single object, but pyatoa
    will look for response information from individual components and individual
    stations. Distrubute this dataless information into the necessary components
    """
    inner_folder = '{STA}.{NET}'
    fid_template = 'RESP.{NET}.{STA}.{LOC}.{CHA}'
    full_template = os.path.join(path_to_response,inner_folder,fid_template)
    for net in inv:
        for sta in net:
            try:
                os.mkdir(os.path.join(path_to_response,inner_folder.format(
                    STA=sta.code, NET=net.code))
                    )
            except FileExistsError:
                pass
            for cha in sta:
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
    nested dictionaries, UTCDateTime objects. So remake the pyflex window
    json dictionary into something that will sit well in a pyasdf object
    """
    win_dict = window._get_json_content()

    # change UTCDateTime objects into strings
    win_dict['absolute_endtime'] = str(win_dict['absolute_endtime'])
    win_dict['absolute_starttime'] = str(win_dict['absolute_starttime'])
    win_dict['time_of_first_sample'] = str(win_dict['time_of_first_sample'])

    phase_arrivals = win_dict['phase_arrivals']
    for phase in phase_arrivals:
        win_dict['phase_arrival_{}'.format(phase['name'])] = phase['time']
    win_dict.pop('phase_arrivals')

    return win_dict
