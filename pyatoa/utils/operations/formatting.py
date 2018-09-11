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


def moment_tensor_list_to_objects(mtlist):
    """
    UNFINISHED
    event objects fetched by obspy do not natively come with any moment tensor
    or nodal plane information, that is stored separately in a .csv file
    located on github. For the tomography problem we need this information, so
    this functino will append information from the .csv file onto the obspy
    event object so that all the information can be located in a single object
    :param event:
    :param geonet_moment_tensor_list:
    :return:
    """
    from obspy.core.event import ResourceIdentifier
    import obspy.core.event.source as eventcore

    id_template = "smi:local/geonetcsv/{0}/{1}".format(mtlist['PublicID'],'{}')
    if len(mtlist) != 32:
        print("geonet moment tensor list does not have the correct number"
              "of requisite components, should have 32")
        return
    nodal_plane_1 = eventcore.NodalPlane(strike=mtlist['strike1'],
                                         dip=mtlist['dip1'],
                                         rake=mtlist['rake1']
                                         )
    nodal_plane_2 = eventcore.NodalPlane(strike=mtlist['strike2'],
                                         dip=mtlist['dip2'],
                                         rake=mtlist['rake2']
                                         )
    nodal_planes = eventcore.NodalPlanes(nodal_plane_1, nodal_plane_2,
                                         preferred_plane=1)
    moment_tensor = eventcore.MomentTensor(
        resource_id=id_template.format('momenttensor'),
        derived_origin_id=id_template.format('origin#ristau'),
        scalar_moment=mtlist['Mo']*1E-7
        )

    focal_mechanism = eventcore.FocalMechanism(
        resource_id=ResourceIdentifier(
            "smi:local/geonetcsv/{}/focal_mechanism".format(mtlist['PublicID'])
            ),
        nodal_planes=nodal_planes, moment_tensor=moment_tensor
        )