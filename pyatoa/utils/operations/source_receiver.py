"""
Functions used in determining source receiver information
"""


def gcd_and_baz(event, inv):
    """
    Calculate great circle distance and backazimuth values for a given
    station and event configuration

    :type event: obspy.core.event.Event
    :param event: event object
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory object, assumed only 1 station which takes up the
        initial index of the inventory
    :rtype gcdist: float
    :return gcdist: great circle distance in km
    :rtype baz: float
    :return baz: backazimuth in degrees
    """
    from obspy.geodetics import gps2dist_azimuth
    gcdist, _, baz = gps2dist_azimuth(lat1=event.preferred_origin().latitude,
                                      lon1=event.preferred_origin().longitude,
                                      lat2=inv[0][0].latitude,
                                      lon2=inv[0][0].longitude
                                      )
    return gcdist*1E-3, baz


def theoretical_p_arrival(source_depth_in_km, distance_in_degrees,
                          model='iasp91'):
    """
    calculate theoretical arrivals
    """
    from obspy.taup import TauPyModel
    model = TauPyModel(model=model)
    return model.get_travel_times(source_depth_in_km, distance_in_degrees,
                                  phase_list=["P"])


def parse_inventory(inv, event=None):
    """
    return information from an inventory object. check data availability based
    on event origin time and return only the available stations
    """
    network_codes, station_codes, latitudes, longitudes, starts, ends =\
        [], [], [], [], [], []
    for net in inv:
        for sta in net:
            if event is not None:
                if sta.is_active(time=event.preferred_origin().time):
                    network_codes.append(net.code)
                    station_codes.append(sta.code)
                    latitudes.append(float(sta.latitude))
                    longitudes.append(float(sta.longitude))
                    # starts.append(sta.start_date)
                    # ends.append(sta.end_date)
            else:
                network_codes.append(net.code)
                station_codes.append(sta.code)
                latitudes.append(float(sta.latitude))
                longitudes.append(float(sta.longitude))

    return network_codes, station_codes, latitudes, longitudes


def half_duration_from_m0(moment):
    """
    Empirical formula for half duration used by Harvard CMT, stated in
    Daheln and Tromp (1998, p.178).

    :type moment: float
    :param moment: seismic moment in N*m
    :rtype: float
    :return: empirically scaled half duration
    """
    return 2.4E-6 * moment**(1/3)


def event_by_distance(cat, filter_type=False, filter_bounds=None, random=False):
    """
    Sort through an obspy catalog by interevent distance. If we have a lot of
    events in a catalog, it's best to spatially vary them such that we don't
    redundantly oversample a spatial region.
    Returns an index list for events that are most distant from one another,
    without repeating any used events.

    example call:
    >> index_list, event_list = event_by_distance(cat, filter_type="magnitude",
                                                  filter_bounds=[5.0,6.0])

    :param cat: obspy.event.Catalog
    :param filter_type: str
    :param filter_bounds: list of floats
    :param random: bool
    :return:
    """
    from obspy.geodetics import locations2degrees
    # filter the catalog by e.g. magnitude before sorting
    if filter_type:
        cat = cat.filter("{f} >= {lb}".format(
            f=filter_type, lb=filter_bounds[0]))
        cat = cat.filter("{f} <= {lb}".format(
            f=filter_type, lb=filter_bounds[1]))

    # determine interevent distances
    latitudes, longitudes = [], []
    for event in cat:
        latitudes.append(event.preferred_origin().latitude)
        longitudes.append(event.preferred_origin().longitude)

    # determine starting point
    if random:
        from random import randint
        starting_index = randint(0, len(cat)-1)
    else:
        starting_index = 0
    lat = cat[starting_index].preferred_origin().latitude
    lon = cat[starting_index].preferred_origin().longitude

    # create list
    index_list = [starting_index]
    while len(index_list) < len(cat):
        dist_list = []
        for i in range(len(cat)):
            # skip if we've already indexed this event
            if i in index_list:
                dist_list.append(0)
                continue
            lat_ = cat[i].preferred_origin().latitude
            lon_ = cat[i].preferred_origin().longitude
            dist_list.append(locations2degrees(lat1=lat, long1=lon,
                                               lat2=lat_, long2=lon_)
                             )
        # were assuming that the distances are unique, which they should be
        # if we are using floating points
        index_list.append(dist_list.index(max(dist_list)))
        lat = cat[dist_list.index(max(dist_list))].preferred_origin().latitude
        lon = cat[dist_list.index(max(dist_list))].preferred_origin().longitude

    # event id list
    event_list = []
    for i in index_list:
        event_list.append(cat[i].resource_id.id.split('/')[1])

    return index_list, event_list


def generate_focal_mechanism(mtlist, event=None):
    """
    Event objects fetched by obspy do not natively come with any moment tensor
    or nodal plane information, that is stored separately in a .csv file
    located on github. For the tomography problem we need this information, so
    this functino will append information from the .csv file onto the obspy
    event object so that all the information can be located in a single object

    :type mtlist: dict
    :param mtlist; row values from the GeoNet moment tensor csv file
    :type event: obspy.core.event.Event
    :param event: event to append focal mechanism to
    :rtype focal_mechanism: obspy.core.event.FocalMechanism
    :return focal_mechanism: generated focal mechanism
    """
    import obspy.core.event
    import obspy.core.event.source as eventcore
    from pyatoa.utils.operations.conversions import mt_transform

    id_template = "smi:local/geonetcsv/{0}/{1}".format(mtlist['PublicID'],'{}')
    if len(mtlist) != 32:
        print("geonet moment tensor list does not have the correct number"
              "of requisite components, should have 32")
        return

    # NODAL PLANES
    nodal_plane_1 = eventcore.NodalPlane(strike=mtlist['strike1'],
                                         dip=mtlist['dip1'],
                                         rake=mtlist['rake1']
                                         )
    nodal_plane_2 = eventcore.NodalPlane(strike=mtlist['strike2'],
                                         dip=mtlist['dip2'],
                                         rake=mtlist['rake2']
                                         )
    nodal_planes = eventcore.NodalPlanes(nodal_plane_1, nodal_plane_2,
                                         preferred_plane=1
                                         )
    # PRINCIPAL AXES
    tension_axis = eventcore.Axis(azimuth=mtlist['Taz'], plunge=mtlist['Tpl'],
                                  length=mtlist['Tva']
                                  )
    null_axis = eventcore.Axis(azimuth=mtlist['Naz'], plunge=mtlist['Npl'],
                               length=mtlist['Nva']
                               )
    pressure_axis = eventcore.Axis(azimuth=mtlist['Paz'], plunge=mtlist['Ppl'],
                                   length=mtlist['Pva']
                                   )
    principal_axes = eventcore.PrincipalAxes(t_axis=tension_axis,
                                             p_axis=pressure_axis,
                                             n_axis=null_axis)

    # MOMENT TENSOR
    cv = 1E20 * 1E-7  # non-units to dyne*cm to N*m
    seismic_moment_in_nm = mtlist['Mo']*1E-7
    rtp = mt_transform(mt={"m_xx": mtlist['Mxx']*cv, "m_yy": mtlist['Myy']*cv,
                           "m_zz": mtlist['Mzz']*cv, "m_xy": mtlist['Mxy']*cv,
                           "m_xz": mtlist['Mxz']*cv, "m_yz": mtlist['Myz']*cv
                           }, method="xyz2rtp"
                       )
    tensor = eventcore.Tensor(m_rr=rtp['m_rr'], m_tt=rtp['m_tt'],
                              m_pp=rtp['m_pp'], m_rt=rtp['m_rt'],
                              m_rp=rtp['m_rp'], m_tp=rtp['m_tp']
                              )
    source_time_function = eventcore.SourceTimeFunction(
        duration=2 * half_duration_from_m0(seismic_moment_in_nm)
        )
    comment = obspy.core.event.base.Comment(
        force_resource_id=False,
        text="Automatically generated by Pyatoa via GeoNet MT .csv"
    )
    moment_tensor = eventcore.MomentTensor(
        force_resource_id=False, tensor=tensor,
        source_time_function=source_time_function,
        derived_origin_id=id_template.format('origin#ristau'),
        scalar_moment=seismic_moment_in_Nm, double_couple=mtlist['DC']/100,
        variance_reduction=mtlist['VR'], comment=comment
        )

    # FOCAL MECHANISM
    focal_mechanism = eventcore.FocalMechanism(
        force_resource_id=False, nodal_planes=nodal_planes,
        moment_tensor=moment_tensor, principal_axes=principal_axes,
        comments=[comment]
        )

    if event is not None:
        event.focal_mechanisms = [focal_mechanism]

    return focal_mechanism
