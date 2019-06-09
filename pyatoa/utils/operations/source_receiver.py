"""
Functions used in determining source receiver information
"""


def seismogram_length(distance_km, slow_wavespeed_km_s=2, binsize=50,
                      minimum_length=100):
    """
    Dynamically determine the length of the seismogram based on source-receiver
    distance. Bin into lengths to keep some uniformity in window lengths
    :return:
    """
    from pyatoa.utils.operations.calculations import myround

    rough_length_s = distance_km / slow_wavespeed_km_s
    if rough_length_s < minimum_length:
        rough_length_s = minimum_length
    return myround(rough_length_s, base=binsize, choice='up')


def sort_by_backazimuth(ds, clockwise=True):
    """
    Its illustrative to show misfit information for an event, sorted by
    backazimuth. Stations with misfit information are generally sorted
    alphabetically, so this function just calcualtes backazimuth and returns a
    sorted list of station names. Can go clockwise or counter clockwise.
    :type ds: pyasdf.ASDFDataSet()
    :param ds: dataset containing event and station information
    :type clockwise: bool
    :param clockwise: False = counter clockwise
    :rytpe: list
    :return: list of stations in order from 0deg to 360deg in direction
    """
    from obspy.geodetics import gps2dist_azimuth

    def baz(event, sta_stats):
        """
        same as gcd_and_baz below, but avoid repeat import and only returns baz
        """
        _, _, baz = gps2dist_azimuth(
            lat1=event.preferred_origin().latitude,
            lon1=event.preferred_origin().longitude,
            lat2=sta_stats.latitude,
            lon2=sta_stats.longitude
            )
        return baz

    station_names, list_of_baz = [], []
    event = ds.events[0]
    for sta_name in ds.waveforms.list():
        sta = ds.waveforms[sta_name].StationXML[0][0]
        list_of_baz.append(baz(event, sta))
        station_names.append(sta_name)
    list_of_baz, station_names = zip(*sorted(zip(list_of_baz, station_names)))

    if not clockwise:
        station_names.reverse()

    return station_names


def lonlat_utm(lon_or_x, lat_or_y, utm_zone=60, inverse=False):
    """convert latitude and longitude coordinates to UTM projection
    from mesh_gen_helper.py (personal code)

    :type lon_or_x: float or int
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
    :type lat_or_y: float or int
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
    :type utm_zone: int
    :param utm_zone: UTM zone for conversion from WGS84
    :type inverse: bool
    :param inverse: if inverse == False, latlon => UTM, vice versa.
    :rtype x_or_lon: float
    :return x_or_lon: x coordinate in UTM or longitude in WGS84
    :rtype y_or_lat: float
    :return y_or_lat: y coordinate in UTM or latitude in WGS84
    """
    from pyproj import Proj
    projstr = ("+proj=utm +zone={}, +south +ellps=WGS84 +datum=WGS84"
               " +units=m +no_defs".format(utm_zone))
    myProj = Proj(projstr)
    x_or_lon, y_or_lat = myProj(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat


# TO DO: this function was moved to utils.operations.file_generation
# def generate_srcrcv_vtk_file(h5_fid, fid_out, model="m00", utm_zone=60,
#                              event_fid_out=None):
#     """
#     It's useful to visualize source receiver locations in Paraview, alongside
#     sensitivity kernels. VTK files are produced by Specfem, however they are for
#     all receivers, and a source at depth which is sometimes confusing. This
#     function will create source_receiver vtk files using the pyasdf h5 files,
#     with only those receiver that were used in the misfit analysis, and only
#     an epicentral source location, such that the source is visible on a top
#     down view from Paraview.
#     Gives the option to create an event vtk file separate to receivers, for
#     more flexibility in the visualization.
#     :type h5_fid: str
#     :param h5_fid: path to pyasdf h5 file outputted by pyatoa
#     :type fid_out: str
#     :param fid_out: output path and filename to save vtk file e.g. 'test.vtk'
#     :type model: str
#     :param model: h5 is split up by model iteration, e.g. 'm00'
#     :type utm_zone: int
#     :param utm_zone: the utm zone of the mesh, 60 for NZ
#     :type event_fid: str
#     :param event_fid: if event vtk file to be made separately
#     """
# 
#     # lazy import pyasdf
#     import pyasdf
# 
#     # get receiver location information in UTM_60 coordinate system from
#     # pyasdf auxiliary_data. make sure no repeat stations
#     ds = pyasdf.ASDFDataSet(h5_fid)
#     sta_x, sta_y, sta_elv, sta_ids = [], [], [], []
#     if bool(ds.auxiliary_data):
#         for adjsrc in ds.auxiliary_data.AdjointSources[model].list():
#             sta = ds.auxiliary_data.AdjointSources[model][adjsrc]
#             station_id = sta.parameters["station_id"]
#             if station_id in sta_ids:
#                 continue
#             latitude = sta.parameters["latitude"]
#             longitude = sta.parameters["longitude"]
#             elevation_in_m = sta.parameters["elevation_in_m"]
# 
#             x, y = lonlat_utm(lon_or_x=longitude, lat_or_y=latitude,
#                               utm_zone=utm_zone, inverse=False)
#             sta_x.append(x)
#             sta_y.append(y)
#             sta_elv.append(elevation_in_m)
#             sta_ids.append(station_id)
# 
#     # get event location information in UTM_60. set depth at 100 for epicenter
#     ev_x, ev_y = lonlat_utm(lon_or_x=ds.events[0].preferred_origin().longitude,
#                             lat_or_y=ds.events[0].preferred_origin().latitude,
#                             utm_zone=utm_zone, inverse=False
#                             )
#     ev_elv = 100.0
# 
#     # write header for vtk file and then print values for source receivers
#     with open(fid_out, "w") as f:
#         f.write("# vtk DataFile Version 2.0\n"
#                 "Source and Receiver VTK file from Pyatoa\n"
#                 "ASCII\n"
#                 "DATASET POLYDATA\n"
#                 )
#         # num points equal to number of stations plus 1 event
#         f.write("POINTS\t{} float\n".format(len(sta_x)+1))
#         f.write("{X:18.6E}{Y:18.6E}{E:18.6E}\n".format(
#             X=ev_x, Y=ev_y, E=ev_elv)
#         )
#         for x, y, e in zip(sta_x, sta_y, sta_elv):
#             f.write("{X:18.6E}{Y:18.6E}{E:18.6E}\n".format(X=x, Y=y, E=e))
# 
#     # make a separate vtk file for the source
#     if event_fid_out:
#         with open(event_fid_out, "w") as f:
#             f.write("# vtk DataFile Version 2.0\n"
#                     "Source and Receiver VTK file from Pyatoa\n"
#                     "ASCII\n"
#                     "DATASET POLYDATA\n"
#                     )
#             f.write("POINTS\t1 float\n".format(len(sta_x) + 1))
#             f.write("{X:18.6E}{Y:18.6E}{E:18.6E}\n".format(
#                 X=ev_x, Y=ev_y, E=ev_elv)
#             )


def gcd_and_baz(event, sta):
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
                                      lat2=sta.latitude,
                                      lon2=sta.longitude
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
        scalar_moment=seismic_moment_in_nm, double_couple=mtlist['DC']/100,
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
