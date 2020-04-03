"""
Functions used in determining information related to sources and receivers,
or their corresponding representations in ObsPy
"""
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core.event.source import Tensor


def seismic_moment(mt):
    """
    Return the seismic moment based on a moment tensor.
    Can take a list of tensor components, or a Tensor object from ObsPy.

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
        assert (len(mt_temp) == 6), "Moment tensor should have 6 components"
        mt = mt_temp
    return 1 / np.sqrt(2) * np.sqrt(sum([_ ** 2 for _ in mt]))


def moment_magnitude(moment):
    """
    Return the moment magitude based on a seismic moment, from
    Hanks & Kanamori (1979)

    :type moment: float
    :param moment: the seismic moment, in units of N*m
    :rtype: float
    :return: moment magnitude M_w
    """
    return 2 / 3 * np.log10(moment) - 10.7


def half_duration_from_m0(moment):
    """
    Empirical formula for half duration used by Harvard CMT, stated in
    Dahlen and Tromp (1998, p.178).

    :type moment: float
    :param moment: seismic moment in N*m
    :rtype: float
    :return: empirically scaled half duration
    """
    return 2.4E-6 * moment**(1/3)


def mt_transform(mt, method):
    """
    Transform moment tensor between XYZ and RTP coordinates

    Acceptable formats for the parameter mt:
        1) [m11,m22,m33,m12,m13,m23]
        2) [mxx,myy,mzz,mxy,mxz,myz]
        3) [mrr,mtt,mpp,mrt,mrp,mtp]

    Based on equation ?? from Aki and Richards Quantitative Seismology
    TO DO: find the correct equation number

    :type mt: dict
    :param mt: moment tensor in format above
    :type method: str
    :param method: type of conversion, "rtp2xyz" or "xyz2rtp"
    :rtype: dict
    :return: converted moment tensor dictionary
    """
    if method == "xyz2rtp":
        if "m_xx" not in mt.keys():
            print("for xyz2rtp, dict must have keys in xyz")
        m_rr = mt["m_zz"]
        m_tt = mt["m_xx"]
        m_pp = mt["m_yy"]
        m_rt = mt["m_xz"]
        m_rp = -1 * mt["m_yz"]
        m_tp = -1 * mt["m_xy"]
        return {"m_rr": m_rr, "m_tt": m_tt, "m_pp": m_pp, "m_rt": m_rt,
                "m_rp": m_rp, "m_tp": m_tp}
    elif method == "rtp2xyz":
        if "m_tt" not in mt.keys():
            print("for rtp2xyz, dict must have keys in rtp")
        m_xx = mt["m_tt"]
        m_yy = mt["m_pp"]
        m_zz = mt["m_rr"]
        m_xy = -1 * mt["m_tp"]
        m_xz = mt["m_rt"]
        m_yz = -1 * mt["m_rp"]
        return {"m_xx": m_xx, "m_yy": m_yy, "m_zz": m_zz, "m_xy": m_xy,
                "m_xz": m_xz, "m_yz": m_yz}
    else:
        print("Invalid transformation method, xyz2rtp or rtp2xyz")
        return None


def lonlat_utm(lon_or_x, lat_or_y, utm_zone=-60, inverse=False):
    """
    Convert latitude and longitude coordinates to UTM projection

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
    # Determine if the projection is north or south
    if utm_zone < 0:
        direction = "south"
    else:
        direction = "north"
    # Proj doesn't accept negative zones
    utm_zone = abs(utm_zone)

    projstr = (f"+proj=utm +zone={utm_zone}, +{direction} +ellps=WGS84"
               " +datum=WGS84 +units=m +no_defs")
    projection = Proj(projstr)

    x_or_lon, y_or_lat = projection(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat


def gcd_and_baz(event, sta):
    """
    Calculate great circle distance and backazimuth values for a given
    station and event configuration

    :type event: obspy.core.event.Event
    :param event: event object
    :type sta: obspy.core.inventory.station.Station
    :param sta: station object
    :rtype gcdist: float
    :return gcdist: great circle distance in km
    :rtype baz: float
    :return baz: backazimuth in degrees
    """
    gcdist, _, baz = gps2dist_azimuth(lat1=event.preferred_origin().latitude,
                                      lon1=event.preferred_origin().longitude,
                                      lat2=sta.latitude,
                                      lon2=sta.longitude
                                      )
    return gcdist*1E-3, baz


def theoretical_p_arrival(source_depth_in_km, distance_in_degrees,
                          model='iasp91'):
    """
    Calculate theoretical arrivals

    :type source_depth_in_km: float
    :param source_depth_in_km: source depth in units of km
    :type distance_in_degrees: float
    :param distance_in_degrees: distance between source and receiver in degrees
    :type model: str
    :param model: model to be used to get travel times from
    """
    from obspy.taup import TauPyModel
    model = TauPyModel(model=model)
    return model.get_travel_times(source_depth_in_km, distance_in_degrees,
                                  phase_list=["P"])


def merge_inventories(inv_a, inv_b):
    """
    Adding inventories together duplicates network and station codes, which is
    kind of annoying for looping. This function will add two inventories
    together while minimizing the amount of redundant networks, stations,
    channels inside the merged inventory.

    :type inv_a: obspy.core.inventory.Inventory
    :param inv_a: inventory to merge into, will be returned
    :type inv_b: obspy.core.inventory.Inventory
    :param inv_b: inventory to merge into inv_a
    :rtype: obspy.core.inventory.Inventory
    :return: merged inventories
    """
    # Network, loop through 'a' and look for matches in 'b'
    for net_a in inv_a.networks:
        for net_b in inv_b.networks:
            # If the network codes match, need to go deeper
            if net_a.code == net_b.code:
                # Station, Loop through 'a' and look for matches in 'b'
                for sta_a in net_a.stations:
                    for sta_b in net_b.stations:
                        # If station codes match, go deeper
                        if sta_a.code == sta_b.code:
                            # Channel, loop through 'a', look in 'b'
                            for i, cha_a in enumerate(sta_a.channels):
                                for cha_b in sta_b.channels:
                                    # If channel codes match, these
                                    # inventories are probably the same.
                                    # This won't work if they contain the
                                    # same channel but different time
                                    # ranges. But we assume that
                                    # case won't arise.
                                    if cha_a.code == cha_b.code:
                                        continue
                                    # If channels don't match, merge and return
                                    else:
                                        sta_a.channels += sta_b.channels
                                        return inv_a
                        # If station codes don't match, merge and return
                        else:
                            net_a.stations += net_b.stations
                            return inv_a
            # If the network codes don't match, merge and return
            else:
                inv_a.networks += inv_b.networks
                return inv_a
    return inv_a + inv_b


def seismogram_length(distance_km, slow_wavespeed_km_s=2, binsize=50,
                      minimum_length=100):
    """
    Dynamically determine the length of the seismogram based on source-receiver
    distance. Bin into lengths to keep some uniformity in window lengths

    :type distance_km: float
    :param distance_km: source-receiver distance in km
    :type slow_wavespeed_km_s: int
    :param slow_wavespeed_km_s: slowest wavespeed in model, in km/s
    :type binsize: int
    :param binsize: bin size for rounding the length to the nearest value
    :type minimum_length: int
    :param minimum_length: the shortest a seismogram can be
    :rtype: int
    :return: expected seismogram length
    """
    from pyatoa.utils.calculate import myround

    # determine based on slowest wavespeed travelling from source to receiver
    rough_length_s = distance_km / slow_wavespeed_km_s
    if rough_length_s < minimum_length:
        rough_length_s = minimum_length

    return myround(rough_length_s, base=binsize, choice='up')


def sort_by_backazimuth(ds, clockwise=True):
    """
    Its illustrative to show misfit information for an event, sorted by
    backazimuth. Stations with misfit information are generally sorted
    alphabetically, so this function just calcualtes backazimuth and returns a
    sorted list of station names. Can go clockwise or counter clockwise,
    starting from 0 degrees.

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
        To avoid repeat imports and unnecessary returns from gcd_and_baz()
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


def event_by_distance(cat, filter_type=False, filter_bounds=None, random=False):
    """
    Sort through an obspy catalog by interevent distance. If we have a lot of
    events in a catalog, it's best to spatially vary them such that we don't
    redundantly oversample a spatial region.
    Returns an index list for events that are most distant from one another,
    without repeating any used events.

    Catalog filter parameters can be found here:
    https://docs.obspy.org/packages/autogen/obspy.core.event.Catalog.filter.html

    example call:
    >> index_list, event_list = event_by_distance(cat, filter_type="magnitude",
                                                  filter_bounds=[5.0,6.0])

    :type cat: obspy.event.Catalog
    :param cat: catalog to sort through
    :type filter_type: str
    :param filter_type: filter to be passed to the Catalog filter
    :type filter_bounds: list of floats
    :param filter_bounds: (min filter bound, max filter bound)
    :type random: bool
    :param random: randomly determined starting point
    :rtype: obspy.event.Catalog
    :return: filtered catalog object
    """
    from obspy.geodetics import locations2degrees
    from obspy.core.event import Catalog

    # filter the catalog by e.g. magnitude before sorting
    if filter_type:
        cat = cat.filter(f"{filter_type} >= {filter_bounds[0]}")
        cat = cat.filter(f"{filter_type} <= {filter_bounds[1]}")

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
    events = []
    for i in index_list:
        events.append(cat[i])

    # Return a new catalog
    return Catalog(events=events)


def generate_focal_mechanism(mtlist, event=None):
    """
    For the New Zealand Tomography Problem

    Focal mechanisms created by John Ristau are written to a .csv file
    located on Github. This function will append information from the .csv file
    onto the Obspy event object so that all the information can be located in a
    single object

    :type mtlist: dict
    :param mtlist; row values from the GeoNet moment tensor csv file
    :type event: obspy.core.event.Event
    :param event: event to append focal mechanism to
    :rtype focal_mechanism: obspy.core.event.FocalMechanism
    :return focal_mechanism: generated focal mechanism
    """
    from obspy.core.event import source
    from obspy.core.event.base import Comment

    # Match the identifier with Goenet
    id_template = f"smi:local/geonetcsv/{mtlist['PublicID']}/{{}}"

    # Check that the input list is properly formatted
    if len(mtlist) != 32:
        print("geonet moment tensor list does not have the correct number"
              "of requisite components, should have 32")
        return

    # Generate the Nodal Plane objects to append
    nodal_plane_1 = source.NodalPlane(
        strike=mtlist['strike1'], dip=mtlist['dip1'], rake=mtlist['rake1']
    )
    nodal_plane_2 = source.NodalPlane(
        strike=mtlist['strike2'], dip=mtlist['dip2'], rake=mtlist['rake2']
    )
    nodal_planes = source.NodalPlanes(
        nodal_plane_1, nodal_plane_2, preferred_plane=1
    )

    # Create the Principal Axes as Axis objects
    tension_axis = source.Axis(
        azimuth=mtlist['Taz'], plunge=mtlist['Tpl'], length=mtlist['Tva']
    )
    null_axis = source.Axis(
        azimuth=mtlist['Naz'], plunge=mtlist['Npl'], length=mtlist['Nva']
    )
    pressure_axis = source.Axis(
        azimuth=mtlist['Paz'], plunge=mtlist['Ppl'], length=mtlist['Pva']
    )
    principal_axes = source.PrincipalAxes(
        t_axis=tension_axis, p_axis=pressure_axis, n_axis=null_axis
    )

    # Create the Moment Tensor object with correct units and scaling
    cv = 1E20 * 1E-7  # convert non-units, to dyne*cm, to N*m
    seismic_moment_in_nm = mtlist['Mo'] * 1E-7

    # Convert the XYZ coordinate system of GeoNet to an RTP coordinate system
    # expected in the CMTSOLUTION file of Specfem
    rtp = mt_transform(mt={"m_xx": mtlist['Mxx']*cv, "m_yy": mtlist['Myy']*cv,
                           "m_zz": mtlist['Mzz']*cv, "m_xy": mtlist['Mxy']*cv,
                           "m_xz": mtlist['Mxz']*cv, "m_yz": mtlist['Myz']*cv
                           },
                       method="xyz2rtp"
                       )
    tensor = source.Tensor(m_rr=rtp['m_rr'], m_tt=rtp['m_tt'],
                           m_pp=rtp['m_pp'], m_rt=rtp['m_rt'],
                           m_rp=rtp['m_rp'], m_tp=rtp['m_tp']
                           )
    # Create the source time function
    source_time_function = source.SourceTimeFunction(
        duration=2 * half_duration_from_m0(seismic_moment_in_nm)
    )

    # Generate a comment for provenance
    comment = Comment(force_resource_id=False,
                      text="Automatically generated by Pyatoa via GeoNet MT CSV"
                      )

    # Fill the moment tensor object
    moment_tensor = source.MomentTensor(
        force_resource_id=False, tensor=tensor,
        source_time_function=source_time_function,
        derived_origin_id=id_template.format('origin#ristau'),
        scalar_moment=seismic_moment_in_nm, double_couple=mtlist['DC']/100,
        variance_reduction=mtlist['VR'], comment=comment
        )

    # Finally, assemble the Focal Mechanism. Force a resource id so that
    # the event can identify its preferred focal mechanism
    focal_mechanism = source.FocalMechanism(
        force_resource_id=True, nodal_planes=nodal_planes,
        moment_tensor=moment_tensor, principal_axes=principal_axes,
        comments=[comment]
        )

    # Append the focal mechanisms to the event object. Set the preferred
    # focal mechanism so that this attribute can be used in the future
    if event:
        event.focal_mechanisms = [focal_mechanism]
        event.preferred_focal_mechanism_id = focal_mechanism.resource_id
        return event, focal_mechanism
    # If no event is given, just return the focal mechanism
    else:
        return None, focal_mechanism




