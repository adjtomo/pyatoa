"""
Functions used in determining information related to sources and receivers,
or their corresponding representations in ObsPy
"""
import warnings
import numpy as np
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.core.event.source import Tensor


def seismic_moment(mt):
    """
    Return the seismic moment based on a moment tensor.
    Can take a list of tensor components, or a Tensor object from ObsPy.

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


def moment_magnitude(moment, c=10.7):
    """
    Return the moment magitude based on a seismic moment, from
    Hanks & Kanamori (1979)

    :type c: float
    :param c: correction factor for conversion, 10.7 for units of N*m,
        16.1 for units of dyne*cm
    :type moment: float
    :param moment: the seismic moment, in units of N*m
    :rtype: float
    :return: moment magnitude, M_w
    """
    return 2 / 3 * np.log10(moment) - c


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
    Based on Equation from 'Aki and Richards: Quantitative Seismology' book

    .. note::
        Acceptable formats for the parameter `mt`:

        1) [m11, m22, m33, m12, m13, m23]
        2) [mxx, myy, mzz, mxy, mxz, myz]
        3) [mrr, mtt, mpp, mrt, mrp, mtp]


    :type mt: dict
    :param mt: moment tensor in format above
    :type method: str
    :param method: type of conversion, "rtp2xyz" or "xyz2rtp"
    :rtype: dict
    :return: converted moment tensor dictionary
    """
    assert(method in ["xyz2rtp", "rtp2xyz"]), \
        "method must be 'xyz2rtp' or 'rtp2xyz'"

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


def lonlat_utm(lon_or_x, lat_or_y, utm_zone=None, inverse=False):
    """
    Convert latitude and longitude coordinates to UTM projection using PyProj

    :type lon_or_x: float or int
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
    :type lat_or_y: float or int
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
    :type utm_zone: int
    :param utm_zone: UTM zone for conversion from WGS84
    :type inverse: bool
    :param inverse: if inverse == False, latlon => UTM, vice versa.
    :rtype: tuple (float, float)
    :return: (x in UTM or longitude in WGS84, y in UTM or latitude in WGS84)
    """
    from pyproj import Proj

    # If converting latlon to utm and no utm zone given, calculate utm zone
    if utm_zone is None and not inverse:
        utm_zone = utm_zone_from_lat_lon(lat_or_y, lon_or_x)
    elif utm_zone is None and inverse:
        raise TypeError(
            "lonlat_utm() missing 1 required positional argument: 'utm_zone'"
        )
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


def utm_zone_from_lat_lon(lat, lon):
    """
    Calculate the UTM zone longitude value using quick maffs.
    Get the sign of the UTM zone based on the latitude value.

    :type lat: float
    :param lat: latitude coordinate in degrees
    :type lon: float
    :param lon: longitude coordinate in degrees
    :rtype: int
    :return: UTM zone number
    """
    try:
        sign = lat / abs(lat)  # silly way to figure out if lat is +/-
    except ZeroDivisionError as e:
        raise Exception("latitude is 0, UTM zone is ambigious") from e
    return sign * np.ceil((lon + 180) / 6)


def gcd_and_baz(event, sta):
    """
    Calculate great circle distance and backazimuth values for a given
    station and event configuration

    :type event: obspy.core.event.Event
    :param event: event object
    :type sta: obspy.core.inventory.station.Station
    :param sta: station object
    :rtype: tuple (float, float)
    :return: (great circle distance in km, backazimuth in degrees)
    """
    gcdist, _, baz = gps2dist_azimuth(lat1=event.preferred_origin().latitude,
                                      lon1=event.preferred_origin().longitude,
                                      lat2=sta.latitude,
                                      lon2=sta.longitude
                                      )
    return gcdist * 1E-3, baz


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
    :rtype: list
    :return: list of stations in order from 0deg to 360deg in direction
    """
    station_names, list_of_baz = [], []
    event = ds.events[0]
    for sta_name in ds.waveforms.list():
        try:
            sta = ds.waveforms[sta_name].StationXML[0][0]
        except AttributeError:
            warnings.warn(f"station {sta_name} has no attribute StationXML",
                          UserWarning)
            continue
        station_names.append(sta_name)

        _, baz = gcd_and_baz(event, sta)
        list_of_baz.append(baz)

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


class Source:
    """
    A generic Source object to characterize FORCESOLUTION files and SPECFEM2D
    SOURCE files without breaking the architechture of Pyatoa which was built
    around CMTSOLUTIONs and ObsPy Event objects

    Essentially this class tries to mimic the ObsPy Event object and return
    required information that is queried throughout a Pyatoa workflow
    """
    def __init__(self, resource_id, origin_time, longitude, latitude, depth):
        """
        Only define the essential values required of a source

        :type resource_id: str
        :param resource_id: unique label for the event
        :type time: str or UTCDateTime
        :param: origin time for the event
        :type longitude: float
        :param longitude: longitude or X value of the event in the domain
        :type latitude: float
        :param latitude: latitude or Y value of the event in the domain
        :type depth: float
        :param depth: depth in km, inverted Z axis, positive values means deeper
        """
        self.id = resource_id
        self.time = UTCDateTime(origin_time)
        self.longitude = float(longitude)
        self.latitude = float(latitude)
        self.depth = float(depth)

    def preferred_origin(self):
        """
        Convenience function to mimic behavior of ObsPy Event object
        """
        return self

    @property
    def resource_id(self):
        """
        Convenenience function to mimic behavior of ObsPy Event object
        """
        return self







