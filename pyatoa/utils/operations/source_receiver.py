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

