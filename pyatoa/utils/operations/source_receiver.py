"""functions used in determining source receiver information
"""

def gcd_and_baz(inventory,event):
    """calculate great circle distance and backazimuth values for a given
    station and event configuration
    """
    from obspy.geodetics import gps2dist_azimuth
    GCDist,Az,BAz = gps2dist_azimuth(lat1=event.origins[0].latitude,
        lon1=event.origins[0].longitude,lat2=inv[0][0].latitude,
        lon2=inv[0][0].longitude)
    return GCDist*1E-3, BAz

def theoretical_p_arrival(source_depth_in_km,distance_in_degrees,model='iasp91'):
    """calculate theoreitcal arrivals
    """
    from obspy.taup import TauPyModel
    model = TauPyModel(model=model)
    return model.get_travel_times(source_depth_in_km,distance_in_degrees,
        phase_list=["P"])
