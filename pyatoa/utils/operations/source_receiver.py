"""functions used in determining source receiver information
"""


def gcd_and_baz(event, inv):
    """calculate great circle distance and backazimuth values for a given
    station and event configuration
    """
    from obspy.geodetics import gps2dist_azimuth
    gcdist, _, baz = gps2dist_azimuth(lat1=event.origins[0].latitude,
                                      lon1=event.origins[0].longitude,
                                      lat2=inv[0][0].latitude,
                                      lon2=inv[0][0].longitude
                                      )
    return gcdist*1E-3, baz


def theoretical_p_arrival(source_depth_in_km, distance_in_degrees,
                          model='iasp91'):
    """calculate theoreitcal arrivals
    """
    from obspy.taup import TauPyModel
    model = TauPyModel(model=model)
    return model.get_travel_times(source_depth_in_km, distance_in_degrees,
                                  phase_list=["P"])


def source_receiver_dictionary(event, moment_tensor, inv=None):
    """determine source receiver parameters such as great circle distance,
    backazimuth etc., plot the source and receiver with a line. have the ability
    to not include an inventory, which will just plot the event then
    source = A, receiver = B
    """
    event_lat, event_lon = (event.origins[0].latitude,
                            event.origins[0].longitude)
    depth = event.origins[0].depth*1E-3
    origintime = event.origins[0].time
    origintime.precision = 0
    for magni in event.magnitudes:
        if magni.magnitude_type == "M":
            magnitude = magni.mag
            magnitude_type = magni.magnitude_type

    if inv:
        station = inv[0][0].code
        sta_lat, sta_lon = (inv[0][0].latitude, inv[0][0].longitude)
        gcdist, baz = gcd_and_baz(inv, event)
    else:
        station, sta_lat, sta_lon = None, None, None
        gcdist, az, baz = None, None, None

    # dictionary output for use in annotations
    srcrcvdict = {"station": station, "sta_lat": sta_lat, "sta_lon": sta_lon,
                  "ev_lat": event_lat,
                  "ev_lon": event_lon,
                  "event_id": event_id,
                  "distance": gcdist,
                  "backazimuth": baz,
                  "date": origintime,
                  "depth": depth,
                  "magnitude":magnitude,
                  "magnitude_type":magnitude_type
                  }

    # connect source receiever with line and color receiver, plot event
    event_x,event_y = m(event_lon,event_lat)
    if inv:
        sta_x,sta_y = m(sta_lon,sta_lat)
        X,Y = [event_x,sta_x],[event_y,sta_y]
        m.plot(X,Y,'--',linewidth=1.1,c='k',zorder=2)
    else:
        X,Y = event_x,event_y
    m.scatter(X,Y,marker='v',color='r',edgecolor='k',s=80,zorder=6)

    beachballcheck = event_beachball(m,MT)
    if not beachballcheck:
        m.scatter(event_x,event_y,marker='o',color='r',edgecolor='k',
                                                        s=200,zorder=6)

    return srcrcvdict