"""
Pyflex requires configuration objects to run;
the functions here will set them properly and return the necessary parameters

For variable descriptions see:
    https://krischer.github.io/pyflex/_modules/pyflex/config.html
"""
from pyflex import Config as pyflexConfig


def set_pyflex_config(min_period, max_period, choice=None, **kwargs):
    """
    Overwriting the default Pyflex parameters with User-defined criteria

    To Do: add the config options from Maggi et al. 2009 for Japan, SoCal, Globe

    :type choice: str
    :param choice: name of map to choose the Pyflex config options, if None,
        default values are used. Kwargs can still overload default values.
    :type min_period: float
    :param min_period: min period of the data
    :type max_period: float
    :param max_period: max period of the data
    :rtype: pyflex.Config
    :return: the pyflex Config option to use when running Pyflex
    """
    # Instantiate the pyflex Config object
    pfconfig = pyflexConfig(min_period=min_period, max_period=max_period)

    # From the example on the Pyflex website
    if choice == "example":
        setattr(pfconfig, "stalta_waterlevel", 0.08)
        setattr(pfconfig, "tshift_acceptance_level", 15.0)
        setattr(pfconfig, "dlna_acceptance_level", 1.0)
        setattr(pfconfig, "cc_acceptance_level", 0.8)
        setattr(pfconfig, "c_0", 0.7)
        setattr(pfconfig, "c_1", 4.0)
        setattr(pfconfig, "c_3a", 1.0)
        setattr(pfconfig, "c_3b", 2.0)
        setattr(pfconfig, "c_4a", 3.0)
        setattr(pfconfig, "c_4b", 10.0)
    # From the UAF group doing regional studies of Alaska
    elif choice == "alaska":
        setattr(pfconfig, "stalta_waterlevel", 0.18)
        setattr(pfconfig, "tshift_acceptance_level", 4.0)
        setattr(pfconfig, "dlna_acceptance_level", 1.5)
        setattr(pfconfig, "cc_acceptance_level", 0.71)
        setattr(pfconfig, "c_0", 0.7)
        setattr(pfconfig, "c_1", 2.0)
        setattr(pfconfig, "c_3a", 3.0)
        setattr(pfconfig, "c_3b", 2.0)
        setattr(pfconfig, "c_4a", 2.5)
        setattr(pfconfig, "c_4b", 12.0)
    # From the VUW/GNS group doing regional studies of North Island, New Zealand
    elif choice == "hikurangi":
        setattr(pfconfig, "stalta_waterlevel", 0.18)
        setattr(pfconfig, "tshift_acceptance_level", 8.0)  # based on sign flip
        setattr(pfconfig, "dlna_acceptance_level", 1.5)
        setattr(pfconfig, "cc_acceptance_level", 0.7)
        setattr(pfconfig, "s2n_limit", 3.)
        setattr(pfconfig, "max_time_before_first_arrival", 0.)  # min strt wind.
        setattr(pfconfig, "c_0", 0.7)
        setattr(pfconfig, "c_1", 2.5)  # min win len = c_1 * min_period
        setattr(pfconfig, "c_3a", 3.0)
        setattr(pfconfig, "c_3b", 2.0)
        setattr(pfconfig, "c_4a", 2.5)
        setattr(pfconfig, "c_4b", 12.0)
    elif choice == "CUSTOM_CHOICE_1":
        # SET ATTRIBUTES FOR CUSTOM PYFLEX CONFIGS HERE
        pass

    # Kwargs can also be passed from the pyatoa.Config object to avoid having to
    # define pre-set values. Kwargs will override preset values
    for key, item in kwargs.items():
        if hasattr(pfconfig, key):
            setattr(pfconfig, key, item)

    return pfconfig


def set_pyflex_station_event(inv, event):
    """
    DEPRECATED
    This isn't necessary anymore because Pyflex can take Inventory and Event
    objects so we skip this level of abstraction

    Set event and station objects expected by Pyflex.

    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing station of interest
    :type event: obspy.core.event.Event
    :param event: event information for relevant earthquake
    :rtype pf_station: pyflex.Station
    :return pf_station: pyflex station object specificying location
    :rtype pf_event pyflex.Event
    :return pf_event: pyflex event object specifyin event location and time
    """
    from pyflex import Station, Event
    pf_station = Station(latitude=inv[0][0].latitude,
                         longitude=inv[0][0].longitude
                         )
    pf_event = Event(latitude=event.preferred_origin().latitude,
                     longitude=event.preferred_origin().longitude,
                     depth_in_m=event.preferred_origin().depth,
                     origin_time=event.preferred_origin().time
                     )

    return pf_station, pf_event

