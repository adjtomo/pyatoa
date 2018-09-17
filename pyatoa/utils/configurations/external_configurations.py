"""
Pyflex and Pyadjoint require configuration objects to run; functions here will
set them properly and return the necessary parameters so that these
packages can run properly
"""


def set_pyflex_configuration(config, inv, event):
    """
    TODO check if pyflex config has been set and don't do it again?
    :type config: pyatoa.core.config.Config
    :param config: config object which sets the values of pyflex config
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory containing station of interest
    :type event: obspy.core.event.Event
    :param event: event information for relevant earthquake
    :rtype pf_config: pyflex.Config
    :return pf_config: a pyflex configuration object which helps run pyflex
    :rtype pf_station: pyflex.Station
    :return pf_station: pyflex station object specificying location
    :rtype pf_event pyflex.Event
    :return pf_event: pyflex event object specifyin event location and time
    """
    import pyflex

    pf_config = pyflex.Config(
        min_period=config.min_period, max_period=config.max_period,
        stalta_waterlevel=config.pyflex_config[0],
        tshift_acceptance_level=config.pyflex_config[1],
        dlna_acceptance_level=config.pyflex_config[2],
        cc_acceptance_level=config.pyflex_config[3],
        c_0=config.pyflex_config[4], c_1=config.pyflex_config[5],
        c_2=config.pyflex_config[6], c_3a=config.pyflex_config[7],
        c_3b=config.pyflex_config[8], c_4a=config.pyflex_config[9],
        c_4b=config.pyflex_config[10]
        )
    pf_station = pyflex.Station(latitude=inv[0][0].latitude,
                                longitude=inv[0][0].longitude
                                )
    pf_event = pyflex.Event(latitude=event.origins[0].latitude,
                            longitude=event.origins[0].longitude,
                            depth_in_m=event.origins[0].depth,
                            origin_time=event.origins[0].time
                            )

    return pf_config, pf_event, pf_station


def set_pyadjoint_configuration(config):
    """
    :type config: pyatoa.core.config.Config
    :param config: config object which sets the pyadjoint type
    :rtype cfgout: pyadjoint.Config*
    :return cfgout: properly set pyadjoint configuration object
    """
    import pyadjoint

    if config.adj_src_type == "waveform":
        cfgout = pyadjoint.ConfigWaveForm(min_period=config.min_period,
                                          max_period=config.max_period,
                                          taper_type="hann",
                                          taper_percentage=0.15)
    elif config.adj_src_type == "cc_traveltime_misfit":
        cfgout = pyadjoint.ConfigCrossCorrelation(
            min_period=config.min_period, max_period=config.max_period,
            taper_type='hann', taper_percentage=0.3, measure_type='dt',
            use_cc_error=True, dt_sigma_min=1.0, dlna_sigma_min=0.5)
    elif config.adj_src_type == "multitaper_misfit":
        cfgout = pyadjoint.ConfigMultiTaper(
            min_period=config.min_period, max_period=config.max_period,
            lnpt=15, transfunc_waterlevel=1e-10, water_threshold=0.02,
            ipower_costaper=10, min_cycle_in_window=0.5, taper_type='hann',
            taper_percentage=0.3, mt_nw=4.0, num_taper=5, dt_fac=2.0,
            phase_step=1.5, err_fac=2.5, dt_max_scale=3.5, measure_type='dt',
            dt_sigma_min=1.0, dlna_sigma_min=0.5, use_cc_error=True,
            use_mt_error=False)
    else:
        raise KeyError("adjoint source type incorrectly specified")
    return cfgout

