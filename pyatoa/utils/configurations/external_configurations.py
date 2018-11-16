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

    parameters (from Pyflex docs)

    min_period (float)
        Minimum period of the filtered input data in seconds.
    max_period (float)
        Maximum period of the filtered input data in seconds.
    stalta_waterlevel (float or numpy.ndarray)
        Water level on the STA/LTA functional.
    tshift_acceptance_level (float or numpy.ndarray)
        Maximum allowable cross correlation time shift/lag relative to the reference.
    tshift_reference (float)
        Allows a systematic shift of the cross correlation time lag.
    dlna_acceptance_level (float or numpy.ndarray)
        Maximum allowable amplitude ratio relative to the reference.
    dlna_reference (float)
        Reference DLNA level. Allows a systematic shift.
    cc_acceptance_level (float or numpy.ndarray)
        Limit on the normalized cross correlation per window.
    s2n_limit (float or numpy.ndarray)
        Limit of the signal to noise ratio per window.
        If the maximum amplitude of the window over the maximum amplitude of the global noise of the waveforms is smaller than this window,
        then it will be rejected.
    earth_model (str)
        The earth model used for the traveltime calculations. Either "ak135" or "iasp91".
    min_surface_wave_velocity (float)
        The minimum surface wave velocity in km/s. All windows containing data later then this velocity will be rejected.
        Only used if station and event information is available.
    max_time_before_first_arrival (float)
        This is the minimum starttime of any window in seconds before the first arrival.
        No windows will have a starttime smaller than this.
    c_0 (float)
        Fine tuning constant for the rejection of windows based on the height of internal minima.
        Any windows with internal minima lower then this value times the STA/LTA water level at the window peak will be rejected.
    c_1 (float)
        Fine tuning constant for the minimum acceptable window length.
        This value multiplied by the minimum period will be the minimum acceptable window length.
    c_2 (float)
        Fine tuning constant for the maxima prominence rejection.
        Any windows whose minima surrounding the central peak are smaller then this value times the central peak will be rejected.
        This value is set to 0 in many cases as it is hard to control.
    c_3a (float)
        Fine tuning constant for the separation height in the phase separation rejection stage.
    c_3b (float)
        Fine tuning constant for the separation time used in the decay function in the phase separation rejection stage.
    c_4a (float)
        Fine tuning constant for curtailing windows on the left with emergent start/stops and/or codas.
    c_4b (float)
        Fine tuning constant for curtailing windows on the right with emergent start/stops and/or codas.
    check_global_data_quality
        Determines whether or not to check the signal to noise ratio of the whole observed waveform.
        If True, no windows will be selected if the signal to noise ratio is above the thresholds.
    snr_integrate_base (float)
        Minimal SNR ratio. If the squared sum of the signal normalized by its length over the squared sum of the
        noise normalized by its length is smaller then this value, no windows will be chosen for the waveforms.
        Only used if check_global_data_quality is True.
    snr_max_base (float)
        Minimal amplitude SNR ratio. If the maximum amplitude of the signal over the maximum amplitude of the noise
        is smaller than this value no windows will be chosen for the waveforms.
        Only used if check_global_data_quality is True.
    noise_start_index (int)
        Index in the observed data where noise starts for the signal to noise calculations.
    noise_end_index (int)
        Index in the observed data where noise ends for the signal to noise calculations.
        Will be set to the time of the first theoretical arrival minus the minimum period if not set and event and station information is available.
    signal_start_index (int)
        Index where the signal starts for the signal to noise calculations.
        Will be set to to the noise end index if not given.
    signal_end_index (int)
        Index where the signal ends for the signal to noise calculations.
    window_weight_fct (function)
        A function returning the weight for a specific window as a single number.
        Directly passed to the Window ‘s initialization function.
    window_signal_to_noise_type (str)
        The type of signal to noise ratio used to reject windows. If "amplitude",
        then the largest amplitude before the arrival is the noise amplitude and the largest amplitude in the window is the signal amplitude.
        If "energy" the time normalized energy is used in both cases. The later one is a bit more stable when having random wiggles before the first arrival.
    resolution_strategy (str)
        Strategy used to resolve overlaps. Possibilities are "interval_scheduling" and "merge".
        Interval scheduling will chose the optimal subset of non-overlapping windows. Merging will simply merge overlapping windows.

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
    pf_event = pyflex.Event(latitude=event.preferred_origin().latitude,
                            longitude=event.preferred_origin().longitude,
                            depth_in_m=event.preferred_origin().depth,
                            origin_time=event.preferred_origin().time
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

