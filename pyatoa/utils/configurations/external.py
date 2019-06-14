"""
Pyflex and Pyadjoint require configuration objects to run; 
the functions here will set them properly and return the necessary parameters 
so that these packages can run properly

   Pyflex Parameter Descrpitions (from Pyflex docs):

    min_period (float)
        Minimum period of the filtered input data in seconds.
    max_period (float)
        Maximum period of the filtered input data in seconds.
    stalta_waterlevel (float or numpy.ndarray)
        Water level on the STA/LTA functional.
    tshift_acceptance_level (float or numpy.ndarray)
        Maximum allowable cross correlation time shift/lag relative to the
        reference.
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
        If the maximum amplitude of the window over the maximum amplitude of the
        global noise of the waveforms is smaller than this window,
        then it will be rejected.
    earth_model (str)
        The earth model used for the traveltime calculations.
        Either "ak135" or "iasp91".
    min_surface_wave_velocity (float)
        The minimum surface wave velocity in km/s. All windows containing data
        later then this velocity will be rejected.
        Only used if station and event information is available.
    max_time_before_first_arrival (float)
        This is the minimum starttime of any window in seconds before the first
        arrival. No windows will have a starttime smaller than this.
    c_0 (float)
        Fine tuning constant for the rejection of windows based on the height of
        internal minima. Any windows with internal minima lower then this value
        times the STA/LTA water level at the window peak will be rejected.
    c_1 (float)
        Fine tuning constant for the minimum acceptable window length.
        This value multiplied by the minimum period will be the minimum
        acceptable window length.
    c_2 (float)
        Fine tuning constant for the maxima prominence rejection.
        Any windows whose minima surrounding the central peak are smaller then
        this value times the central peak will be rejected.
        This value is set to 0 in many cases as it is hard to control.
    c_3a (float)
        Fine tuning constant for the separation height in the phase separation
        rejection stage.
    c_3b (float)
        Fine tuning constant for the separation time used in the decay function
        in the phase separation rejection stage.
    c_4a (float)
        Fine tuning constant for curtailing windows on the left with emergent
        start/stops and/or codas.
    c_4b (float)
        Fine tuning constant for curtailing windows on the right with emergent
        start/stops and/or codas.
    check_global_data_quality
        Determines whether or not to check the signal to noise ratio of the
        whole observed waveform. If True, no windows will be selected if the
        signal to noise ratio is above the thresholds.
    snr_integrate_base (float)
        Minimal SNR ratio. If the squared sum of the signal normalized by its
        length over the squared sum of the noise normalized by its length is
        smaller then this value, no windows will be chosen for the waveforms.
        Only used if check_global_data_quality is True.
    snr_max_base (float)
        Minimal amplitude SNR ratio. If the maximum amplitude of the signal over
        the maximum amplitude of the noise is smaller than this value no windows
        will be chosen for the waveforms.
        Only used if check_global_data_quality is True.
    noise_start_index (int)
        Index in the observed data where noise starts for the signal to noise
        calculations.
    noise_end_index (int)
        Index in the observed data where noise ends for the signal to noise
        calculations. Will be set to the time of the first theoretical arrival
        minus the minimum period if not set and event and station information is
         available.
    signal_start_index (int)
        Index where the signal starts for the signal to noise calculations.
        Will be set to to the noise end index if not given.
    signal_end_index (int)
        Index where the signal ends for the signal to noise calculations.
    window_weight_fct (function)
        A function returning the weight for a specific window as a single number
        Directly passed to the Window's initialization function.
    window_signal_to_noise_type (str)
        The type of signal to noise ratio used to reject windows. If "amplitude"
        then the largest amplitude before the arrival is the noise amplitude and
        the largest amplitude in the window is the signal amplitude.
        If "energy" the time normalized energy is used in both cases. The later
        one is a bit more stable when having random wiggles before the first
        arrival.
    resolution_strategy (str)
        Strategy used to resolve overlaps.
        Possibilities are "interval_scheduling" and "merge".
        Interval scheduling will chose the optimal subset of non-overlapping
        windows. Merging will simply merge overlapping windows.
"""


def pyflex_configs():
    """
    Manual set of pyflex overwrites
    :return:
    """
    configs = {
        # default, don't overwrite the pyflex parameters
        "default": {},
        # values taken from the Pyflex example script
        "example": {
            "stalta_waterlevel": 0.08,
            "tshift_acceptance_level": 15.0,
            "dlna_acceptance_level": 1.0,
            "cc_acceptance_level": 0.8,
            "c_0": 0.7,
            "c_1": 4.0,
            "c_2": 0.0,
            "c_3a": 1.0,
            "c_3b": 2.0,
            "c_4a": 3.0,
            "c_4b": 10.0
        },
        # values from Carl Tape at University of Alaska, Fairbanks
        "fairbanks": {
            "stalta_waterlevel": 0.18,
            "tshift_acceptance_level": 4.0,
            "dlna_acceptance_level": 1.5,
            "cc_acceptance_level": 0.71,
            "c_0": 0.7,
            "c_1": 2.0,
            "c_2": 0.0,
            "c_3a": 3.0,
            "c_3b": 2.0,
            "c_4a": 2.5,
            "c_4b": 12.0
        },
        # values tested for the Hikurangi tomography problem, strict criteria
        "hikurangi_strict": {
            "stalta_waterlevel": 0.08,
            "tshift_acceptance_level": 10.0,
            "dlna_acceptance_level": 0.75,
            "cc_acceptance_level": 0.9,
        },
    }

    return configs


def set_pyflex_config(config):
    """
    Set the Pyflex configuration based on Pyatoa Config parameter
    
    TO DO:
        1) check if pyflex config has been set and don't do it again?
        2) set Pyatoa config either as a keyword, or as a dictionary object?

    :type config: pyatoa.core.config.Config
    :param config: config object which sets the values of pyflex config
    :rtype pf_config: pyflex.Config
    :return pf_config: a pyflex configuration object which helps run pyflex

    """
    from pyflex import Config

    # Pyflex defaults, so that Pyatoa can overwrite them manually
    # when setting the config object
    # !!! DO NOT CHANGE THESE VALUES, THEY REFLECT THE DEFAULT PYFLEX VALUES !!!
    # !!! instead overwrite them using _pyflex_overwrites() !!!
    cfg = {
        "min_period": config.min_period,    # min period of filtered data (s)
        "max_period": config.max_period,    # max period of filtered data (s)
        "stalta_waterlevel": 0.07,          # water level of STA/LTA
        "tshift_acceptance_level": 10.0,    # maximum cc time shift wrt ref.
        "tshift_reference": 0.0,            # systematic shift of cc time lag
        "dlna_acceptance_level": 1.3,       # max amplitude ratio wrt ref.
        "dlna_reference": 0.0,              # systematic shift of DLNA
        "cc_acceptance_level": 0.7,         # limit on normalized cc per window
        "s2n_limit": 1.5,                   # limit on SNR per window
        "earth_model": "ak135",             # earth model for traveltime calcs
        "min_surface_wave_velocity": 3.0,   # min vel (km/s), for srcrcv calcs
        "max_time_before_first_arrival": 50.0,  # min starttime of window (s)
        "c_0": 1.0,                         # height of internal minima
        "c_1": 1.5,                         # min window length (* min_period)
        "c_2": 0.0,                         # max prominence rejection
        "c_3a": 4.0,                        # phase separatin height
        "c_3b": 2.5,                        # separation time in decay function
        "c_4a": 2.0,                        # emergent start/stops (left)
        "c_4b": 5.0,                        # emergent start/stops (right)
        "check_global_data_quality": False,  # check SNR of observed
        "snr_integrate_base": 3.5,          # minimal SNR ratio
        "snr_max_base": 3.0,                # min amplitude SNR ratio
        "noise_start_index": 0,             # index where noise starts for SNR
        "noise_end_index": None,            # index where noise ends
        "signal_start_index": None,         # index where signal starts
        "signal_end_index": -1,             # index where signal ends
        "window_weight_fct": None,          # Weighting for specific window
        "window_signal_to_noise_type": "amplitude",     # Type of SNR
        "resolution_strategy": "interval_scheduling"    # overlap resolution
    }

    # overwrite default cfg values
    overwrite = pyflex_configs()[config.pyflex_config[0]]
    for key in overwrite.keys():
        cfg[key] = overwrite[key]

    # set pyflex config
    pf_config = Config(
        min_period=config.min_period,
        max_period=config.max_period,
        stalta_waterlevel=cfg["stalta_waterlevel"],
        tshift_acceptance_level=cfg["tshift_acceptance_level"],
        tshift_reference=cfg["tshift_reference"],
        dlna_acceptance_level=cfg["dlna_acceptance_level"],
        dlna_reference=cfg["dlna_reference"],
        cc_acceptance_level=cfg["cc_acceptance_level"],
        s2n_limit=cfg["s2n_limit"],
        earth_model=cfg["earth_model"],
        min_surface_wave_velocity=cfg["min_surface_wave_velocity"],
        max_time_before_first_arrival=cfg["max_time_before_first_arrival"],
        c_0=cfg["c_0"],
        c_1=cfg["c_1"],
        c_2=cfg["c_2"],
        c_3a=cfg["c_3a"],
        c_3b=cfg["c_3b"],
        c_4a=cfg["c_4a"],
        c_4b=cfg["c_4b"],
        check_global_data_quality=cfg["check_global_data_quality"],
        snr_integrate_base=cfg["snr_integrate_base"],
        snr_max_base=cfg["snr_max_base"],
        noise_start_index=cfg["noise_start_index"],
        signal_start_index=cfg["signal_start_index"],
        signal_end_index=cfg["signal_end_index"],
        window_weight_fct=cfg["window_weight_fct"],
        window_signal_to_noise_type=cfg["window_signal_to_noise_type"],
        resolution_strategy=cfg["resolution_strategy"]
        )

    return pf_config


def set_pyflex_station_event(inv, event):
    """
    Set event and station objects expected by Pyflex

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


def set_pyadjoint_config(config):
    """
    Set the Pyadjoint config based on Pyatoa config parameter    

    :type config: pyatoa.core.config.Config
    :param config: config object which sets the pyadjoint type
    :rtype cfgout: pyadjoint.Config*
    :return cfgout: properly set pyadjoint configuration object
    """
    import pyadjoint

    if config.pyadjoint_config[0] == "waveform":
        cfgout = pyadjoint.ConfigWaveForm(min_period=config.min_period,
                                          max_period=config.max_period,
                                          taper_type="hann",
                                          taper_percentage=0.15)
    elif config.pyadjoint_config[0] == "cc_traveltime_misfit":
        cfgout = pyadjoint.ConfigCrossCorrelation(
            min_period=config.min_period, max_period=config.max_period,
            taper_type='hann', taper_percentage=0.3, measure_type='dt',
            use_cc_error=True, dt_sigma_min=1.0, dlna_sigma_min=0.5)
    elif config.pyadjoint_config[0] == "multitaper_misfit":
        cfgout = pyadjoint.ConfigMultiTaper(
            min_period=config.min_period, max_period=config.max_period,
            lnpt=15, transfunc_waterlevel=1e-10, water_threshold=0.02,
            ipower_costaper=10, min_cycle_in_window=0.5, taper_type='hann',
            taper_percentage=0.3, mt_nw=4.0, num_taper=5, dt_fac=2.0,
            phase_step=1.5, err_fac=2.5, dt_max_scale=3.5, measure_type='dt',
            dt_sigma_min=1.0, dlna_sigma_min=0.5, use_cc_error=True,
            use_mt_error=False)
    else:
        raise KeyError("adjoint source type incorrectly specified, "
                       "must be: 'waveform', 'cc_traveltime_misfit', "
                       "'multitaper_misfit'")
    return cfgout

