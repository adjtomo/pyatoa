"""
Pyflex requires configuration objects to run;
the functions here will set them properly and return the necessary parameters

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
import warnings
from pyflex import Config, Station, Event


def set_pyflex_config(choice, min_period, max_period):
    """
    Overwriting the default Pyflex parameters with User-defined criteria

    To Do: add the config options from Maggi et al. 2009 for Japan, SoCal, Globe

    :type choice: str
    :param choice: name of map to choose the Pyflex config options
    :type min_period: float
    :param min_period: min period of the data
    :type max_period: float
    :param max_period: max period of the data
    :rtype: pyflex.Config
    :return: the pyflex Config option to use when running Pyflex
    """
    available_cfgs = {
        # Default config parameters, i.e. don't overwrite anything
        "default": {},
        # Config parameter values from the Pyflex example script
        "example": {
            "stalta_waterlevel": 0.08,
            "tshift_acceptance_level": 15.0,
            "dlna_acceptance_level": 1.0,
            "cc_acceptance_level": 0.8,
            "c_0": 0.7,
            "c_1": 4.0,
            "c_3a": 1.0,
            "c_3b": 2.0,
            "c_4a": 3.0,
            "c_4b": 10.0
        },
        # Config parameter values from researchers at UAF studying Alaska
        "alaska": {
            "stalta_waterlevel": 0.18,
            "tshift_acceptance_level": 4.0,
            "dlna_acceptance_level": 1.5,
            "cc_acceptance_level": 0.71,
            "c_0": 0.7,
            "c_1": 2.0,
            "c_3a": 3.0,
            "c_3b": 2.0,
            "c_4a": 2.5,
            "c_4b": 12.0
        },
        # Config values from researchers at VUW/GNS studying New Zealand NI
        "hikurangi": {
            "stalta_waterlevel": 0.18,
            "tshift_acceptance_level": 8.0,  # based on sign flip waveforms
            "dlna_acceptance_level": 1.5,
            "cc_acceptance_level": 0.7,
            "s2n_limit": 3.,
            "max_time_before_first_arrival": 0.,  # min start of window (s)
            "c_0": 0.7, 
            "c_1": 2.5,  # min win length = c_1 * min_period
            "c_3a": 3,
            "c_3b": 2.,
            "c_4a": 2.5,
            "c_4b": 12.,
            },
        # User can choose to overwrite parameters with a custom map
        "custom": {},
    }

    # instantiate the pyflex Config object
    pf_config = Config(min_period=min_period, max_period=max_period)

    # Check if given key is available, if not set to default
    if choice not in available_cfgs.keys():
        warnings.warn(
            "Pyflex Config choice not in available keys: {}"
            "setting 'default' parameter map.".format(available_cfgs.keys()),
            UserWarning
        )

    # Overwrite the default config values 
    overwrite = available_cfgs[choice]
    for key in overwrite.keys():
        setattr(pf_config, key, overwrite[key])
    
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
    pf_station = Station(latitude=inv[0][0].latitude,
                         longitude=inv[0][0].longitude
                         )
    pf_event = Event(latitude=event.preferred_origin().latitude,
                     longitude=event.preferred_origin().longitude,
                     depth_in_m=event.preferred_origin().depth,
                     origin_time=event.preferred_origin().time
                     )

    return pf_station, pf_event

