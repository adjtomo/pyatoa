"""
Pyflex requires configuration objects to set the number of available variables
that tune the Flexwin algorithm.

This file contains some preset maps for Pyflex that are mostly taken directly
from the Flexwin publication Maggi et al. 2009.

For variable descriptions see:
    https://krischer.github.io/pyflex/_modules/pyflex/config.html

+ Descriptions of a few commonly used parameters that are not self explanatory

    1. Short Term Average Long Term Average water level
        :stalta_waterlevel (float): reject windows where sta/lta waveform dips
            below this threshold value. between 0 and 1
    2. Water level rejection
        :c_0: reject if window.stalta.min() < c_0 * stalta_waterlevel
        :c_1: min_acceptable_window_length = c_1 * T_min
    3. Prominence rejection
        :c_2: reject if window.stalta.min() < c_2 * window.stalta.max()
    4. Separation height in phase separation
        :c_3a: d_stlta > c_3a * d_stalta_center * f_time
            where d_stalta = current max height above min
            and   d_stalta_center = central max height above min
            and   f_time = time decay function
    5. Emergent start/stops and coda wave curtailing
        :c_4a: time_decay_left = T_min * c_4a / dt
        :c_4b: time_decay_right: T_min * c_4b / dt
"""

pyflex_presets = {
    # Empty preset or 'default' will use default init values
    "default": {
    },
    # Example configuration from the Pyflex website
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
    # For use with testing the workflow using a homogeneous halfspace example
    "homogeneous_halfspace": {
        "stalta_waterlevel": 0.05, 
        "tshift_acceptance_level": 15.0,
        "dlna_acceptance_level": 2.0,
        "s2n_limit": 3.,
        "min_surface_wave_velocity": 1.5,
        "c_0": 0.7,
        "c_1": 2., 
        "c_3a": 3.0,
        "c_3b": 2.0,
        "c_4a": 2.5,
        "c_4b": 12.0
    },
    # NZNORTH: From the New Zealand group doing regional North Island studies
    # These are the main parameters used in Chow et al. (2020)
    "nznorth_10-30s_chow_et_al": {
        "stalta_waterlevel": 0.10,
        "tshift_acceptance_level": 8.0,  # based on sign-flip
        "dlna_acceptance_level": 2.0,
        "cc_acceptance_level": 0.7,
        "s2n_limit": 3.,
        "max_time_before_first_arrival": 5.,
        "min_surface_wave_velocity": 1.2,  # Default is 3.0, chow et al.==1.4
        "check_global_data_quality": True,  # Default is False
        "c_0": 0.7,
        "c_1": 2.0,
        "c_3a": 3.0,
        "c_3b": 2.0,
        "c_4a": 2.5,
        "c_4b": 12.0
    },
    # North Island study area, 15-30s bandpass, for long period start waveforms
    "nznorth_15-30s": {
        "stalta_waterlevel": 0.08,
        "tshift_acceptance_level": 12.0,
        "dlna_acceptance_level": 2.5,
        "cc_acceptance_level": 0.7,
        "s2n_limit": 2.5,
        "max_time_before_first_arrival": 10.,
        "min_surface_wave_velocity": 1.2, 
        "check_global_data_quality": True,
        "c_0": 0.7,
        "c_1": 2.0,
        "c_3a": 1.0,
        "c_3b": 2.0,
        "c_4a": 3.0,
        "c_4b": 10.0
    },
    # North Island study area, 10-30s bandpass
    "nznorth_10-30s": {
        "stalta_waterlevel": 0.10, 
        "tshift_acceptance_level": 10.0,  # based on sign-flip
        "dlna_acceptance_level": 2.0,
        "cc_acceptance_level": 0.7,
        "s2n_limit": 3.,
        "max_time_before_first_arrival": 5.,
        "min_surface_wave_velocity": 1.2,  # Default is 3.0, chow et al.==1.4
        "check_global_data_quality": True,  # Default is False
        "c_0": 0.7,
        "c_1": 2.0, 
        "c_3a": 3.0,
        "c_3b": 2.0,
        "c_4a": 2.5,
        "c_4b": 12.0
    },
    # North Island study area, 8-30s bandpass
    # !!! NOT YET TESTED, SIMPLY SCALED FROM 'nznorth_10-30s'
    "nznorth_8-30s": {
        "stalta_waterlevel": 0.10, 
        "tshift_acceptance_level": 8.0,
        "dlna_acceptance_level": 1.5,
        "cc_acceptance_level": 0.7,
        "s2n_limit": 3.,
        "max_time_before_first_arrival": 5.,
        "tshift_reference": 2.,
        "min_surface_wave_velocity": 1.1,
        "check_global_data_quality": True,
        "c_0": 0.7,
        "c_1": 2.5, 
        "c_3a": 4.0,
        "c_3b": 2.0,
        "c_4a": 2.5,
        "c_4b": 12.0
    },
    # North Island study area, 6-30s bandpass
    # !!! NOT YET TESTED, SIMPLY SCALED FROM 'nznorth_10-30s'
    "nznorth_6-30s": {
        "stalta_waterlevel": 0.08,
        "tshift_acceptance_level": 6.0,  
        "dlna_acceptance_level": 1.5,
        "cc_acceptance_level": 0.75,
        "s2n_limit": 3.,
        "max_time_before_first_arrival": 5., 
        "tshift_reference": 2.,
        "min_surface_wave_velocity": 1.1,
        "check_global_data_quality": True,
        "c_0": 0.85,
        "c_1": 3.0, 
        "c_3a": 4.0,
        "c_3b": 2.0,
        "c_4a": 2.5,
        "c_4b": 12.0
    },
    # North Island study area, 3-30s bandpass
    # !!! NOT YET TESTED, LOOSELY BASED ON 'socal_3-30s'
    "nznorth_3-30s": {
        "stalta_waterlevel": 0.07,
        "tshift_acceptance_level": 3.0, 
        "dlna_acceptance_level": 1.,
        "cc_acceptance_level": 0.85,
        "s2n_limit": 4.,
        "max_time_before_first_arrival": 5.,
        "min_surface_wave_velocity": 1.0,
        "check_global_data_quality": True, 
        "c_0": 1.0,
        "c_1": 5.0, 
        "c_3a": 4.0,
        "c_3b": 2.5,
        "c_4a": 2.,
        "c_4b": 6.0
    },
    # For the 1D velocity model of Ristau (2008)
    "nznorth_1D": {
        "stalta_waterlevel": 0.07,
        "tshift_acceptance_level": 10.0,
        "dlna_acceptance_level": 2.0,
        "cc_acceptance_level": 0.7,
        "tshift_reference": 4.,
        "s2n_limit": 3.,
        "max_time_before_first_arrival": 5.,
        "min_surface_wave_velocity": 1.7,  # Default is 3.0
        "check_global_data_quality": True,  # Default is False
        "c_0": 0.7,
        "c_1": 2.,
        "c_3a": 3.0,
        "c_3b": 2.0,
        "c_4a": 2.5,
        "c_4b": 12.0
    },
    # Global scale from Maggi et al. 2009 Table 3 for 50s < T < 150s
    "global": {
        "s2n_limit": 2.5,
        "stalta_waterlevel": 0.08,
        "cc_acceptance_level": 0.85,
        "tshift_acceptance_level": 15.0,
        "dlna_acceptance_level": 1.,
        "tshift_reference": 0.,
        "dlna_reference": 0.,
        "c_0": 0.7,
        "c_1": 4.,
        "c_2": 0.3,
        "c_3a": 1.0,
        "c_3b": 2.0,
        "c_4a": 3.,
        "c_4b": 10.0
    },
    # Japan scale from Maggi et al. 2009 Table 3 for 6 < T < 30s
    "japan_6-30s": {
        "s2n_limit": 3.,
        "stalta_waterlevel": 0.12,
        "cc_acceptance_level": 0.73,
        "tshift_acceptance_level": 3.0,
        "dlna_acceptance_level": 1.5,
        "tshift_reference": 0.,
        "dlna_reference": 0.,
        "c_0": 0.7,
        "c_1": 3.,
        "c_2": 0.6,
        "c_3a": 1.0,
        "c_3b": 2.0,
        "c_4a": 3.,
        "c_4b": 12.0
    },
    # Southern California scale from Maggi et al. 2009 Table 3 for 6 < T < 30s
    "socal_6-30s": {
        "s2n_limit": 3.,
        "stalta_waterlevel": 0.18,
        "cc_acceptance_level": 0.71,
        "tshift_acceptance_level": 8.0,
        "dlna_acceptance_level": 1.5,
        "tshift_reference": 4.,
        "dlna_reference": 0.,
        "c_0": 0.7,
        "c_1": 2.,
        "c_2": 0.,
        "c_3a": 3.0,
        "c_3b": 2.0,
        "c_4a": 2.5,
        "c_4b": 12.0
    },
    # Southern California scale from Maggi et al. 2009 Table 3 for 3 < T < 30s
    "socal_3-30s": {
        "s2n_limit": 4.,
        "stalta_waterlevel": 0.11,
        "cc_acceptance_level": 0.8,
        "tshift_acceptance_level": 4.0,
        "dlna_acceptance_level": 1.,
        "tshift_reference": 2.,
        "dlna_reference": 0.,
        "c_0": 1.3,
        "c_1": 4.,
        "c_2": 0.,
        "c_3a": 4.0,
        "c_3b": 2.5,
        "c_4a": 2.,
        "c_4b": 6.0
    },
    # Southern California scale from Maggi et al. 2009 Table 3 for 2 < T < 30s
    "socal_2-30s": {
        "s2n_limit": 4.,
        "stalta_waterlevel": 0.07,
        "cc_acceptance_level": 0.85,
        "tshift_acceptance_level": 3.0,
        "dlna_acceptance_level": 1.,
        "tshift_reference": 1.,
        "dlna_reference": 0.,
        "c_0": 1.,
        "c_1": 5.,
        "c_2": 0.,
        "c_3a": 4.0,
        "c_3b": 2.5,
        "c_4a": 2.,
        "c_4b": 6.0
    },
    # From the UAF group doing regional studies of Alaska
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
}
