"""
Pyflex requires configuration objects to run;
the functions here will set them properly and return the necessary parameters

For variable descriptions see:
    https://krischer.github.io/pyflex/_modules/pyflex/config.html
"""
from pyflex import Config as pyflexConfig


# Some preset Pyflex Configuration parameters based on previous work,
# current work or literature.
presets = {
    # example configuration from the Pyflex website, different from default
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
    # from the UAF group doing local studies of Alaska
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
    # From the New Zealand group doing local studies of North Island, 10-30s
    "hikurangi_10-30s": {
        "stalta_waterlevel": 0.18,
        "tshift_acceptance_level": 8.0,  # based on sign-flip
        "dlna_acceptance_level": 1.5,
        "cc_acceptance_level": 0.7,
        "s2n_limit": 3.,
        "max_time_before_first_arrival": 0.,  # minimum starting before P-wave
        "c_0": 0.7,
        "c_1": 2.5,  # min_win_len = c_1 * min_period
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
}


def set_pyflex_config(min_period, max_period, choice=None, **kwargs):
    """
    Overwriting the default Pyflex parameters with User-defined criteria

    :type choice: str or dict
    :param choice: name of map to choose the Pyflex config options, if None,
        default values are used. Kwargs can still overload default values.
        Also dicts can be passed in as User-defined preset
    :type min_period: float
    :param min_period: min period of the data
    :type max_period: float
    :param max_period: max period of the data
    :rtype: pyflex.Config
    :return: the pyflex Config option to use when running Pyflex
    """
    # Instantiate the pyflex Config object
    pfconfig = pyflexConfig(min_period=min_period, max_period=max_period)

    # Set preset configuration parameters based on hard-coded presets
    if isinstance(choice, str):
        if choice in presets.keys():
            preset = presets[choice]
            for key, item in preset.items():
                setattr(pfconfig, key, item)
        else:
            raise KeyError(f"'{choice}' does not match any available presets "
                           f"for Pyflex. "
                           f"Presets include {list(presets.keys())}")
    # Allow dictionary object to be passed in as a preset
    elif isinstance(choice, dict):
        for key, item in choice.items():
            setattr(pfconfig, key, item)

    # Kwargs can also be passed from the pyatoa.Config object to avoid having to
    # define pre-set values. Kwargs will override preset values
    unused_kwargs = []
    for key, item in kwargs.items():
        if hasattr(pfconfig, key):
            setattr(pfconfig, key, item)
        else:
            unused_kwargs.append(key)

    return pfconfig, unused_kwargs


