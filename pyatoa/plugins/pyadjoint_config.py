"""
Overwrite some Pyadjoint configurations. Config parameters can be found at:

https://github.com/krischer/pyadjoint/blob/master/src/pyadjoint/config.py
"""
from pyadjoint import Config as pyadjointConfig


def set_pyadjoint_config(min_period, max_period, **kwargs):
    """
    Set the Pyadjoint config based on Pyatoa config parameter.
    Kwargs can be fed to the Pyadjoint Config object.
    Return unnused kwargs.

    :type min_period: float
    :param min_period: min period of the data
    :type max_period: float
    :param max_period: max period of the data
    :rtype cfgout: pyadjoint.Config
    :return cfgout: properly set pyadjoint configuration object
    """
    paconfig = pyadjointConfig(min_period=min_period,
                               max_period=max_period
                               )

    unused_kwargs = []
    for key, item in kwargs.items():
        if hasattr(paconfig, key):
            setattr(paconfig, key, item)
        else:
            unused_kwargs.append(key)

    return paconfig, unused_kwargs


def src_type(choice):
    """
    Pyadjoint requires that a standard adjoint source type is given in the
    calculate_adjoint_source() function. This function acts as simple dictionary
    input to provide the correct input for that function

    :type choice: str
    :param choice: pyatoa.Config.adj_src_type
    :rtype: str
    :return: pyadjoint adj_src_type
    """
    if "cc" in choice:
        adj_src_type = "cc_traveltime_misfit"
    elif "mt" in choice:
        adj_src_type = "multitaper_misfit"
    else:
        adj_src_type = "waveform"
    return adj_src_type


def set_pyadjoint_config_devel(choice, min_period, max_period):
    """
    Set the Pyadjoint config based on Pyatoa config parameter

    For available adjoint source types, see:
    https://github.com/computational-seismology/
                                pyadjoint/blob/dev/src/pyadjoint/config.py

    Available but not listed here:
    ExponentiatedPhase, DoubleDifferenceCrossCorrelation,
    DoubleDifferenceMultiTaper,

    :type choice: str
    :param choice: type of adjoint source to use
    :type min_period: float
    :param min_period: min period of the data
    :type max_period: float
    :param max_period: max period of the data
    :rtype cfgout: pyadjoint.Config*
    :return cfgout: properly set pyadjoint configuration object
    """
    import pyadjoint

    if "waveform" in choice:
        cfgout = pyadjoint.ConfigWaveForm(min_period=min_period,
                                          max_period=max_period)
    elif "cc" in choice:
        cfgout = pyadjoint.ConfigCrossCorrelation(min_period=min_period,
                                                  max_period=max_period)
        # custom cc_traveltime_misfit configuration
        if choice == "cc_hikurangi":
            setattr(cfgout, "use_cc_error", False)  # default = True
    # The default multitaper parameter min_cycle_in_windows=0.5 is wrong,
    # this sets it back to the proper default of 3
    elif "multitaper" in choice:
        cfgout = pyadjoint.ConfigMultiTaper(min_period=min_period,
                                            max_period=max_period,
                                            min_cycle_in_window=3)
        # custom mtm configuration
        if choice == "multitaper_hikurangi":
            setattr(cfgout, "ipower_costaper", 8)  # default = 10
            setattr(cfgout, "min_cycle_in_window", 0)  # default = 3
            setattr(cfgout, "taper_percentage", 0.5)  # default = 0.3
            setattr(cfgout, "use_cc_error", False)  # default = True
            setattr(cfgout, "use_mt_error", False)  # default = True
    else:
        raise KeyError("adjoint source type incorrectly specified, "
                       "must be: 'waveform', 'cc_traveltime_misfit', "
                       "'multitaper_misfit'")

    return cfgout
