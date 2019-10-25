"""
Here we overwrite some Pyadjoint configurations

Pyadjoint Multitaper Config

    :param min_period: Minimum period of the filtered input data in seconds.
    :type min_period: float
    :param max_period: Maximum period of the filtered input data in seconds.
    :type max_period: float
    :param lnpt: power index to determin the time lenght use in FFT (2^lnpt)
    :type lnpt: int
    :param transfunc_waterlevel: Water level on the transfer function
        between data and synthetic.
    :type transfunc_waterlevel: float
    :param ipower_costaper: order of cosine taper, higher the value,
        steeper the shoulders.
    :type ipower_costaper: int
    :param min_cycle_in_window:  Minimum cycle of a wave in time window to
        determin the maximum period can be reliably measured.
    :type min_cycle_in_window: int
    :param taper_percentage: Percentage of a time window needs to be
        tapered at two ends, to remove the non-zero values for adjoint
        source and for fft.
    :type taper_percentage: float
    :param taper_type: Taper type, supports "hann", "cos", "cos_p10" so far
    :type taper_type: str
    :param mt_nw: bin width of multitapers (nw*df is the the half
        bandwidth of multitapers in frequency domain,
        typical values are 2.5, 3., 3.5, 4.0)
    :type mt_nw: float
    :param num_taper: number of eigen tapers (2*nw - 3 gives tapers
        with eigen values larger than 0.96)
    :type num_taper: int
    :param dt_fac
    :type dt_fac: float
    :param err_fac
    :type err_fac: float
    :param dt_max_scale
    :type dt_max_scale: float
    :param phase_step: maximum step for cycle skip correction (?)
    :type phase_step: float
    :param dt_sigma_min: minimum travel time error allowed
    :type dt_sigma_min: float
    :param dlna_sigma_min: minimum amplitude error allowed
    :type dlna_sigma_min: float
    :param measure_type: type of measurements:
                            dt(travel time),
                            am(dlnA),
                            wf(full waveform)
    :param measure_type: string
    :param use_cc_error: use cross correlation errors for
    :type use_cc_error: logic
    :param use_mt_error: use multi-taper error
    :type use_mt_error: logic
"""
import pyadjoint


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
    elif "multitaper" in choice:
        adj_src_type = "multitaper_misfit"
    else:
        adj_src_type = "waveform"
    return adj_src_type


def set_pyadjoint_config(choice, min_period, max_period):
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
