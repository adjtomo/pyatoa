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


def set_pyadjoint_config(choice, min_period, max_period):
    """
    Set the Pyadjoint config based on Pyatoa config parameter

    For available adjoint source types, see:
    https://github.com/computational-seismology/pyadjoint/blob/dev/src/pyadjoint/config.py

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
    # Default parameters
    if choice == "waveform":
        cfgout = (choice, pyadjoint.ConfigWaveForm(min_period=min_period,
                                                   max_period=max_period)
                  )
    elif choice == "cc_traveltime_misfit":
        cfgout = (choice, pyadjoint.ConfigCrossCorrelation(
            min_period=min_period, max_period=max_period)
                  )
    elif choice == "multitaper_misfit":
        cfgout = (choice, pyadjoint.ConfigMultiTaper(min_period=min_period,
                                                     max_period=max_period,
                                                     min_cycle_in_window=3
                                                     )
                  )
    # Custom configurations
    elif choice == "mtm_hikurangi":
        cfgout = ("multitaper_misfit",
                  pyadjoint.ConfigMultiTaper(
                      min_period=min_period,
                      max_period=max_period,
                      lnpt=15,
                      transfunc_waterlevel=1e-10,
                      water_threshold=0.02,
                      ipower_costaper=2,  # default=10, 2=broader smoothing
                      min_cycle_in_window=0,  # default 0.5 but that's incorrect
                      taper_type='hann',
                      taper_percentage=0.5,  # % of window to taper at each end
                      mt_nw=4.0,  # bin width of multitapers
                      num_taper=5,  # number of tapers
                      dt_fac=2.0,
                      phase_step=1.5,
                      err_fac=2.5,
                      dt_max_scale=3.5,
                      measure_type='dt',
                      dt_sigma_min=1.0,
                      dlna_sigma_min=0.5,
                      use_cc_error=False,
                      use_mt_error=False)
                  )
    elif choice == "cc_hikurangi":
        cfgout = ("cc_traveltime_misfit",
                  pyadjoint.ConfigCrossCorrelation(
                      min_period=min_period,
                      max_period=max_period,
                      taper_type='hann',  # hann
                      measure_type='dt',  # dt
                      use_cc_error=False,  # True
                      dt_sigma_min=1.0,  # 1.
                      dlna_sigma_min=0.5)  # .5
                  )
    else:
        raise KeyError("adjoint source type incorrectly specified, "
                       "must be: 'waveform', 'cc_traveltime_misfit', "
                       "'multitaper_misfit'")

    return cfgout
