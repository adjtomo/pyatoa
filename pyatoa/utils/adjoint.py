"""
Additional capabilities for creating adjoint sources not defined by PyAdjoint
"""
import numpy as np
from scipy.integrate import simps
from scipy.signal.windows import tukey


def traveltime_adjoint_source(tr, time_window=None, reverse=True, save=False,
                              zeros=False):
    """
    Define a traveltime adjoint source, used to generate 'Banana-doughtnut'
    kernels. Traveltime adjoint sources are not data dependent, but rather they
    are sensitivity kernels that illuminate the finite-frequency ray path
    of the waveforms.

    Equation and variable naming is based on Tromp et al. (2005) Eq. 45.
    Implementation is based on Pyadjoint's cc_traveltime adjoint source
    Tapering is done with a hanning window.

    :type st_syn: obspy.core.trace.Trace
    :param st_syn: Synthetic data to be converted to traveltime adjoint source
    :type t_window: list of float
    :param t_window: [t_start, t_end] window to cut phases from waveform
    :rtype: np.array
    :return: a numpy array that defines the adjoint source
    """
    s = tr.data
    deltat = tr.stats.delta
    offset = float(tr.stats.starttime)
    times = tr.times() + offset

    # Generate the adjoint source 'fp'
    dsdt = np.gradient(s, deltat)
    nnorm = simps(y=dsdt * dsdt, dx=deltat)
    fp = 1 / nnorm * dsdt[:]

    if zeros:
        # Only write zeros for empty adjoint sources
        fp = np.zeros(tr.stats.npts)
    else:
        # Window the adjoint source based on given start and end times
        if time_window is not None:
            t_start, t_end = time_window
            overlay = np.zeros(tr.stats.npts)
            samp_start = int((t_start - offset) / deltat)
            samp_end = int((t_end - offset) / deltat)
            
            window = tukey(samp_end - samp_start, alpha=0.5)
            overlay[samp_start:samp_end] = window

            fp = np.multiply(overlay, fp) 

        # Adjoint sources need to be time reversed
        if reverse:
            fp = fp[::-1]

    data = np.vstack((times, fp)).T

    if save:
        np.savetxt(save, data, "%14.7f %20.8E")
        
    return data



