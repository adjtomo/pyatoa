"""
Additional capabilities for creating adjoint sources not defined by PyAdjoint
"""
import numpy as np
from scipy.integrate import simps


def traveltime_adjoint_source(tr_syn, t_start=0, t_end=-1, save=False):
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
    :type t_start: float
    :param t_start: optional start time to window data on e.g. specfic phase
    :type t_end: float
    :param t_end: optional end time to window data on e.g. specfic phase
    :rtype: np.array
    :return: a numpy array that defines the adjoint source
    """
    s = tr_syn.data
    deltat = tr_syn.stats.delta

    # Generate the adjoint source 'fp'
    dsdt = np.gradient(s, deltat)
    nnorm = simps(y=dsdt * dsdt, dx=deltat)
    fp = 1 / nnorm * dsdt[:]

    return fp[::-1]




