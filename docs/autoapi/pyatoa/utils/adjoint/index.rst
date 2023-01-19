:py:mod:`pyatoa.utils.adjoint`
==============================

.. py:module:: pyatoa.utils.adjoint

.. autoapi-nested-parse::

   Additional capabilities for creating adjoint sources not defined by PyAdjoint



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.adjoint.traveltime_adjoint_source



.. py:function:: traveltime_adjoint_source(tr, time_window=None, reverse=True, save=False, zeros=False)

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


