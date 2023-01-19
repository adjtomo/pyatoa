:py:mod:`pyatoa.utils.window`
=============================

.. py:module:: pyatoa.utils.window

.. autoapi-nested-parse::

   Auxiliary criteria for weighting misfit windows to enhance or suppress certain
   measurements in order to fairly assess misfit criteria. This is meant to
   compliment the functionalities of Pyflex without having to directly edit the
   Pyflex source code.

   Functions should work in place on a Manager class to avoid having to pass in
   all the different arguments from the Manager.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.window.zero_pad_then_window
   pyatoa.utils.window.reject_on_global_amplitude_ratio



.. py:function:: zero_pad_then_window(ws, pad_by_fraction_of_npts=0.2)

   To address Pyflex throwing ValueErrors when source-receiver distances are

   .. note::
       Sept 1, 2020
       Work in progress, may not actually want to do this to avoid any
       near-source effects?

   :type ws: pyflex.WindowSelector
   :param ws: an already-filled window selector object that should
       be passed in from the Manager object
   :rtype: list of pyflex.Window
   :return: a list of Window objects, or an empty list if no windows found or
       the zero padding didnt work


.. py:function:: reject_on_global_amplitude_ratio(data, windows, ratio=0.2)

   Reject windows where peak amplitude falls below some threshold value.

   This was created in order to suppress windows containing long period direct
   arrivals, which were creating high-frequency adjoint sources.

   :type data: np.ndarray
   :param data: data array to query amplitude values from
   :type windows: list of pyflex.window.Window
   :param windows: list of window objects to check
   :type ratio: float
   :param ratio: percentage threshold of the peak value within a given window
       and the global peak value in the data array. Defaults to 0.2
   :rtype: tuple of lists of pyflex.window.Window
   :return: lists of accepted and rejected windows


