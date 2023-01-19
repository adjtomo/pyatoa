:py:mod:`pyatoa.utils.asdf.load`
================================

.. py:module:: pyatoa.utils.asdf.load

.. autoapi-nested-parse::

   Functions for extracting information from a Pyasdf ASDFDataSet object



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.asdf.load.load_windows
   pyatoa.utils.asdf.load.load_adjsrcs
   pyatoa.utils.asdf.load.dataset_windows_to_pyflex_windows
   pyatoa.utils.asdf.load.previous_windows



.. py:function:: load_windows(ds, net, sta, iteration, step_count, return_previous=False)

   Returns misfit windows from an ASDFDataSet for a given iteration, step,
   network and station, as well as a count of windows returned.

   If given iteration and step are not present in dataset (e.g. during line
   search, new step), will try to search the previous step, which may or
   may not be contained in the previous iteration.

   Returns windows as Pyflex Window objects which can be used in Pyadjoint or
   in the Pyatoa workflow.

   .. note::
       Expects that windows are saved into the dataset at each iteration and
       step such that there is a coherent structure within the dataset

   :type ds: pyasdf.ASDFDataSet
   :param ds: ASDF dataset containing MisfitWindows subgroup
   :type net: str
   :param net: network code used to find the name of the misfit window
   :type sta: str
   :param sta: station code used to find the name of the misfit window
   :type iteration: int or str
   :param iteration: current iteration, will be formatted by the function
   :type step_count: int or str
   :param step_count: step count, will be formatted by the function
   :type return_previous: bool
   :param return_previous: search the dataset for available windows
       from the previous iteration/step given the current iteration/step
   :rtype window_dict: dict
   :return window_dict: dictionary containing misfit windows, in a format
       expected by Pyatoa Manager class


.. py:function:: load_adjsrcs(ds, net, sta, iteration, step_count)

   Load adjoint sources from a pyasdf ASDFDataSet and return in the format
   expected by the Manager class, that is a dictionary of adjoint sources

   :type ds: pyasdf.ASDFDataSet
   :param ds: ASDF dataset containing MisfitWindows subgroup
   :type net: str
   :param net: network code used to find the name of the adjoint source
   :type sta: str
   :param sta: station code used to find the name of the adjoint source
   :type iteration: int or str
   :param iteration: current iteration, will be formatted by the function
   :type step_count: int or str
   :param step_count: step count, will be formatted by the function
   :rtype: dict
   :return: dictionary containing adjoint sources, in a format expected by
       Pyatoa Manager class


.. py:function:: dataset_windows_to_pyflex_windows(windows, network, station)

   Convert the parameter dictionary of an ASDFDataSet MisfitWindow into a
   dictionary of Pyflex Window objects, in the same format as Manager.windows

   Returns empty dict and 0 if no windows are found

   :type windows: pyasdf.utils.AuxiliaryDataAccessor
   :param windows: ds.auxiliary_data.MisfitWindows[iter][step]
   :type network: str
   :param network: network of the station related to the windows
   :type station: str
   :param station: station related to the windows
   :rtype: dict
   :return: dictionary of window attributes in the same format that Pyflex
       outputs


.. py:function:: previous_windows(windows, iteration, step_count)

   Given an iteration and step count, find windows from the previous step
   count. If none are found for the given iteration, return the most recently
   available windows.

   .. note::
       Assumes that windows are saved at each iteration.

   :type windows: pyasdf.utils.AuxiliaryDataAccessor
   :param windows: ds.auxiliary_data.MisfitWindows[iter][step]
   :type iteration: int or str
   :param iteration: the current iteration
   :type step_count: int or str
   :param step_count: the current step count
   :rtype: pyasdf.utils.AuxiliaryDataAccessor
   :return: ds.auxiliary_data.MisfitWindows


