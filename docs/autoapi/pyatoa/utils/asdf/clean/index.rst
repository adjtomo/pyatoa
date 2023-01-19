:py:mod:`pyatoa.utils.asdf.clean`
=================================

.. py:module:: pyatoa.utils.asdf.clean

.. autoapi-nested-parse::

   Convenience functions for removing data from Pyasdf ASDFDataSet objects.
   All functions work with the dataset as an input and act in-place on the
   dataset so no returns



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.asdf.clean.clean_dataset
   pyatoa.utils.asdf.clean.del_synthetic_waveforms
   pyatoa.utils.asdf.clean.del_auxiliary_data



.. py:function:: clean_dataset(ds, iteration=None, step_count=None, fix_windows=False)

   Removes synthetic waveforms and auxiliary data so that only observation
   data remains for new iterations. If no iteration is given, will
   wipe all non-observation data and all auxiliary data

   :type ds: pyasdf.ASDFDataSet
   :param ds: dataset to be cleaned
   :type fix_windows: bool
   :param fix_windows: don't delete MisfitWindows
   :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
   :type step_count: str or int
   :param step_count: step count e.g. "s00". Will be formatted so int ok.


.. py:function:: del_synthetic_waveforms(ds, iteration=None, step_count=None)

   Remove "synthetic_{iter_tag}{step_tag}" tagged waveforms from an asdf
   dataset. If no iter_tag number given, wipes all synthetic data from dataset.

   :type ds: pyasdf.ASDFDataSet
   :param ds: dataset to be cleaned
   :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
   :type step_count: str or int
   :param step_count: step count e.g. "s00". Will be formatted so int ok.


.. py:function:: del_auxiliary_data(ds, iteration=None, step_count=None, retain=None, only=None)

   Delete all items in auxiliary data for a given iter_tag, if iter_tag not
   given, wipes all auxiliary data.

   :type ds: pyasdf.ASDFDataSet
   :param ds: dataset to be cleaned
   :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
   :type step_count: str or int
   :param step_count: step count e.g. "s00". Will be formatted so int ok.
   :type retain: list of str
   :param retain: list of auxiliary data tags to retain, that is: delete all
       auxiliary data EXCEPT FOR the names given in this variable
   :type only: list of str
   :param only: list of auxiliary data tags to remove, that is: ONLY delete
       auxiliary data that matches the names given in this variable.
       Lower in priority than 'retain'


