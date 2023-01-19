:py:mod:`pyatoa.tests.test_asdf_utils`
======================================

.. py:module:: pyatoa.tests.test_asdf_utils

.. autoapi-nested-parse::

   Test the PyASDF ASDFDataSet auxiliary utilities

   !!! TO DO: utils.asdf.write tests



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_asdf_utils.empty_dataset
   pyatoa.tests.test_asdf_utils.dataset
   pyatoa.tests.test_asdf_utils.st_obs
   pyatoa.tests.test_asdf_utils.st_syn
   pyatoa.tests.test_asdf_utils.cat
   pyatoa.tests.test_asdf_utils.event
   pyatoa.tests.test_asdf_utils.inv
   pyatoa.tests.test_asdf_utils.config
   pyatoa.tests.test_asdf_utils.mgmt_pre
   pyatoa.tests.test_asdf_utils.mgmt_post
   pyatoa.tests.test_asdf_utils.test_add_misfit_windows
   pyatoa.tests.test_asdf_utils.test_add_adjoint_sources
   pyatoa.tests.test_asdf_utils.test_clean_dataset
   pyatoa.tests.test_asdf_utils.test_clean_dataset_fix_windows
   pyatoa.tests.test_asdf_utils.test_load_windows
   pyatoa.tests.test_asdf_utils.test_load_previous_windows
   pyatoa.tests.test_asdf_utils.test_load_adjsrcs



.. py:function:: empty_dataset(tmpdir)

   Re-used test data pointing to STATIONS file


.. py:function:: dataset()

   Filled ASDFDataSet


.. py:function:: st_obs()

   Raw observed waveforms from station NZ.BFZ.HH? for New Zealand event
   2018p130600 (GeoNet event id)


.. py:function:: st_syn()

   Synthetic data stream generated using Specfem3D and focal mechanism for
   2018p130600. Minimum resolved period roughly 10s.


.. py:function:: cat()

   ObsPy Event Catalog for New Zealand based event with
   GeoNet Event ID: 2018p130600


.. py:function:: event(cat)

   Event from Catalog


.. py:function:: inv()

   StationXML information for station NZ.BFZ.HH?


.. py:function:: config()

   Default Pyatoa Config object


.. py:function:: mgmt_pre(config, event, st_obs, st_syn, inv)

   A manager filled with data but pre-workflow


.. py:function:: mgmt_post(mgmt_pre)

   A manager that has completed the full workflow


.. py:function:: test_add_misfit_windows(empty_dataset, mgmt_post)

   Test adding misfit windows to an ASDFDataSet
   :return:


.. py:function:: test_add_adjoint_sources(empty_dataset, mgmt_post)

   Test adding adjoint sources to an ASDFDataSet
   :return:


.. py:function:: test_clean_dataset(empty_dataset, mgmt_pre)

   Test dataset clean functions. Need to perform tasks on a dataset we create
   here, otherwise we may permanently affect test data if we use a pre-built
   ASDFDataSet
   :return:


.. py:function:: test_clean_dataset_fix_windows(empty_dataset, mgmt_pre)

   Test cleaning a dataset but retaining windows
   :return:


.. py:function:: test_load_windows(dataset)

   Test the function that returns windows in the Pyflex output format from
   an ASDFDataSet


.. py:function:: test_load_previous_windows(empty_dataset, mgmt_pre)

   Test the function that returns windows in the Pyflex output format from
   an ASDFDataSet

   :param dataset:
   :return:


.. py:function:: test_load_adjsrcs(dataset)

   Test the function that returns windows in the Pyflex output format from
   an ASDFDataSet

   :param dataset:
   :return:


