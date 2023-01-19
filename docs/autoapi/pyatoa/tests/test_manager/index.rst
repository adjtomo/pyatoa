:py:mod:`pyatoa.tests.test_manager`
===================================

.. py:module:: pyatoa.tests.test_manager

.. autoapi-nested-parse::

   Test the functionalities of the Pyaflowa Manager class



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_manager.st_obs
   pyatoa.tests.test_manager.st_syn
   pyatoa.tests.test_manager.cat
   pyatoa.tests.test_manager.event
   pyatoa.tests.test_manager.inv
   pyatoa.tests.test_manager.config
   pyatoa.tests.test_manager.mgmt_pre
   pyatoa.tests.test_manager.mgmt_post
   pyatoa.tests.test_manager.test_read_write_from_asdfdataset
   pyatoa.tests.test_manager.test_standardize_to_synthetics
   pyatoa.tests.test_manager.test_standardize_raises_manager_error
   pyatoa.tests.test_manager.test_preprocess_dont_rotate
   pyatoa.tests.test_manager.test_preprocess_rotate_to_rtz
   pyatoa.tests.test_manager.test_preprocess_overwrite
   pyatoa.tests.test_manager.test_select_window
   pyatoa.tests.test_manager.test_save_and_retrieve_windows
   pyatoa.tests.test_manager.test_save_adjsrcs
   pyatoa.tests.test_manager.test_format_windows
   pyatoa.tests.test_manager.test_flow_multiband



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_manager.propogate


.. py:data:: propogate
   :annotation: = False

   

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


.. py:function:: test_read_write_from_asdfdataset(tmpdir, mgmt_pre, config)

   Write a Manager into an ASDFDataSet and then read it back


.. py:function:: test_standardize_to_synthetics(mgmt_pre)

   Ensure that standardizing streams performs three main tasks, trimming
   origin times, matching sampling rates, and matching number of points.


.. py:function:: test_standardize_raises_manager_error(mgmt_pre)

   Asser that Manager will raise an error if user tries to standardize with
   no traces present


.. py:function:: test_preprocess_dont_rotate(mgmt_pre)

   Standard preprocessing, dont rotate components, just filter, check if
   filtering worked


.. py:function:: test_preprocess_rotate_to_rtz(mgmt_pre)

   Standard preprocessing but rotate components based on the backazimuth


.. py:function:: test_preprocess_overwrite(mgmt_pre)

   Apply an overwriting preprocessing function to ensure functionality works


.. py:function:: test_select_window(mgmt_pre)

   Ensure windows functionality works as advertised


.. py:function:: test_save_and_retrieve_windows(tmpdir, mgmt_post)

   Test retrieve_windows() and save_windows() by saving windows into a
   scratch dataset and retrieving them back. Window criteria will be
   recalculated but since the waveforms are the same, the values will be the
   same as before.


.. py:function:: test_save_adjsrcs(tmpdir, mgmt_post)

   Checks that adjoint sources can be written to dataset and will match the
   formatting required by Specfem3D


.. py:function:: test_format_windows(mgmt_post)

   Basic check that format windows returns as formatted lists expected


.. py:function:: test_flow_multiband(mgmt_pre)

   Test that the workflow for multiple period bands returns a single
   adjoint source


