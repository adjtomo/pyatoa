:py:mod:`pyatoa.tests.test_process_util`
========================================

.. py:module:: pyatoa.tests.test_process_util

.. autoapi-nested-parse::

   Test the functionalities of the suite of processing functions in the utilities



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_process_util.st_obs
   pyatoa.tests.test_process_util.st_syn
   pyatoa.tests.test_process_util.cat
   pyatoa.tests.test_process_util.event
   pyatoa.tests.test_process_util.inv
   pyatoa.tests.test_process_util.config
   pyatoa.tests.test_process_util.mgmt_pre
   pyatoa.tests.test_process_util.test_filters
   pyatoa.tests.test_process_util.test_zero_pad
   pyatoa.tests.test_process_util.test_trim_streams_and_match_npts
   pyatoa.tests.test_process_util.test_is_preprocessed
   pyatoa.tests.test_process_util.test_stf_convolve



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


.. py:function:: test_filters(st_obs)

   Make sure st_obs works with various input frequencies/ periods


.. py:function:: test_zero_pad(st_obs)

   Ensure that zero padding adds the same number of data points each time


.. py:function:: test_trim_streams_and_match_npts(st_obs, st_syn)

   Ensure that forcing number of points standardization works


.. py:function:: test_is_preprocessed(st_obs)

   Test the check function that determines if a stream is preprocessed


.. py:function:: test_stf_convolve()

   !!! TO DO


