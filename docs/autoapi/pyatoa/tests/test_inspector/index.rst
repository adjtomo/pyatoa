:py:mod:`pyatoa.tests.test_inspector`
=====================================

.. py:module:: pyatoa.tests.test_inspector

.. autoapi-nested-parse::

   Test the Inspector class and its ability to generate dataframes for bulk
   analyses of an inversion



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_inspector.test_data
   pyatoa.tests.test_inspector.asdf_dataset_fid
   pyatoa.tests.test_inspector.seisflows_inspector
   pyatoa.tests.test_inspector.inspector
   pyatoa.tests.test_inspector.test_append
   pyatoa.tests.test_inspector.test_discover
   pyatoa.tests.test_inspector.test_extend
   pyatoa.tests.test_inspector.test_read_write_csv
   pyatoa.tests.test_inspector.test_isolate
   pyatoa.tests.test_inspector.test_nwin
   pyatoa.tests.test_inspector.test_misfit
   pyatoa.tests.test_inspector.test_stats
   pyatoa.tests.test_inspector.test_minmax
   pyatoa.tests.test_inspector.test_compare
   pyatoa.tests.test_inspector.test_compare_no_data
   pyatoa.tests.test_inspector.test_compare_windows
   pyatoa.tests.test_inspector.test_no_step_information
   pyatoa.tests.test_inspector.test_get_models
   pyatoa.tests.test_inspector.test_get_unique_models



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_inspector.propagate


.. py:data:: propagate
   :annotation: = False

   

.. py:function:: test_data()


.. py:function:: asdf_dataset_fid()


.. py:function:: seisflows_inspector(test_data)

   An Inspector filled in by a SeisFlows workflow used to test some of the
   more complex functions


.. py:function:: inspector(asdf_dataset_fid)


.. py:function:: test_append(asdf_dataset_fid)

   Make sure that the Inspector can append a single event and put the correct
   data in the correct place


.. py:function:: test_discover(test_data)

   Make sure Inspector can find HDF5 files generally and read them in.


.. py:function:: test_extend(inspector)

   Make sure that you can extend the current inspector with the widnows of
   another. In this example we just use the same inspector twice
   Also tests the copy function to make sure that the extension doesn't
   affect both inspectors


.. py:function:: test_read_write_csv(tmpdir, inspector)

   Test the read and write functions with both an empty and a filled inspector


.. py:function:: test_isolate(seisflows_inspector)

   Test the isolate function to grab specific data from a filled inspector.
   Test by checking number of windows at each isolation call


.. py:function:: test_nwin(seisflows_inspector)

   Test the number of windows function


.. py:function:: test_misfit(seisflows_inspector)

   Test the misfit calculation function


.. py:function:: test_stats(seisflows_inspector)

   Test the per-level stats calculations


.. py:function:: test_minmax(seisflows_inspector)

   Test the minmax printing function


.. py:function:: test_compare(seisflows_inspector)

   Test inter-event comparisons


.. py:function:: test_compare_no_data()

   Test that compare with no data returns NoneType


.. py:function:: test_compare_windows(seisflows_inspector)

   TODO Need fixed window inversion results to make this work


.. py:function:: test_no_step_information(seisflows_inspector)

   TODO Test that when no step information is provided (only iteration),
   TODO inspector can still handle data


.. py:function:: test_get_models(seisflows_inspector)

   Test the model state tracker


.. py:function:: test_get_unique_models(seisflows_inspector)

   Test convenience function that finds accepted models only


