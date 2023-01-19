:py:mod:`pyatoa.tests.test_utils`
=================================

.. py:module:: pyatoa.tests.test_utils

.. autoapi-nested-parse::

   Test all the various utility functions except for the processing functions,
   which get their own test suite.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_utils.station_fid
   pyatoa.tests.test_utils.sem_fid
   pyatoa.tests.test_utils.ds
   pyatoa.tests.test_utils.test_traveltime_adjoint_source
   pyatoa.tests.test_utils.test_calc_utils
   pyatoa.tests.test_utils.test_form_utils
   pyatoa.tests.test_utils.test_images_utils
   pyatoa.tests.test_utils.test_write_utils



.. py:function:: station_fid()

   Re-used test data pointing to STATIONS file


.. py:function:: sem_fid(tmpdir)

   Re-used test data pointing to a two-column ascii file representing
   SPECFEM3D seismograms


.. py:function:: ds()

   Test ASDFDataSet


.. py:function:: test_traveltime_adjoint_source()

   !!! TO DO


.. py:function:: test_calc_utils()

   Wrapping all the calculate utilities into a single function because they're
   all pretty basic, but it's good to know they function the same each time


.. py:function:: test_form_utils()

   Test the string formatting utilities which ensures that IO is standardized
   throughout the package. Again, place them all in the same function.


.. py:function:: test_images_utils()

   Test image processing utilities


.. py:function:: test_write_utils(tmpdir, station_fid, sem_fid, ds)

   Test write utilities


