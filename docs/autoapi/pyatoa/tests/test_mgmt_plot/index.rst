:py:mod:`pyatoa.tests.test_mgmt_plot`
=====================================

.. py:module:: pyatoa.tests.test_mgmt_plot

.. autoapi-nested-parse::

   Test plotting capabilities by comparing baseline images for example data



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_mgmt_plot.mgmt
   pyatoa.tests.test_mgmt_plot.reset_matplotlib
   pyatoa.tests.test_mgmt_plot.images_are_identical
   pyatoa.tests.test_mgmt_plot.test_waveform_plot
   pyatoa.tests.test_mgmt_plot.test_map_plot
   pyatoa.tests.test_mgmt_plot.test_combined_plot



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_mgmt_plot.propogate
   pyatoa.tests.test_mgmt_plot.IMGDIR


.. py:data:: propogate
   :annotation: = False

   

.. py:data:: IMGDIR
   :annotation: = ./test_data/baseline_images

   

.. py:function:: mgmt()

   Fully realized manager object which will be used to make plots


.. py:function:: reset_matplotlib()

   Reset matplotlib to a common default. Copied from Pyflex


.. py:function:: images_are_identical(created_img)

   Modified from Pyflex which was partially copied from ObsPy.
   Used to check images for equality.


.. py:function:: test_waveform_plot(tmpdir, mgmt)

   Test that plotting waveforms, windows, adjsrc etc. by themselves works


.. py:function:: test_map_plot(tmpdir, mgmt)

   Test that map plotting with moment tensor works


.. py:function:: test_combined_plot(tmpdir, mgmt)

   Test that a combined gridspec waveform + map figure works


