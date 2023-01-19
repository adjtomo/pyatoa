:py:mod:`pyatoa.scripts.load_example_data`
==========================================

.. py:module:: pyatoa.scripts.load_example_data

.. autoapi-nested-parse::

   It's useful to generate a fully loaded Manager object for testing purposes.
   Load data from the test directory and run the Manager workflow to achieve this.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.scripts.load_example_data.load_example_data
   pyatoa.scripts.load_example_data.load_example_inspector
   pyatoa.scripts.load_example_data.load_example_asdfdataset
   pyatoa.scripts.load_example_data.generate_example_asdfdataset



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.scripts.load_example_data._root_dir
   pyatoa.scripts.load_example_data._test_data_dir


.. py:data:: _root_dir
   

   

.. py:data:: _test_data_dir
   

   

.. py:function:: load_example_data()

   Returns example test data that can be used to preload a Manager


.. py:function:: load_example_inspector()

   Returns an example Inspector loaded with some waveform measurements


.. py:function:: load_example_asdfdataset()

   Returns an example ASDFDataSet with a few waveforms etc


.. py:function:: generate_example_asdfdataset()

   Create the test_ASDFDataSet file, which sometimes needs to be re-made if
   the package or dependencies change


