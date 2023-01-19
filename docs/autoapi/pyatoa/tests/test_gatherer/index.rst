:py:mod:`pyatoa.tests.test_gatherer`
====================================

.. py:module:: pyatoa.tests.test_gatherer

.. autoapi-nested-parse::

   Test the functionalities of the Pyatoa Gatherer class



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_gatherer.code
   pyatoa.tests.test_gatherer.event_id
   pyatoa.tests.test_gatherer.dataset_fid
   pyatoa.tests.test_gatherer.cat
   pyatoa.tests.test_gatherer.event
   pyatoa.tests.test_gatherer.origintime
   pyatoa.tests.test_gatherer.config
   pyatoa.tests.test_gatherer.gatherer
   pyatoa.tests.test_gatherer.test_asdf_event_fetch
   pyatoa.tests.test_gatherer.test_asdf_station_fetch
   pyatoa.tests.test_gatherer.test_asdf_waveform_fetch
   pyatoa.tests.test_gatherer.test_fetch_event_by_dir
   pyatoa.tests.test_gatherer.test_fetch_inv_by_dir
   pyatoa.tests.test_gatherer.test_fetch_observed_by_dir
   pyatoa.tests.test_gatherer.test_fetch_synthetic_by_dir
   pyatoa.tests.test_gatherer.test_gather_event



.. py:function:: code()

   Example NZ station code


.. py:function:: event_id()

   Example NZ event identifier


.. py:function:: dataset_fid()

   The name for the dataset file id used for temporary storage
   :return:


.. py:function:: cat()

   ObsPy Event Catalog for New Zealand based event with
   GeoNet Event ID: 2018p130600


.. py:function:: event(cat)

   Event from Catalog


.. py:function:: origintime(event)

   The origin time of the example event
   :return:


.. py:function:: config(event_id)

   Default Pyatoa Config object


.. py:function:: gatherer(config, origintime)

   The Gatherer which is responsible for gathering data.


.. py:function:: test_asdf_event_fetch(gatherer, dataset_fid)

   Get event from an ASDFDataSet.


.. py:function:: test_asdf_station_fetch(gatherer, dataset_fid, code)

   Get station from an ASDFDataSet


.. py:function:: test_asdf_waveform_fetch(gatherer, dataset_fid, code, config)

   Get waveforms from an ASDFDataSet


.. py:function:: test_fetch_event_by_dir(gatherer, event_id)

   Get event information based on given directory structure. Test the various
   types of input sources that are allowable by Pyatoa


.. py:function:: test_fetch_inv_by_dir(gatherer, code)

   Get response based on given directory structure


.. py:function:: test_fetch_observed_by_dir(gatherer, code)

   Get waveforms based on given directory strucutre


.. py:function:: test_fetch_synthetic_by_dir(gatherer, code)

   Get synthetics based on given directory strucutre


.. py:function:: test_gather_event(gatherer, dataset_fid)

   Ensure gatherer can get an event from the correct sources


