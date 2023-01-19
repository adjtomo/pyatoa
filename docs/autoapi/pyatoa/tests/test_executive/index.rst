:py:mod:`pyatoa.tests.test_executive`
=====================================

.. py:module:: pyatoa.tests.test_executive

.. autoapi-nested-parse::

   Test the functionalities of the Executive class which runs many Managers



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_executive.events
   pyatoa.tests.test_executive.stations
   pyatoa.tests.test_executive.config
   pyatoa.tests.test_executive.test_executive_single_event_single_station_no_concurrent
   pyatoa.tests.test_executive.test_executive_single_event_single_station
   pyatoa.tests.test_executive.test_executive_single_event_multi_station



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_executive.propogate


.. py:data:: propogate
   :annotation: = False

   

.. py:function:: events()

   A list of GeoNet event ids used for gathering metadata from GeoNet client


.. py:function:: stations()

   A list of GeoNet station codes used for gathering metadata and waveforms


.. py:function:: config(events)

   A preset Config object that specifies where to grab data from, which
   already exists in the test data directory


.. py:function:: test_executive_single_event_single_station_no_concurrent(tmpdir, config, events, stations)

   Attempt a single event single station processing without using concurrency


.. py:function:: test_executive_single_event_single_station(tmpdir, config, events, stations)

   Attempt a single event single station processing with concurrency


.. py:function:: test_executive_single_event_multi_station(tmpdir, config, events, stations)

   Attempt a single event multi station processing with concurrency


