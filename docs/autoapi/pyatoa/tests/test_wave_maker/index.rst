:py:mod:`pyatoa.tests.test_wave_maker`
======================================

.. py:module:: pyatoa.tests.test_wave_maker

.. autoapi-nested-parse::

   Test the functionalities of the Pyaflowa Manager class



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_wave_maker.st_obs
   pyatoa.tests.test_wave_maker.st_syn
   pyatoa.tests.test_wave_maker.cat
   pyatoa.tests.test_wave_maker.event
   pyatoa.tests.test_wave_maker.inv
   pyatoa.tests.test_wave_maker.config
   pyatoa.tests.test_wave_maker.mgmt
   pyatoa.tests.test_wave_maker.wm
   pyatoa.tests.test_wave_maker.test_setup_plot
   pyatoa.tests.test_wave_maker.test_plot_waveforms
   pyatoa.tests.test_wave_maker.test_plot_stalta
   pyatoa.tests.test_wave_maker.test_plot_windows
   pyatoa.tests.test_wave_maker.test_plot_adjsrcs
   pyatoa.tests.test_wave_maker.test_plot_rejected_windows
   pyatoa.tests.test_wave_maker.test_plot



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


.. py:function:: mgmt(config, event, st_obs, st_syn, inv)

   A manager filled with data and that has progressed through the workflow


.. py:function:: wm(mgmt)

   WaveMaker object to create waveform figures


.. py:function:: test_setup_plot(wm)

   Make sure the correct number of axes are established


.. py:function:: test_plot_waveforms(wm)

   Plot waveforms by themselves


.. py:function:: test_plot_stalta(wm)

   Plot STALTA waveforms by themselves


.. py:function:: test_plot_windows(wm)

   Test plotting windows


.. py:function:: test_plot_adjsrcs(wm)

   Test plotting adjoint sources


.. py:function:: test_plot_rejected_windows(wm)

   Test plotting rejected windows


.. py:function:: test_plot(wm)

   Test the main plotting function


