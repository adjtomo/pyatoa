:py:mod:`pyatoa.visuals.comp_wave`
==================================

.. py:module:: pyatoa.visuals.comp_wave

.. autoapi-nested-parse::

   A stripped down (arguably better) version of the Waveform Improvement class,
   used to simply compare two synthetic waveforms from a given PyASDF DataSet.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.visuals.comp_wave.CompWave



Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.comp_wave.main
   pyatoa.visuals.comp_wave.main_plot_specific



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.comp_wave.a


.. py:class:: CompWave(dsfid, station, min_period, max_period, dsfid_final=None)

   A class to plot waveform improvement between two models for a given dataset

   .. py:method:: _gather_model_from_dataset(dsfid, model=None, init_or_final=None, save_windows=False)

      Gather data from an ASDFDataSet based on the given model (iter/step)

      :type dsfid: str
      :param dsfid: file identifier for the dataset
      :type model: str
      :param model: iteration/step, e.g. 'i01/s00'
      :type init_or_final: str
      :param init_or_final: for choosing default values if model is None
          * 'init': choose the first iteration/step for the initial model
          * 'final': choose final iteration/step for final model
      :return:


   .. py:method:: gather(m_init=None, m_final=None, save_windows=False)

      Gather data from the correct dataset. If no m_init or m_final given,
      will gather the first and last models

      :type m_init: str
      :param m_init: initial iteration/step, e.g. 'i01/s00'
      :type m_final: str
      :param m_final: final iteration/step


   .. py:method:: calculate_vrl(init_or_final)

      Caclulate the logarithmic variance reduction to look at how waveforms
      imrpove from m_init to m_final. Following Eq. 8 of Tape et al. (2010).


   .. py:method:: setup_plot(nrows, ncols, **kwargs)

      Dynamically set up plots according to number_of given
      Returns a list of lists of axes objects
      e.g. axes[i][j] gives the ith column and the jth row

      :type nrows: int
      :param nrows: number of rows in the gridspec
      :type ncols: int
      :param ncols: number of columns in the gridspec
      :rtype axes: matplotlib axes
      :return axes: axis objects


   .. py:method:: _xlim_from_envelope(data, dt)

      Get rough bounds for the xlimits by looking at waveform envelopes
      :return:


   .. py:method:: plot(component_list=None, show=True, save=False, **kwargs)

      Plot waveforms iterative based on model updates

      :type show: bool
      :param show: Show the plot or do not
      :type save: str
      :param save: if given, save the figure to this path


   .. py:method:: plot_with_map(corners=None, dpi=100, figsize=None, show=True, save=False, **kwargs)

      Similar to Manager plotter, plot the waveform comparisons next to a
      source receiver map. Wraps the internal plotting functionality with
      a gridspec



.. py:function:: main(event_id=None, station=None, component=None, xmin=None, xmax=None, cfg='plot')

   Main call script to choose event and station based on what's available


.. py:function:: main_plot_specific()

   


.. py:data:: a
   

   

