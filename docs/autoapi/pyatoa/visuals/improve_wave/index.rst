:py:mod:`pyatoa.visuals.improve_wave`
=====================================

.. py:module:: pyatoa.visuals.improve_wave

.. autoapi-nested-parse::

   Functions to create a figure showing progressive waveform improvement over
   the course of a seismic inversion.

   Show the changes in synthetic waveforms with progressive model updates.
   Each individual model gets its on row in the plot.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.visuals.improve_wave.ImproveWave




Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.improve_wave.pairs
   pyatoa.visuals.improve_wave.event_id


.. py:class:: ImproveWave

   A class to plot waveform improvement for a given ASDFDataSet

   .. code:: python

       ds = pyasdf.ASDFDataSet("dataset.h5")
       wi = WaveformImprovement(ds)
       wi.gather("NZ.BFZ", 10, 30)
       wi.plot()
       wi.plot("NZ.KNZ", 8, 30)

   .. py:method:: get_models()

      Figure out which step goes to which iteration to get model numbers


   .. py:method:: gather(sta, min_period, max_period, rotate_to_rtz=False, fix_windows=False, pyflex_preset=False)

      Parse dataset for given station, gather observed and synthetic data,
      preprocess data and return as stream objects.

      :type sta: str
      :param sta: station to gather data for
      :type min_period: float
      :param min_period: minimum filter period in seconds
      :type max_period: float
      :param max_period: maximum filter period in seconds
      :type rotate_to_rtz: bool
      :param rotate_to_rtz: rotate components from NEZ to RTZ
          Config. if False, instrument response will be removed from obs.
      :type fix_windows: bool
      :param fix_windows: dont recalculate windows when gathering
      :type pyflex_preset: str
      :param pyflex_preset: overwrite the pyflex preset provided in the
          Config object


   .. py:method:: gather_simple(event, sta, min_period, max_period, path_dict=None, component=None)

      Manually set the model values based on inspection of the Inspector
      Don't return windows or anything, keep it simple


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


   .. py:method:: plot(sta=None, event_id=None, min_period=None, max_period=None, plot_windows=False, trace_length=None, show=True, save=False, **kwargs)

      Plot waveforms iterative based on model updates

      :type sta: str
      :param sta: station to gather data for, if None, skips gathering
          assuming data has already been gathered
      :type min_period: float
      :param min_period: minimum filter period for waveforms
      :type max_period: float
      :param max_period: maximum filter period for waveforms
      :type plot_windows: bool
      :param plot_windows: plot misfit windows above waveforms
      :type trace_length: list of floats
      :param trace_length: [trace_start, trace_end] will be used to set the x
          limit on the waveform data. If none, no xlim will be set
      :type show: bool
      :param show: Show the plot or do not
      :type save: str
      :param save: if given, save the figure to this path


   .. py:method:: gather_simple(models, event_id, sta, component, min_period, max_period)

      Gather waveforms from manually input model values, usually determined
      by using the Inspector class

      :type models: dict of tuples
      :param models: model values as keys, (iter/step, tag) as tuple value.
          Tags allow multiple datasets to be used, e.g. if an inversion
          spans over multiple legs and more than one dataset is used to
          store waveform data
      :type event_id: str
      :param event_id: name of the event, used to access the ASDFDataSet
      ;type sta: str
      :param sta: station id to gather data for
      :type min_period: float
      :param min_period: period to filter data at
      :type max_period: float
      :param max_period: period to filter data at



.. py:data:: pairs
   :annotation: = [['2013p617227', 'NZ.TOZ', 'Z']]

   

.. py:data:: event_id
   :annotation: = 2013p507880

   

