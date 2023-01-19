:py:mod:`pyatoa.visuals.insp_plot`
==================================

.. py:module:: pyatoa.visuals.insp_plot

.. autoapi-nested-parse::

   The plotting functionality of the Inspector class. Used to generate statistics
   and basemap like plots from the Inspector DataFrame objects.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.visuals.insp_plot.InspectorPlotter



Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.insp_plot.default_axes
   pyatoa.visuals.insp_plot.colormap_colorbar
   pyatoa.visuals.insp_plot.hover_on_plot
   pyatoa.visuals.insp_plot.get_histogram_stats
   pyatoa.visuals.insp_plot.annotate_txt



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.insp_plot.common_labels


.. py:data:: common_labels
   

   

.. py:class:: InspectorPlotter

   A class of methods for plotting statistics from an Inspector.
   Should not be called on its own, these functions will be inherited by
   the Inspector class automatically.

   .. py:method:: map(event=None, network=None, station=None, show=True, save=False, **kwargs)

      Plot source and receiver locations with map view. Optional arguments
      for only plotting certain stations or events.

      :type event: str or list
      :param event: particular event or list of events to plot
      :type network: str or list
      :param network: particular network or list of networks to plot
      :type station: str or list
      :param station: particular station or list of stations to plot
      :type show: bool
      :param show: Show the plot
      :type save: str
      :param save: fid to save the given figure


   .. py:method:: scatter(x, y, iteration=None, step_count=None, save=None, show=True, **kwargs)

      Create a scatter plot between two chosen keys in the windows attribute

      :type x: str
      :param x: key to choose for the x axis of the plot
      :type y: str
      :param y: key to chooose for the y axis of the plot
      :type iteration: str
      :param iteration: the chosen iteration to plot for, if None will default
          to the latest iteration available
      :type step_coutn: str
      :param step_count: chosen step count. If None, defaults to latest


   .. py:method:: travel_times(iteration=None, step_count=None, component=None, constants=None, t_offset=0, hist=False, hist_max=None, plot_end=False, save=None, show=True, **kwargs)

      Plot relative window starttime (proxy for phase arrival) against
      source-receiver distance, to try to convey which phases are included
      in the measurement.

      Similar to Figure 4.18 in Shearer's Intro to Seismology.

      :type iteration: str
      :param iteration: the chosen iteration to plot for, if None will default
          to the latest iteration available
      :type step_count: str
      :param step_count: chosen step count. If None, defaults to latest
      :type component: str
      :param component: optional specify a measurement component to isolate
          only e.g., 'Z' components to look at Rayleigh waves
      :type constants: list of floats
      :param constants: plot lines of constant velocity to estimate the
          average wavespeed that leads to some of the linear trends
      :type t_offset: float
      :param t_offset: if the synthetic offset time in SPECFEM is set then
          the constant lines will need to be offset by the same amount to
          match the measurements.
      :type hist: bool
      :param hist: create a histogram binning the approximate seismic
          velocities
      :type plot_end: bool
      :param plot_end: if True, plots the beginning and end of the misfit
          window as a vertical line. If False, plots only the beginning of
          the misfit window


   .. py:method:: event_depths(xaxis='longitude', show=True, save=None, **kwargs)

      Create a scatter plot of events at depth. Compresses all events onto a
      single slice, optional choice of showing the x-axis or the y-axis

      :type xaxis: str
      :param xaxis: variable to use as the x-axis on the plot
          'latitude' or 'longitude'
      :type show: bool
      :param show: show the plot
      :type save: str
      :param save: fid to save the figure


   .. py:method:: raypaths(iteration=None, step_count=None, color_by=None, show=True, save=False, vmin=None, vmax=None, **kwargs)

      Plot rays connecting sources and receivers based on the availability
      of measurements. Useful for getting an approximation of resolution.

      :type iteration: int
      :param iteration: iteration to retrieve data from
      :type step_count: int
      :param step_count: step count to retrieve data from
      :type color_by: str
      :param color_by: allow rays to be colored based on a normalized value.
          nwin: color rays by the number of windows available for that path
          misfit: color rays by total misfit
      :type show: bool
      :param show: show the plot
      :type save: str
      :param save: fid to save the figure


   .. py:method:: raypath_density(iteration=None, step_count=None, point_spacing_km=0.5, bin_spacing_km=8, cmap='viridis', show=True, save=False, **kwargs)

      Create a raypath density plot to provide a more deatiled illustration of
      raypath gradients, which may be interpreted alongside tomographic
      inversion results as a preliminary resolution test.

      The idea behind this is to partition each individual raypath line into
      discrete points and then create a 2D histogram with all points

      :type point_spacing_km: float
      :param point_spacing_km: approximate discretization interval for each
          raypath line. Smaller numbers will lead to higher resolution but
          also longer computation time.
      :type bin_spacing_km: float
      :param bin_spacing_km: the bin size in km of the 2d histogram. If
          the same as 'point_spacing_km' then you'll probably just see the
          lines. Should be larger than 'point_spacing_km' for a more
          contour plot looking feel.


   .. py:method:: event_hist(choice, show=True, save=None)

      Make a histogram of event information
      :return:


   .. py:method:: measurement_hist(iteration=None, step_count=None, choice='event', show=True, save=False)

      Make histograms of measurements for stations or events to show the
      distribution of measurements.

      :type iteration: str
      :param iteration: iteration number e.g. 'i00'
      :type step_count: str
      :param step_count: step count e.g. 's00'
      :type choice: str
      :param choice: choice of making hist by 'event' or 'station'
      :type show: bool
      :param show: Show the plot
      :type save: str
      :param save: fid to save the given figure


   .. py:method:: station_event_misfit_map(station, iteration, step_count, choice, show=True, save=False, **kwargs)

      Plot a single station and all events that it has measurements for.
      Events will be colored by choice of value: misfit or nwin (num windows)

      :type station: str
      :param station: specific station to use for map
      :type iteration: str
      :param iteration: iteration number e.g. 'i00'
      :type step_count: str
      :param step_count: step count e.g. 's00'
      :type choice: str
      :param choice: choice of misfit value, either 'misfit' or 'nwin'
      :type show: bool
      :param show: Show the plot
      :type save: str
      :param save: fid to save the given figure


   .. py:method:: event_station_misfit_map(event, iteration, step_count, choice, show=True, save=False, **kwargs)

      Plot a single event and all stations with measurements. Stations are
      colored by choice of value: misfit or nwin (number of windows)

      :type event: str
      :param event: specific event to use for map
      :type iteration: str
      :param iteration: iteration number e.g. 'i00'
      :type step_count: str
      :param step_count: step count e.g. 's00'
      :type choice: str
      :param choice: choice of misfit value, either 'misfit' or 'nwin'
      :type show: bool
      :param show: Show the plot
      :type save: str
      :param save: fid to save the given figure


   .. py:method:: event_misfit_map(choice=None, iteration=None, step_count=None, show=True, save=False, **kwargs)

      Plot all events on a map and their corresponding scaled misfit value

      :type iteration: str
      :param iteration: iteration number e.g. 'i00'
      :type step_count: str
      :param step_count: step count e.g. 's00'
      :type choice: str
      :param choice: choice of misfit value, either 'misfit' or 'nwin' or
          'unscaled_misfit'
      :type show: bool
      :param show: Show the plot
      :type save: str
      :param save: fid to save the given figure


   .. py:method:: hist(iteration=None, step_count=None, iteration_comp=None, step_count_comp=None, f=None, ax=None, event=None, station=None, choice='cc_shift_in_seconds', binsize=None, show=True, save=None, **kwargs)

      Create a histogram of misfit information for either time shift or
      amplitude differences. Option to compare against different iterations,
      and to look at different choices.

      Choices are any column value in the Inspector.windows attribute

      :type iteration: str
      :param iteration: iteration to choose for misfit
      :type step_count: str
      :param step_count: step count to query, e.g. 's00'
      :type iteration_comp: str
      :param iteration_comp: iteration to compare with, will be plotted in
          front of `iteration`
      :type step_count_comp: str
      :param step_count_comp: step to compare with
      :type f: matplotlib.figure
      :param f: plot to an existing figure
      :type ax: matplotlib.axes._subplots.AxesSubplot
      :param ax: plot to an existing axis e.g. to string together histograms
      :type event: str
      :param event: filter for measurements for a given event
      :type station: str
      :param station: filter for measurements for a given station
      :type choice: str
      :param choice: choice of 'cc_shift_s' for time shift, or 'dlnA' as
          amplitude
      :type binsize: float
      :param binsize: size of the histogram bins
      :type show: bool
      :param show: show the plot
      :type save: str
      :param save: fid to save the figure


   .. py:method:: plot_windows(iteration=None, step_count=None, iteration_comp=None, step_count_comp=None, choice='cc_shift_in_seconds', event=None, network=None, station=None, component=None, no_overlap=True, distances=False, annotate=False, bounds=False, show=True, save=False, **kwargs)

      Show lengths of windows chosen based on source-receiver distance, akin
      to Tape's Thesis or to the LASIF plots. These are useful for showing
      which phases are chosen, and window choosing behavior as distance
      increases and (probably) misfit increases.

      :type iteration: str
      :param iteration: iteration to analyze
      :type step_count: str
      :param step_count: step count to query, e.g. 's00'
      :type iteration_comp: str
      :param iteration_comp: Optional, if provided, difference the 'choice'
          values with the chosen 'iteration/step'. Useful for easily checking
          for improvement. Only works if the windows are the same.
      :type step_count_comp: str
      :param step_count_comp: associated step count for 'iteration_comp'
      :type event: str
      :param event: filter for measurements for a given event
      :type network: str
      :param network: filter for measurements for a given network
      :type station: str
      :param station: filter for measurements for a given station
      :type component: str
      :param component: choose a specific component to analyze
      :type choice: str
      :param choice: choice of value to define the colorscale by. These relate
          to the keys of Inspector.windows. Default is 'cc_shift_in_seconds'
      :type no_overlap: bool
      :param no_overlap: If real distances are used, many src-rcv pairs are
          at the same or very similar distances, leading to overlapping
          rectangles. If this is set to True, to minimize overlap, the
          function will try to shift the distance to a value that hasn't yet
          been plotted. It will alternate larger positive and negative values
          until something is found. Will lead to non-real distances.
      :type distances: bool
      :param distances: If set False, just plot one window atop the other,
          which makes for more concise, easier to view plots, but
          then real distance information is lost, only relative distance
          kept.
      :type annotate: bool
      :param annotate: If True, will annotate event and station information
          for each window. May get messy if `distances == True` and
          `no_overlap == False` because you will get many overlapping
          annotations. Works ideally if `distances == False`.
      :type bounds: bool or list of float
      :param bounds:
          * (bool) False: set default bounds based on the min and max of data
          * (bool) True: set default bounds equal, based on abs max of data
          * (list) Manually set the bounds of the colorbar
      :type show: bool
      :param show: show the plot after generating
      :type save: str
      :param save: save the plot to the given filename

      :keyword float alpha: The opacity of the rectangles, defaults to 0.25
      :keyword str cmap: The colormap used to plot the values of `choice`
      :keyword str cbar_label: The label for the colorbar
      :keyword float rectangle_height: The vertical size of the rectangles, defaults to 1.
      :keyword float anno_shift: The distance in seconds to shift the plot to accomodate
                                 annotations. This needs to be played as its based on the length
                                 of the strings that are used in the annotations.


   .. py:method:: convergence(windows='length_s', trials=False, show=True, save=None, normalize=False, float_precision=3, annotate=False, restarts='default', restart_annos=None, xvalues='model', **kwargs)

      TO DO:
      Separate the sorting functionality from the plotting functionality,
      this function is too confusing.

      Plot the convergence rate over the course of an inversion.
      Scatter plot of total misfit against iteration number, or by step count

      .. note::
          Because misfits are floats, they wont be exactly equal, so we need
          to set some small tolerance in which they can differ

      :type windows: str or bool
      :param windows: parameter to use for Inspector.measurements() to
          determine how to illustrate measurement number, either by:

          * length_s: cumulative window length in seconds
          * nwin: number of misfit windows
          * None: will not plot window information
      :type trials: str
      :param trials: plot the discarded trial step function evaluations from
          the line searches. Useful for understanding optimization efficiency

          * marker: plot trial steps as red x's at their respective misfit val
          * text: annotate the number of trial steps but not their misfit val
      :type normalize: bool
      :param normalize: normalize the objective function values between [0, 1]
      :type float_precision: int
      :param float_precision: acceptable floating point precision for
          comparisons of misfits. Defaults to 3 values after decimal
      :type restarts: list of int
      :param restarts: If the inversion was restarted, e.g. for parameter
          changes, then the convergence figure should separate two line plots.
          This list allows the User to tell the function where to separate
          the convergence plot. The integers should correspond to indices of
          the Inspector.models attribute.
      :type annotate: bool
      :param annotate: annotate misfit values next to markers
      :type restart_annos: list of str
      :param restart_annos: if restarts is not None, allow annotating text
          next to each restart. Useful for annotating e.g. parameter changes
          that accompany each restart
      :type xvalues: str
      :param xvalues: How the x-axis should be labelled, available:

          * model: plot the model number under each point
          * eval: number sequentially from 1
      :type show: bool
      :param show: show the plot after making it
      :type save: str
      :param save: file id to save the figure to



.. py:function:: default_axes(ax, cbar=None, **kwargs)

   Ensure that all plots have the same default look. Should be more flexible
   than setting rcParams or having a style sheet. Also allows the same kwargs
   to be thrown by all functions so that the function calls have the same
   format.



.. py:function:: colormap_colorbar(cmap, vmin=0.0, vmax=1.0, dv=None, cbar_label='', extend='neither')

   Create a custom colormap and colorbar

   :type cmap: matplotlib.colors.ListedColormap
   :param cmap: colormap to use, called like plt.cm.viridis
   :type vmin: float
   :param vmin: min value for colormap
   :type vmax: float
   :param vmax: max value for colormap
   :type dv: float
   :param dv: colormap boundary separations, if None, continuous colorbar
   :type cbar_label: str
   :param cbar_label: label for colorbar
   :rtype:
   :return:


.. py:function:: hover_on_plot(f, ax, obj, values, dissapear=True)

   Allow for hover on a plot for custom annotated information

   .. note::
       This functionality is copied from StackOverflow:
       https://stackoverflow.com/questions/7908636/possible-to-make-labels-appear-when-hovering-over-a-point-in-matplotlib

   :type f: matplotlib.figure.Figure
   :param f: figure object for hover
   :type ax: matplotlib.axes._subplot.AxesSubplot
   :param ax: axis object for hover
   :type obj: matplotlib.collections.PathCollection or
               matplotlib.lines.Line2D
   :param obj: scatter plot, returned from plt.scatter() or plt.plot()
   :type values: list of str
   :param values: list of annotations
   :type dissapear: bool
   :param dissapear: annotations dissapear when mouse moves off
   :rtype hover: function
   :return hover: the hover function to be passed to matplotlib


.. py:function:: get_histogram_stats(n, bins)

   Get mean, variance and standard deviation from a histogram

   :type n: array or list of arrays
   :param n: values of histogram bins
   :type bins: array
   :param bins: edges of the bins


.. py:function:: annotate_txt(ax, txt, anno_location='lower-right', **kwargs)

   Convenience function to annotate some information

   :type ax: matplot.axes._subplots.AxesSubplot
   :param ax: axis to annotate onto
   :type txt: str
   :param txt: text to annotate
   :type anno_location: str
   :param anno_location: location on the figure to annotate
       available: bottom-right


