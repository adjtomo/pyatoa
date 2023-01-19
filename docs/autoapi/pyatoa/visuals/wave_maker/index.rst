:py:mod:`pyatoa.visuals.wave_maker`
===================================

.. py:module:: pyatoa.visuals.wave_maker

.. autoapi-nested-parse::

   Waveform plotting functionality

   Produces waveform plots from stream objects and plots misfit windows
   outputted by Pyflex as well as adjoint sources from Pyadjoint.
   Flexible to allow for only waveform plots, or for the addition of objects
   based on inputs.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.visuals.wave_maker.WaveMaker



Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.wave_maker.align_yaxes
   pyatoa.visuals.wave_maker.pretty_grids
   pyatoa.visuals.wave_maker.format_axis



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.wave_maker.rejected_window_colors


.. py:data:: rejected_window_colors
   

   

.. py:class:: WaveMaker(mgmt, **kwargs)

   Standardized waveform figures featuring observed and synthetic traces,
   STA/LTA waveforms, misfit windows, rejected windows, adjoint sources and
   auxiliary information collected within the workflow.

   WAVEFORM KEYWORD ARGUMENTS:
       FIGURE:
       figsize (tuple): size of the figure, defaults (800, 200) pixels
           per channel
       dpi (float): dots per inch of the figure, defaults 100
       legend_fontsize (int): size
       linewidth (float): line width for all lines on plot, default 1.6
       axes_linewidth (float): line width for axis spines, default 1
       xlim_s (list): time-axis bounds of the plot in seconds, def full trace
       figure (pyplot.Figure): an existing figure object to plot to, rather
           than generating a new figure
       subplot_spec (gridspec.GridSpec): an overlying grid that waveforms
           will be plotted into. Useful for combining waveform plots

       FONTSIZE:
       fontsize (int): font size of the title, axis labels, def 8
       axes_fontsize (int): font size of the tick labels, def 8
       rejected_window_fontsize (int): fontsize for the annotations that
           describe the rejected windows, default 6
       window_anno_fontsize (str): fontsize for window annotation, def 7

       COLORS:
       obs_color (str): color for observed waveform, defaults 'k'
       syn_color (str): color for synthetic waveform, defaults 'r
       stalta_color (str): color of stalta waveform, default 'gray'
       window_color (str): color for misfit windows, default 'orange'
       adj_src_color (str): color for adjoint sources, default 'g'

       ADJOINT SOURCE
       adj_src_linestyle (str, tuple): adjoint souce style, default tight dash
       adj_src_alpha (float): opacity of adjoint source, default 0.4

       WINDOW ANNOTATIONS:
       window_anno (str): a custom string which can contain the optional
           format arguemnts: [max_cc, cc_shift, dlnA, left, length]. None,
           defaults to formatting all arguments
       window_anno_alternate (str): custom string for all windows that
           aren't the first window, useful for dropping the labels for
           parameters, allows for cleaner annotations without
           compromising readability
       window_anno_height (float): annotation height, percentage of y axis,
           default 0.7
       alternate_anno_height (float): optional, shift the annotation height
           each window to prevent overlapping annotations
       window_anno_rotation (float): rotation of annotation (deg), def 0
       window_anno_fontcolor (str): color of annotation text, def 'k'
       window_anno_fontweight (str): weight of font, default 'normal'
       window_anno_bbox (dict): bbox dict for window annotations, None means
           no bounding box

       TOGGLES:
       plot_xaxis (bool): toggle the labels and ticks on the x axis, def True
       plot_yaxis (bool): toggle the labels and ticks on the y axis, def True
       plot_windows (bool): toggle window plotting, default True
       plot_rejected_windows (bool): toggle rejected window plot, default T
       plot_window_annos (bool): toggle window annotations, default True
       plot_staltas (bool): toggle stalta plotting, default True
       plot_adjsrcs (bool): toggle adjoint source plotting, default True
       plot_waterlevel (bool): toggle stalta waterlevel plotting, def True
       plot_arrivals (bool): toggle phase arrival plotting, default True
       plot_legend (bool): toggle legend, default True

       MISC:
       normalize (bool): normalize waveform data before plotting
       set_title (bool or str): create a default title using workflow
           parameters, if str given, overwrites all title
       append_title(str): User appended string to the end of the title.
           useful to get extra information on top of the default title


   .. py:method:: setup_plot(dpi, figsize, twax_off=False)

      Dynamically set up plots according to number_of given

      Calculate the figure size based on DPI, (800, 250) pixels per channel

      :type dpi: float
      :param dpi: dots per inch, to be set by plot()
      :type figsize: tuple
      :param figsize: size of the figure, set by plot()
      :type twax_off: bool
      :param twax_off: if True, dont instantiate a twin-x axis
      :rtype (tw)axes: matplotlib axes
      :return (tw)axes: axis objects


   .. py:method:: plot_waveforms(ax, obs, syn, normalize=False)

      Plot observed and synthetic data on the same axis, label them according
      to their trace ID

      :type ax: matplotlib.axes.Axes
      :param ax: axis object on which to plot
      :type obs: obspy.core.trace.Trace
      :param obs: observed waveform data to plot
      :type syn: obspy.core.trace.Trace
      :param syn: synthetic waveform data to plot
      :type normalize: bool
      :param normalize: option to normalize the data traces between [-1, 1],
          defaults to False, do not normalize


   .. py:method:: plot_stalta(ax, stalta, plot_waterlevel=True)

      Plot the Short-term-average/long-term-average waveform to help visually
      identify the peaks/troughs used to determine windows

      :type ax: matplotlib.axes.Axes
      :param ax: axis object on which to plot
      :type stalta: numpy.ndarray
      :param stalta: data array containing the sta/lta waveform
      :type plot_waterlevel: bool
      :param plot_waterlevel: plot a horizontal line showing the relative
          waterlevel of the sta/lta which is used in determining windows
      :rtype: list of matplotlib.lines.Lines2D objects
      :return: list containing the lines plotted on the axis


   .. py:method:: plot_windows(ax, windows, plot_window_annos=True, plot_phase_arrivals=True)

      Plot misfit windows, add annotations to each window related to
      information contained in the Window object.

      .. note::
          The keyword argument 'window_anno_height' should be given as a
          percentage of visible y-axis, e.g. 0.25 means 25% of the y-axis

      :type ax: matplotlib.axes.Axes
      :param ax: axis object on which to plot
      :type windows: list of pyflex.Window objects
      :param windows: list of windows to plot
      :type plot_window_annos: bool
      :param plot_window_annos: annotate window information onto windows
      :type plot_phase_arrivals: bool
      :param plot_phase_arrivals: make small tick mark if P or S phase arrival
          within the window


   .. py:method:: plot_rejected_windows(ax, rejwin, windows=None, skip_tags=None)

      Plot rejected windows as transparent lines at the bottom of the axis.
      Hardcoded color dictionary (defined at top) used as a way to visually
      identify why certain windows were rejected

      The function performs some array manipulation to exclude rejected
      windows that fall within already chosen time windows to avoid redundant
      plotting.

      :type ax: matplotlib.axes.Axes
      :param ax: axis object on which to plot
      :type rejwin: list of pyflex.Window objects
      :param rejwin: list of rejected windows to plot
      :type windows: list of pyflex.Window objects
      :param windows: list of windows to use for exclusion
      :type skip_tags: list of str
      :param skip_tags: an optional list of tags that can be used to skip
          specific rejected window tags


   .. py:method:: plot_adjsrcs(ax, adjsrc)

      Plot adjoint sources behind streams, time reverse the adjoint source
      that is provided by Pyadjoint so that it lines up with waveforms
      and windows.

      .. note::

         The unit of adjoint source is based on Eq. 57 of Tromp et al. 2005,
         which is a traveltime adjoint source, and includes the units of the
         volumentric delta function since it's assumed this is happening in a
         3D volume.

      :type ax: matplotlib.axes.Axes
      :param ax: axis object on which to plot
      :type adjsrc: pyadjoint.adjoint_source.AdjointSource objects
      :param adjsrc: adjsrc object containing data to plot
      :rtype: list of matplotlib.lines.Lines2D objects
      :return: list containing the lines plotted on the axis


   .. py:method:: plot_amplitude_threshold(ax, obs)

      Plot a line to show the amplitude threshold criteria used by Pyatoa

      :type ax: matplotlib.axes.Axes
      :param ax: axis object on which to plot
      :type obs: obspy.core.trace.Trace
      :param obs: observed trace plotted on the current axis, used to
          determine the peak amplitude value


   .. py:method:: create_title(normalized=False, append_title=None)

      Create the title based on information provided to the class

      :type normalized: bool
      :param normalized: whether or not the data was normalized
      :type append_title: str
      :param append_title: append extra information to title
      :rtype: str
      :return: title string composed of configuration parameters and source
          receiver information


   .. py:method:: plot(show=True, save=False, **kwargs)

      High level plotting function that plots all parts of the class and
      formats the axes nicely



.. py:function:: align_yaxes(ax1, ax2)

   Plotting tool to adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1

   :type ax1: matplotlib axis
   :param ax1: axes to adjust
   :type ax2: matplotlib axis
   :param ax2: axes to adjust


.. py:function:: pretty_grids(input_ax, twax=False, grid=False, fontsize=8, linewidth=1, sci_format=True)

   Standard plot skeleton formatting, thick lines and internal tick marks etc.

   :type input_ax: matplotlib axis
   :param input_ax: axis to prettify
   :type twax: bool
   :param twax: If twax (twin axis), do not set grids
   :type grid: bool
   :param grid: turn on grids of the axes, default grids off
   :type fontsize: float
   :param fontsize: fontsize of the axis tick labels
   :type linewidth: float
   :param linewidth: line width of the axis spines or boundign box
   :type sci_format: bool
   :param sci_format: turn on/off scientific formatting of tick labels
       default scientific format on/True.


.. py:function:: format_axis(input_ax, percent_over=0.125)

   Ensure that the upper and lower y-bounds are the same value

   :type input_ax: matplotlib axis
   :param input_ax: axis to prettify
   :type percent_over: float
   :param percent_over: the percentage above the peak value that the bounds
       of the axis will be set


