"""
Wraps the functionalities of WaveMaker and MapMaker into a single class
that can be used to plot waveforms, maps or both together
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa.visuals.map_maker import MapMaker
from pyatoa.visuals.wave_maker import WaveMaker


class ManagerPlotter:
    """
    A simple class used to generate plots for the Manager class. The numerous
    keyword arguments that are accessible to change the style of plotting
    are given below:

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

    MAPPING KEYWORD ARUGMENTS:
        write me.

    """
    def __init__(self, mgmt):
        """
        Manager plotter only requires reference to the Manager object

        :type mgmt: pyatoa.Manager
        :param mgmt: Manager object to be plotted
        """
        self.mgmt = mgmt

    def plot_wav(self, save=None, show=True, **kwargs):
        """
        Plot observed and synthetics waveforms, misfit windows, STA/LTA and
        adjoint sources for all available components. Append information
        about misfit, windows and window selection.

        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        """
        # Call on window making function to produce waveform plots
        wm = WaveMaker(mgmt=self.mgmt, **kwargs)
        wm.plot(show=show, save=save)
        plt.close()

    def plot_map(self, corners=None, save=None, show=True, **kwargs):
        """
        Generate a basemap for a given target region. Plot station and receiver
        from internal station and event attributes.

        :type corners: dict
        :param corners: {lat_min, lat_max, lon_min, lon_max}
            corners to cut the map to, otherwise a global map is provided
        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        """
        # Call external function to generate map
        mm = MapMaker(inv=self.mgmt.inv, cat=self.mgmt.event, **kwargs)
        mm.plot(corners=corners, show=show, save=save)
        plt.close()

    def plot(self, corners=None, dpi=100, figsize=None, show=True, save=False,
             **kwargs):
        """
        Plot the Waveform and Map figures together
        """
        # Default figure size
        if figsize is None:
            figsize = (1400 / dpi, 600 / dpi)

        # Create an overlying GridSpec that will contain both plots
        gs = mpl.gridspec.GridSpec(1, 2, wspace=0.25, hspace=0.)
        fig = plt.figure(figsize=figsize, dpi=dpi)

        # Plot the waveform on the left
        wm = WaveMaker(mgmt=self.mgmt)
        wm.plot(figure=fig, subplot_spec=gs[0], show=False, save=False,
                **kwargs)

        # Plot the map on the right
        mm = MapMaker(inv=self.mgmt.inv, cat=self.mgmt.event, **kwargs)
        ax = fig.add_subplot(gs[1])
        mm.plot(corners=corners, figure=fig, ax=ax, show=False, save=False)

        if save:
            plt.savefig(save)
        if show:
            plt.show()
        else:
            plt.close()

