"""
Waveform plotting functionality

Produces waveform plots from stream objects and plots misfit windows
outputted by Pyflex as well as adjoint sources from Pyadjoint.
Flexible to allow for only waveform plots, or for the addition of objects
based on inputs.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pyatoa.utils.calculate import normalize_a_to_b, abs_max


# Hardcoded colors that represent rejected misfit windows. Description from 
# pyflex.WindowSelector function that rejected the window
rejected_window_colors = {
    "water_level": "C0",  # reject_on_minima_water_level()
    "prominence": "C1",   # reject_on_prominence_of_central_peak()
    "phase_sep": "C2",   # reject_on_phase_separation()
    "curtail": "C3",   # curtail_length_of_windows()
    "min_length": "C4",  # reject_windows_based_on_minimum_length()
    "cc": "C5",   # reject_based_on_data_fit_criteria()
    "tshift": "C6",   # reject_based_on_data_fit_criteria()
    "dlna": "C7",   # reject_based_on_data_fit_criteria()
    "s2n": "C8",  # reject_based_on_signal_to_noise_ratio()
    "traveltimes": "C9",  # reject_on_traveltimes()
    "amplitude": "C10"  # pyatoa reject_on_global_amplitude_ratio()
    }


class WaveMaker:
    """
    Standardized waveform figures featuring observed and synthetic traces, 
    STA/LTA waveforms, misfit windows, rejected windows, adjoint sources and
    auxiliary information collected within the workflow.
    """
    def __init__(self, mgmt, **kwargs):
        """
        Introduce class-wide objects that accessed for plotting the figure
        Accessible keyword arguments can be found in ManagerPlotter()
        """
        self.st_obs = mgmt.st_obs.copy()
        self.st_syn = mgmt.st_syn.copy()
        self.config = mgmt.config

        # If auxiliary data is None, initialize as empty dictionary so that 
        # waveforms can still be plotted
        self.windows = mgmt.windows or {}
        self.staltas = mgmt.staltas or {}
        self.adjsrcs = mgmt.adjsrcs or {}
        self.rejwins = mgmt.rejwins or {}

        self.fig = None
        self.axes = None
        self.twaxes = None
        self.kwargs = kwargs

        self.time_axis = self.st_obs[0].times() + mgmt.stats.time_offset_sec

    def setup_plot(self, dpi, figsize, twax_off=False):
        """
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
        """
        # Optional kwargs to override the figure and axis instantiation
        figure = self.kwargs.get("figure", None)
        subplot_spec = self.kwargs.get("subplot_spec", None)

        # Axis related kwargs
        fontsize = self.kwargs.get("axes_fontsize", 8)
        axes_linewidth = self.kwargs.get("axes_linewidth", 1)

        # Initiate the figure and fill it up with grids
        if figure is None:
            self.fig = plt.figure(figsize=figsize, dpi=dpi)
        else:
            self.fig = figure

        nrows, ncols = len(self.st_obs), 1
        heights = [1] * len(self.st_obs)
        if subplot_spec is None:
            gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=heights,
                                       hspace=0)
        else:
            # gridspeception!
            gs = mpl.gridspec.GridSpecFromSubplotSpec(nrows, ncols,
                                                      subplot_spec=subplot_spec,
                                                      height_ratios=heights,
                                                      hspace=0)
        axes, twaxes = [], []
        for i in range(gs.get_geometry()[0]):
            if i == 0:
                ax = plt.subplot(gs[i])
            else:
                ax = plt.subplot(gs[i], sharex=axes[0])
            twinax = ax.twinx()

            pretty_grids(twinax, twax=True, fontsize=fontsize, 
                         linewidth=axes_linewidth)
            pretty_grids(ax, fontsize=fontsize, linewidth=axes_linewidth)
            twaxes.append(twinax)
            axes.append(ax)

        # remove x-tick labels except for last axis
        for ax in axes[0:-1]:
            plt.setp(ax.get_xticklabels(), visible=False)

        # option to turn off twin axis
        if twax_off:
            for twax in twaxes:
                twax.axes.get_xaxis().set_visible(False)
                twax.axes.get_yaxis().set_visible(False)

        self.axes = axes
        self.twaxes = twaxes

    def plot_waveforms(self, ax, obs, syn, normalize=False):
        """
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
        """
        obs_color = self.kwargs.get("obs_color", "k")
        syn_color = self.kwargs.get("syn_color", "r")
        linewidth = self.kwargs.get("linewidth", 1.6)

        # Option to normalize the data traces
        if normalize:
            obs = obs.copy()
            syn = syn.copy()
            obs.data /= np.abs(obs.data.max())
            syn.data /= np.abs(syn.data.max())

        # Change the legend label depending on syn-syn or obs-syn case
        if self.config.synthetics_only:
            obs_tag = "TRUE"
        else:
            obs_tag = "OBS"

        # Convention of black for obs, red for syn
        a1, = ax.plot(self.time_axis, obs.data, obs_color, zorder=11, 
                      label=f"{obs.id} ({obs_tag})", linewidth=linewidth)
        a2, = ax.plot(self.time_axis, syn.data, syn_color, zorder=10, 
                      label=f"{syn.id} (SYN)", linewidth=linewidth)

        return [a1, a2]

    def plot_stalta(self, ax, stalta, plot_waterlevel=True):
        """
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
        """
        linewidth = self.kwargs.get("linewidth", 1.6)
        stalta_color = self.kwargs.get("stalta_color", "gray")
        fontsize = self.kwargs.get("stalta_wl_fontsize", 6)

        # Get waterlevel from the Pyflex config object
        stalta_wl = self.config.pyflex_config.stalta_waterlevel

        # Bounds for use in setting positions
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        # Normalize the bounds of the sta/lta to the y-axis because we only
        # care about the phase and waterlevel
        stalta = normalize_a_to_b(stalta, ymin, ymax)
        waterlevel = (ymax - ymin) * stalta_wl + ymin

        b2, = ax.plot(self.time_axis, stalta, stalta_color, alpha=0.4, 
                      linewidth=linewidth, zorder=9, label=f"STA/LTA")

        if plot_waterlevel:
            # Plot the waterlevel of the STA/LTA defined by Pyflex Config
            ax.axhline(y=waterlevel, xmin=self.time_axis[0], 
                       xmax=self.time_axis[-1], alpha=0.4, zorder=8, 
                       linewidth=linewidth, c=stalta_color, linestyle='--')
            ax.annotate(s=f"stalta_waterlevel = {stalta_wl}", alpha=0.7, 
                        fontsize=fontsize,
                        xy=(0.75 * (xmax - xmin) + xmin, waterlevel)
                        )

        return [b2]

    def plot_windows(self, ax, windows, plot_window_annos=True, 
                     plot_phase_arrivals=True):
        """
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
        """
        window_anno = self.kwargs.get("window_anno", None)
        window_color = self.kwargs.get("window_color", "orange")
        window_anno_fontsize = self.kwargs.get("window_anno_fontsize", 6)
        window_anno_height = self.kwargs.get("window_anno_height", 1)
        window_anno_rotation = self.kwargs.get("window_anno_rotation", 90)
        window_anno_fontcolor = self.kwargs.get("window_anno_fontcolor", "k")
        window_anno_fontweight = self.kwargs.get("window_anno_fontweight", 
                                                 "normal")
        window_anno_ha = self.kwargs.get("window_anno_ha", "left")
        window_anno_va = self.kwargs.get("window_anno_va", "top")
        window_anno_bbox = self.kwargs.get("window_anno_box", None)
        
        # Determine heights for the annotations, allow alternating heights so 
        # that adjacent windows don't write over one another
        ymin, ymax = ax.get_ylim()
        y_anno = window_anno_height * (ymax - ymin) + ymin

        # Default window annotation string
        if window_anno is None:
            window_anno = "cc={max_cc:.2f} / dT={cc_shift:.2f} / dA={dlnA:.2f}"

        for j, window in enumerate(windows):
            tleft = window.left * window.dt + self.time_axis[0]
            tright = window.right * window.dt + self.time_axis[0]

            # Misfit windows as rectangle; taken from Pyflex
            ax.add_patch(Rectangle(xy=(tleft, ymin), width=tright - tleft, 
                                   height=(ymax + np.abs(ymin)), 
                                   fc=window_color, ec="k", 
                                   alpha=(window.max_cc_value ** 2) * 0.25,
                                   zorder=10
                                   )
                         )

            # Outline the rectangle with solid black lines
            for x_ in [tleft, tright]:
                ax.axvline(x=x_, ymin=0, ymax=1, color="k", alpha=1., zorder=11)

            if plot_window_annos:
                # Annotate window information into each window
                t_anno = (tright - tleft) * 0.025 + tleft
                s_anno = window_anno.format(
                                i=j+1,
                                max_cc=window.max_cc_value,
                                cc_shift=window.cc_shift * window.dt,
                                dlnA=window.dlnA,
                                left=tleft,
                                length=tright - tleft)
                ax.annotate(s=s_anno, ha=window_anno_ha, va=window_anno_va,
                            xy=(t_anno, y_anno),
                            zorder=12, fontsize=window_anno_fontsize,
                            rotation=window_anno_rotation, 
                            color=window_anno_fontcolor,
                            fontweight=window_anno_fontweight,
                            bbox=window_anno_bbox
                            )

            if plot_phase_arrivals:
                for phase_arrivals in window.phase_arrivals:
                    if phase_arrivals["name"] in ["p", "s"]:
                        ax.axvline(x=phase_arrivals["time"], ymin=0, ymax=0.05,
                                   color='b', alpha=0.5
                                   )
                        ax.annotate(s=phase_arrivals["name"],
                                    xy=(0.975 * phase_arrivals["time"], 
                                        0.05 * (ymax-ymin) + ymin),
                                    fontsize=8
                                    )

    def plot_rejected_windows(self, ax, rejwin, windows=None, skip_tags=None):
        """
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
        """
        fontsize = self.kwargs.get("rejected_window_fontsize", 6)

        # By default skip the water_level tag, from experience those windows
        # usually just cover the entire window search range, i.e. uninformative
        if skip_tags is None:
            skip_tags = ["water_level"]

        ymin, ymax = ax.get_ylim()
        dy = 0.075 * (ymax - ymin)  # increment for bar location

        # Chosen time windows - collapse adjacent windows into single windows by 
        # looking for repeated values and removing them, then recombining array
        if windows is not None:
            win_arr = np.array(
                [[_.left * _.dt, _.right * _.dt] for _ in windows])
            win_arr, count = np.unique(win_arr, return_counts=True)
            idx_vals_repeated = np.where(count > 1)[0]
            win_arr = np.delete(win_arr, idx_vals_repeated)
            win_arr = win_arr.reshape(len(win_arr) // 2, 2)
        else:
            win_arr = None

        for tag in rejwin.keys():
            # Skip plotting certain window rejects
            if tag in skip_tags:
                continue

            # We will compare the start and endtimes using a boolean array
            rwin_arr = np.array([[_.left * _.dt, _.right * _.dt] for _ in
                                 rejwin[tag]]
                                )

            # Check if rejected windows are contained within the window bounds
            if win_arr is not None:
                bool_arr = None
                for wa in win_arr:
                    bool_arr_ = np.logical_and(rwin_arr[:, 0] >= wa[0],  # start
                                               rwin_arr[:, 1] <= wa[1]   # end
                                               )
                    if bool_arr is not None:
                        bool_arr = np.logical_and(bool_arr, bool_arr_)
                    else:
                        bool_arr = bool_arr_
                rwin_arr = rwin_arr[~bool_arr]

            # Negate the booleans to exclude rej windows within bounds
            if rwin_arr.any():
                for rw in rwin_arr:
                    # Shift rejected windows by the proper time offset
                    rw += self.time_axis[0]
                    # Plot as rectangle, shorter windows get larger zorder
                    ax.add_patch(Rectangle(xy=(rw[0], ymin),
                                           width=rw[1] - rw[0],
                                           height=dy, ec="k", alpha=0.25,
                                           fc=rejected_window_colors[tag],
                                           zorder=15 + 1 / (rw[1] - rw[0])
                                           )
                                 )

                # Annotate the leftmost rejected window point with the tag
                ax.annotate(xy=(rwin_arr[:, 0].min(), ymin), 
                            s=tag.replace("_", " "), fontsize=fontsize,
                            zorder=14)
                ymin -= dy

        # Reset ylimits based on the extent of the rejected windows
        ax.set_ylim([-abs(ymin), ymax])

    def plot_adjsrcs(self, ax, adjsrc):
        """
        Plot adjoint sources behind streams, time reverse the adjoint source
        that is provided by Pyadjoint so that it lines up with waveforms 
        and windows.

        Note:
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
        """
        linewidth = self.kwargs.get("linewidth", 1.6)
        linestyle = self.kwargs.get("adj_src_linestyle", (0, (5, 1)))
        color = self.kwargs.get("adj_src_color", "g")
        alpha = self.kwargs.get("adj_src_alpha", 0.4)

        # Time reverse adjoint source; line up with waveforms
        b1, = ax.plot(self.time_axis, adjsrc.adjoint_source[::-1], color, 
                      alpha=alpha, linewidth=linewidth, linestyle=linestyle, 
                      zorder=9,
                      label=fr"Adjoint Source ($\chi$={adjsrc.misfit:.2f})"
                      )
        return [b1]

    def plot_amplitude_threshold(self, ax, obs):
        """
        Plot a line to show the amplitude threshold criteria used by Pyatoa

        :type ax: matplotlib.axes.Axes
        :param ax: axis object on which to plot    
        :type obs: obspy.core.trace.Trace
        :param obs: observed trace plotted on the current axis, used to 
            determine the peak amplitude value
        """            
        xmin, xmax = ax.get_xlim()
        threshold_amp = abs(self.config.win_amp_ratio * abs_max(obs.data))

        # Plot both negative and positive bounds
        for sign in [-1, 1]:
            ax.axhline(y=sign * threshold_amp, xmin=self.time_axis[0], 
                       xmax=self.time_axis[-1], alpha=0.35, zorder=6, 
                       linewidth=1.25, c='k', linestyle=':')

        # Annotate window amplitude ratio
        ax.annotate(s=f"{self.config.win_amp_ratio * 100:.0f}% peak amp. obs.", 
                    alpha=0.7, xy=(0.85 * (xmax-xmin) + xmin, threshold_amp), 
                    fontsize=8)

    def create_title(self, normalized=False, append_title=None):
        """
        Create the title based on information provided to the class

        :type normalized: bool
        :param normalized: whether or not the data was normalized
        :type append_title: str
        :param append_title: append extra information to title
        :rtype: str
        :return: title string composed of configuration parameters and source
            receiver information
        """
        title = f"{self.st_obs[0].stats.network}.{self.st_obs[0].stats.station}"

        # Event id may not be given, if it is, append to title
        if self.config.event_id is not None:
            title += f" {self.config.event_id}"

        # Filter bounds to plot title
        title += f" [{self.config.min_period}s, {self.config.max_period}s]"

        # Tell the User if the data has been normalized
        if normalized:
            title += " (normalized) "

        # Add information about the iteration, windowing and misfit measurement
        title += "\n"
        if self.config.iter_tag is not None:
            title += self.config.iter_tag
        if self.config.step_tag is not None:
            title += self.config.step_tag

        # Add information about the Pyflex and Pyadjoint parameters used
        if self.kwargs.get("plot_stalta", True):
            title += f" pyflex={self.config.pyflex_preset}, "
        if self.kwargs.get("plot_adjsrc", True):
            title += f" pyadjoint={self.config.adj_src_type}, "

        # User appended title information
        if append_title is not None:
            title = " ".join([title, append_title])

        return title

    def plot(self, show=True, save=False, **kwargs):
        """
        High level plotting function that plots all parts of the class and 
        formats the axes nicely
        """
        # Allow additional kwargs to be passed in to the plot argument
        self.kwargs.update(kwargs)

        # Distribute some kwargs before starting
        dpi = self.kwargs.get("dpi", 100)
        figsize = self.kwargs.get("figsize",
                                  (800 / dpi, 200 * len(self.st_obs) / dpi))
        fontsize = self.kwargs.get("fontsize", 8)
        legend_fontsize = self.kwargs.get("legend_fontsize", 6)
        append_title = self.kwargs.get("append_title", None)
        normalize = self.kwargs.get("normalize", False)
        xlim_s = self.kwargs.get("xlim_s", None)
        percent_over = self.kwargs.get("percent_over", 0.125)
        set_title = self.kwargs.get("set_title", True)
        plot_xaxis = self.kwargs.get("plot_xaxis", True)
        plot_yaxis = self.kwargs.get("plot_yaxis", True)

        plot_windows = self.kwargs.get("plot_windows", True)
        plot_rejected_windows = self.kwargs.get("plot_rejected_windows", True)
        plot_window_annos = self.kwargs.get("plot_window_anno", True)
        plot_staltas = self.kwargs.get("plot_stalta", True)
        plot_adjsrcs = self.kwargs.get("plot_adjsrc", True)
        plot_waterlevel = self.kwargs.get("plot_waterlevel", True)
        plot_arrivals = self.kwargs.get("plot_arrivals", True)
        plot_legend = self.kwargs.get("legend", True)

        # If nothing on the twin axis, this will turn off tick marks
        twax_off = bool(not plot_staltas or not plot_adjsrcs)

        # Plot per component in the same fashion, only if observed data exists
        self.setup_plot(dpi=dpi, figsize=figsize, twax_off=twax_off)
        for i, obs in enumerate(self.st_obs):
            comp = obs.stats.channel[-1]

            # Try to retrieve auxiliary data by component name
            syn = self.st_syn.select(component=comp)[0]
            windows = self.windows.get(comp, None)
            adjsrc = self.adjsrcs.get(comp, None)
            stalta = self.staltas.get(comp, None)
            rejwin = self.rejwins.get(comp, None)

            # Begin plotting by distributing axis objects
            ax = self.axes[i]
            twax = self.twaxes[i]
            lines = []  # List of lines for making the legend
            lines += self.plot_waveforms(obs=obs, syn=syn, ax=ax, 
                                         normalize=normalize)

            if rejwin is not None and plot_rejected_windows:
                self.plot_rejected_windows(ax=ax, rejwin=rejwin,
                                           windows=windows)

            # Format now as windows use the y-limits for height of windows
            format_axis(ax, percent_over)

            if windows is not None and plot_windows: 
                self.plot_windows(ax=ax, windows=windows, 
                                  plot_window_annos=plot_window_annos, 
                                  plot_phase_arrivals=plot_arrivals)

            if adjsrc is not None and plot_adjsrcs:
                lines += self.plot_adjsrcs(ax=twax, adjsrc=adjsrc)
                if i == len(self.st_obs) // 2:  
                    # middle trace: append units of the adjoint source on ylabel
                    twax.set_ylabel("adjoint source [m$^{-4}$ s]", rotation=270, 
                                    labelpad=20, fontsize=fontsize)
            else:
                twax.set_yticklabels([])  # turn off yticks if no adjsrc

            # Format twax because stalta will use y-limits for its waveforms
            format_axis(twax)

            if stalta is not None and plot_staltas:
                lines += self.plot_stalta(ax=twax, stalta=stalta,
                                          plot_waterlevel=plot_waterlevel
                                          )
                if self.config.win_amp_ratio > 0:
                    self.plot_amplitude_threshold(ax=ax, obs=obs)

            # Finish with axes formatting
            ax.set_ylabel(comp.upper(), fontsize=fontsize)
            if i == len(self.st_obs) // 2:
                # Middle trace will carry the units of the waveforms
                units = {"DISP": "displacement [m]", 
                         "VEL": "velocity [m/s]",
                         "ACC": "acceleration [m/s^2]",
                         "none": ""}[self.config.unit_output]
                ax.set_ylabel(f"{units}\n{ax.get_ylabel()}", fontsize=fontsize)

            if plot_legend:
                labels = [l.get_label() for l in lines]
                ax.legend(lines, labels, prop={"size": legend_fontsize}, 
                          loc="upper right")

            align_yaxes(ax, twax)

        # Final touch ups for the figure
        if isinstance(set_title, bool):
            if set_title:
                self.axes[0].set_title(self.create_title(
                    append_title=append_title, normalized=normalize),
                    fontsize=fontsize
                )
        else:
            self.axes[0].set_title(set_title, fontsize=fontsize)

        if xlim_s is not None:
            self.axes[0].set_xlim(xlim_s)
        else:
            self.axes[0].set_xlim([self.time_axis[0], self.time_axis[-1]])
        self.axes[-1].set_xlabel("time [s]", fontsize=fontsize)

        # Option to turn off tick labels and axis labels
        if not plot_xaxis or not plot_yaxis:
            for axis in [self.axes, self.twaxes]:
                for ax in axis:
                    if not plot_xaxis:
                        ax.axes.xaxis.set_ticklabels([])
                        ax.set_xlabel("")
                    if not plot_yaxis:
                        ax.axes.yaxis.set_ticklabels([])
                        ax.set_ylabel("")

        if save:
            plt.savefig(save, figsize=figsize, dpi=dpi)
        if show:
            plt.show()


def align_yaxes(ax1, ax2):
    """
    Plotting tool to adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1

    :type ax1: matplotlib axis
    :param ax1: axes to adjust
    :type ax2: matplotlib axis
    :param ax2: axes to adjust
    """
    ymin_a1, ymax_a1 = ax1.get_ylim()
    ymin_a2, ymax_a2 = ax2.get_ylim()

    _, y1 = ax1.transData.transform((0, (ymax_a1+ymin_a1)/2))
    _, y2 = ax2.transData.transform((0, (ymax_a2+ymin_a2)/2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))

    ax2.set_ylim(ymin_a2+dy, ymax_a2+dy)


def pretty_grids(input_ax, twax=False, grid=False, fontsize=8, linewidth=1,
                 sci_format=True):
    """
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
    """
    input_ax.set_axisbelow(True)
    if sci_format:
        input_ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # Make sure the scientific notation also has the same fontsize
        input_ax.yaxis.get_offset_text().set_fontsize(fontsize)
    # Ensure the twin axis doesn't clash with ticks of the main axis
    if twax:
        left = False
        right = True
    else:
        left = True
        right = False

    input_ax.tick_params(which='major', direction='in', top=True, right=right,
                         left=left, length=8, labelsize=fontsize, 
                         width=2*linewidth/3)
    input_ax.tick_params(which='minor', direction='in', top=True, right=right,
                         left=left, length=4, labelsize=fontsize, 
                         width=2*linewidth/3)

    for axis in ["top", "bottom", "left", "right"]:
        input_ax.spines[axis].set_linewidth(linewidth)

    # Set the grids 'on' only if main axis
    if not twax:
        input_ax.minorticks_on()
        if grid:
            for axis_ in ['major', 'minor']:
                input_ax.grid(which=axis_, linestyle=':', linewidth='0.5',
                              color='k', alpha=0.25)


def format_axis(input_ax, percent_over=0.125):
    """
    Ensure that the upper and lower y-bounds are the same value

    :type input_ax: matplotlib axis
    :param input_ax: axis to prettify
    :type percent_over: float
    :param percent_over: the percentage above the peak value that the bounds
        of the axis will be set
    """
    ymin, ymax = input_ax.get_ylim()
    maxvalue = max([abs(_) for _ in input_ax.get_ylim()])
    percent_over = maxvalue * percent_over
    if abs(round(ymin / ymax)) != 0:
        # Set bounds to be the same positive and negative
        bounds = (-1 * (maxvalue + percent_over), (maxvalue + percent_over))
    else: 
        bounds = (-0.05, 1.05)
    input_ax.set_ylim(bounds)

