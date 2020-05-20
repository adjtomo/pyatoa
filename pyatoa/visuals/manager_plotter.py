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
from pyatoa.utils.form import format_model_number, format_step_count
from pyatoa.utils.calculate import normalize_a_to_b, abs_max


# Hardcoded colors that represent rejected misfit windows. Description from 
# pyflex.WindowSelector function that rejected the window
rejected_window_colors = {
    "schedule": "C0",  # Blue: schedule_weighted_interval()
    "water_level": "C1",  # Orange: reject_on_minima_water_level()
    "prominence": "C2",   # Green: reject_on_prominence_of_central_peak()
    "phase_sep": "C3",   # Red: reject_on_phase_separation()
    "curtail": "C4",   # Violet: curtail_length_of_windows()
    "min_length": "C5",  # Brown: reject_windows_based_on_minimum_length()
    "data_fit": "C6",   # Pink: reject_based_on_data_fit_criteria()
    "s2n": "C7",  # Gray: reject_based_on_signal_to_noise_ratio()
    "traveltimes":"C8",  # Olive: reject_on_traveltimes()
    "amplitude": "C9"  # Cyan: Pyatoa.utils.window.reject_on_global_amplitude_ratio()
    }



class ManagerPlotter:
    """
    Standardized waveform figures featuring observed and synthetic traces, 
    STA/LTA waveforms, misfit windows, rejected windows, adjoint sources and
    auxiliary information collected within the workflow.
    """
    def __init__(self, mgmt, show=True, save=None, **kwargs):
        """
        Introduce class-wide objects that accessed for plotting the figure

        :Keyword Arguments:
            figsize (tuple): size of the figure, defaults A4 (11.689, 8.27)
            dpi (float): dots per inch of the figure, defaults 100
            linewidth (float): line width for all lines on plot, default 1.6

            obs_color (str): color for observed waveform, defaults 'k'
            syn_color (str): color for synthetic waveform, defaults 'r
            stalta_color (str): color of stalta waveform, default 'gray'
            window_color (str): color for misfit windows, default 'orange'
            adj_src_color (str): color for adjoint sources, default 'g'

            window_anno_fontsize (str): fontsize for window annotation, def 8
            window_anno_height (float): annotation height, percentage of y axis,
                default 0.7
            window_anno_rotation (float): rotation of annotation (deg), def 0
            window_anno_fontcolor (str): color of annotation text, def 'k'
            window_anno_fontweight (str): weight of font, default 'normal' 

            plot_windows (bool): toggle window plotting, default True
            plot_rejected_windows (bool): toggle rejected window plot, default T
            plot_window_annos (bool): toggle window annotations, default True
            plot_staltas (bool): toggle stalta plotting, default True
            plot_adjsrcs (bool): toggle adjoint source plotting, default True
            plot_waterlevel (bool): toggle stalta waterlevel plotting, def True
            plot_arrivals (bool): toggle phase arrival plotting, default True
            plot_legend (bool): toggle legend, default True

            normalize (bool): normalize waveform data before plotting
            append_title(str): User appended string to the end of the title
        """
        self.st_obs = mgmt.st_obs.copy()
        self.st_syn = mgmt.st_syn.copy()
        self.config = mgmt.config
        self.windows = mgmt.windows
        self.staltas = mgmt.staltas
        self.adjsrcs = mgmt.adj_srcs
        self.rejected_windows = mgmt._rej_win

        self.show = show
        self.save = save
        self.kwargs = kwargs

        self.time_axis = self.st_obs[0].times() + mgmt.stats.time_offset_sec

    def setup_plot(self):
        """
        Dynamically set up plots according to number_of given

        :type number_of: int
        :param number_of: number of subplots to generate using gridspec
        :type twax: bool
        :param twax: whether or not twin x axes are required
        :rtype (tw)axes: matplotlib axes
        :return (tw)axes: axis objects
        """
        # general kwargs
        figsize = self.kwargs.get("figsize", (11.689, 8.27))  # A4
        dpi = self.kwargs.get("dpi", 100)

        # Initiate the figure and fill it up with grids
        f = plt.figure(figsize=figsize, dpi=dpi)

        nrows, ncols = len(self.st_obs), 1
        height_ratios = [1] * len(self.st_obs)
        gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=height_ratios,
                                   hspace=0)
        axes, twaxes = [], []
        for i in range(len(self.st_obs)):
            if i == 0:
                ax = plt.subplot(gs[i])
            else:
                ax = plt.subplot(gs[i], sharex=axes[0])
            twinax = ax.twinx()

            pretty_grids(twinax, twax=True)
            pretty_grids(ax)
            twaxes.append(twinax)
            axes.append(ax)

        # remove x-tick labels except for last axis
        for ax in axes[0:-1]:
            plt.setp(ax.get_xticklabels(), visible=False)

        return f, axes, twaxes

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
        """
        obs_color = self.kwargs.get("obs_color", "k")
        syn_color = self.kwargs.get("syn_color", "r")
        linewidth = self.kwargs.get("linewdith", 1.6)

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
        a1, = ax.plot(self.time_axis, obs.data, obs_color, zorder=10, 
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
        :type normalize: bool
        :param normalize: normalize the trace to the bounds of the axis before
            plotting. used if adjoint source already plotted on the axis, stalta
            will have different units so normalize against the axis because 
            we only care about phases not amplitudes
        :type plot_waterlevel: bool
        :param plot_waterlevel: plot a horizontal line showing the relative
            waterlevel of the sta/lta which is used in determining windows
        :rtype: list of matplotlib.lines.Lines2D objects
        :return: list containing the lines plotted on the axis
        """
        linewidth = self.kwargs.get("linewdith", 1.6)
        stalta_color = self.kwargs.get("stalta_color", "gray")

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
                      zorder=9, label=f"STA/LTA")

        if plot_waterlevel:
            # Plot the waterlevel of the STA/LTA defined by Pyflex Config
            ax.axhline(y=waterlevel, xmin=self.time_axis[0], 
                       xmax=self.time_axis[-1], alpha=0.4, zorder=8, 
                       linewidth=linewidth, c=stalta_color, linestyle='--')
            ax.annotate(s=f"waterlevel = {stalta_wl}", alpha=0.7,
                        xy=(0.85 * (xmax - xmin) + xmin, waterlevel), fontsize=8
                        )

        return [b2]

    def plot_windows(self, ax, windows, plot_window_annos=True, 
                     plot_phase_arrivals=True):
        """
        Plot misfit windows, add annotations to each window related to 
        information contained in the Window object.

        Note:
         -kwarg 'window_anno_height' should be given as a percentage of 
          visible y-axis

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
        window_color = self.kwargs.get("window_color", "orange")
        window_anno_fontsize = self.kwargs.get("window_anno_fontsize", 8)
        window_anno_height = self.kwargs.get("window_anno_height", 0.65)  
        window_anno_rotation = self.kwargs.get("window_anno_rotation", 0)
        window_anno_fontcolor = self.kwargs.get("window_anno_fontcolor", "k")
        window_anno_fontweight = self.kwargs.get("window_anno_fontweight", 
                                                 "normal")

        ymin, ymax = ax.get_ylim()

        for j, window in enumerate(windows):
            tleft = window.left * window.dt + self.time_axis[0]
            tright = window.right * window.dt + self.time_axis[0]

            # Misfit windows as rectangle; taken from Pyflex
            ax.add_patch(Rectangle(xy=(tleft, ymin), width=tright - tleft, 
                                   height=(ymax + np.abs(ymin)), 
                                   fc=window_color, ec="k", 
                                   alpha=(window.max_cc_value ** 2) * 0.25
                )
            )
            if plot_window_annos:
                # Annotate window information into each window
                t_anno = (tright - tleft) * 0.025 + tleft
                y_anno = window_anno_height * (ymax - ymin) + ymin

                window_anno = (
                    f"max_cc={window.max_cc_value:.2f}\n"
                    f"cc_shift={window.cc_shift * window.dt:.2f}s\n"
                    f"dlnA={window.dlnA:.3f}\n"
                    f"left={tleft:.1f}s\n"
                    f"length={tright - tleft:.1f}s\n"
                    )
                ax.annotate(s=window_anno, xy=(t_anno, y_anno), zorder=11, 
                            fontsize=window_anno_fontsize,
                            rotation=window_anno_rotation, 
                            color=window_anno_fontcolor,
                            fontweight=window_anno_fontweight
                            )

            if plot_phase_arrivals:
                for phase_arrivals in window.phase_arrivals:
                    if phase_arrivals["name"] in ["p", "s"]:
                        ax.axvline(x=phase_arrivals["time"], ymin=0, ymax=0.15,
                                   color='b', alpha=0.2
                                   )
                        ax.annotate(s=phase_arrivals["name"],
                                    xy=(0.975 * phase_arrivals["time"], 
                                        0.10 * (ymax-ymin) + ymin),
                                    fontsize=8
                                    )

    def plot_rejected_windows(self, ax, windows):
        """
        Plot rejected windows as solid transparent lines out of the way of the 
        plot. Hardcoded color dictionary used as a way to visually identify why
        certain windows were rejected

        :type ax: matplotlib.axes.Axes
        :param ax: axis object on which to plot
        :type windows: list of pyflex.Window objects
        :param windows: list of windows to plot
        """
        ymin, ymax = ax.get_ylim()
        dy = 0.025 * (ymax - ymin)  # increment for window location

        for tag, window in windows.items():
            for win in window:
                tleft = win.left * win.dt + self.time_axis[0]
                tright = win.right * win.dt + self.time_axis[0]
                ax.hlines(y=ymin, xmin=tleft, xmax=tright, linewidth=1.75,
                          colors=rejected_window_colors[tag], alpha=1.)
                # ax.annotate(xy=(tleft, ymin), s=tag, fontsize=3)
                ymin -= dy

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
        linewidth = self.kwargs.get("linewdith", 1.6)
        adj_src_color = self.kwargs.get("adj_src_color", "g")

        # Time reverse adjoint source; line up with waveforms
        b1, = ax.plot(self.time_axis, adjsrc.adjoint_source[::-1], 
                      adj_src_color, alpha=0.55, linewidth=linewidth,
                      linestyle="-.", zorder=9,
                      label=fr"Adjoint Source ($\chi$={adjsrc.misfit:.2E})"
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
        ax.annotate(s=f"{config.win_amp_ratio * 100:.0f}% peak amp. obs.", 
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
        if self.config.model is not None:
            title += format_model_number(self.config.model)
        if self.config.step is not None:
            title += format_step_count(self.config.step)

        # Add information about the Pyflex and Pyadjoint parameters used
        if self.kwargs.get("plot_stalta", True):
            title += f" pyflex={self.config.pyflex_preset} "
        if self.kwargs.get("plot_adjsrc", True):
            title += f" pyadjoint={self.config.adj_src_type} "

        # User appended title information
        if append_title is not None:
            title = " ".join([title, append_title])

        return title

    def plot(self):
        """
        High level plotting function that plots all parts of the class and 
        formats the axes nicely
        """
        # Distribute some kwargs before starting
        figsize = self.kwargs.get("figsize", (11.689, 8.27))  # A4
        dpi = self.kwargs.get("dpi", 100)
        append_title = self.kwargs.get("append_title", None)
        normalize = self.kwargs.get("normalize", False)

        plot_windows = self.kwargs.get("plot_windows", True)
        plot_rejected_windows = self.kwargs.get("plot_rejected_windows", True)
        plot_window_annos = self.kwargs.get("plot_window_anno", True)
        plot_staltas = self.kwargs.get("plot_stalta", True)
        plot_adjsrcs = self.kwargs.get("plot_adjsrc", True)
        plot_waterlevel = self.kwargs.get("plot_waterlevel", True)
        plot_arrivals = self.kwargs.get("plot_arrivals", True)
        plot_legend = self.kwargs.get("legend", True)

        # Plot per component in the same fashion, only if observed data exists
        f, axes, twaxes = self.setup_plot()
        for i, obs in enumerate(self.st_obs):
            comp = obs.stats.component

            # Try to retrieve auxiliary data by component name
            syn = self.st_syn.select(component=comp)[0]
            windows = self.windows.get(comp, None)
            adjsrc = self.adjsrcs.get(comp, None)
            stalta = self.staltas.get(comp, None)
            rejwins = self.rejected_windows.get(comp, None)

            # Begin plotting by distributing axis objects
            ax = axes[i]
            twax = twaxes[i]
            lines = []  # List of lines for making the legend
            lines += self.plot_waveforms(obs=obs, syn=syn, ax=ax, 
                                         normalize=normalize)

            if rejwins is not None and plot_rejected_windows:
                self.plot_rejected_windows(ax=ax, windows=rejwins)

            # Format now as windows use the y-limits for height of windows
            format_axis(ax)

            if windows is not None and plot_windows: 
                self.plot_windows(ax=ax, windows=windows, 
                                  plot_window_annos=plot_window_annos, 
                                  plot_phase_arrivals=plot_arrivals)

            if adjsrc is not None and plot_adjsrcs:
                lines += self.plot_adjsrcs(ax=twax, adjsrc=adjsrc)
                if i == len(self.st_obs) // 2:  
                    # middle trace: append units of the adjoint source on ylabel
                    twax.set_ylabel("adjoint source [m$^{-4}$ s]", rotation=270, 
                                    labelpad=42.5)
            else:
                twax.set_yticklabels([])  # turn off yticks if no adjsrc

            # Format twax because stalta will use y-limits for its waveforms
            format_axis(twax)

            if stalta is not None and plot_staltas:
                lines += self.plot_stalta(ax=twax, stalta=stalta)
                if self.config.win_amp_ratio > 0:
                    plot_amplitude_threshold(ax=ax, obs=obs)

            # Finish with axes formatting
            if i == len(self.st_obs) // 2:
                # Middle trace will carry the units of the waveforms
                units = {"DISP": "displacement [m]", 
                         "VEL": "velocity [m/s]",
                         "ACC": "acceleration [m/s^2]",
                         "none": ""}[self.config.unit_output]
                ax.set_ylabel(f"{units}\n{comp}")
            else:
                ax.set_ylabel(comp)

            if plot_legend:
                labels = [l.get_label() for l in lines]
                ax.legend(lines, labels, prop={"size": 10}, loc="upper right")

            align_yaxes(ax, twax)

        # Final touch ups for the figure
        axes[0].set_title(self.create_title(append_title=append_title,
                                            normalized=normalize))
        axes[0].set_xlim([np.maximum(self.time_axis[0], -10), 
                          self.time_axis[-1]]
                          )
        axes[-1].set_xlabel("time [s]")
        f.tight_layout()

        if self.save:
            plt.savefig(self.save, figsize=figsize, dpi=dpi)
        if self.show:
            plt.show()
        else:
            plt.close()


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


def pretty_grids(input_ax, twax=False, grid=False):
    """
    Standard plot skeleton formatting, thick lines and internal tick marks etc.

    :type input_ax: matplotlib axis
    :param input_ax: axis to prettify
    :type twax: bool
    :param twax: If twax (twin axis), do not set grids
    """
    input_ax.set_axisbelow(True)
    input_ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    input_ax.tick_params(which='both', direction='in', top=True, right=True)

    # Set the grids 'on' only if main axis
    if not twax:
        input_ax.minorticks_on()
        if grid:
            for axis_ in ['major', 'minor']:
                input_ax.grid(which=axis_, linestyle=':', linewidth='0.5',
                              color='k', alpha=0.25)


def format_axis(input_ax):
    """
    Ensure that the upper and lower y-bounds are the same value

    :type input_ax: matplotlib axis
    :param input_ax: axis to prettify
    """
    ymin, ymax = input_ax.get_ylim()
    maxvalue = max([abs(_) for _ in input_ax.get_ylim()])
    percentover = maxvalue * 0.125
    if abs(round(ymin / ymax)) != 0:
        # Set bounds to be the same positive and negative
        bounds = (-1 * (maxvalue + percentover), (maxvalue + percentover))
    else: 
        bounds = (-0.05, 1.05)
    input_ax.set_ylim(bounds)

