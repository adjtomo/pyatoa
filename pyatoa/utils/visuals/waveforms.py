"""
Waveform plotting functionality

Produces waveform plots from stream objects and plots misfit windows
outputted by pyflex as well as adjoint sources from pyadjoint
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.signal import detrend

from pyatoa.utils.operations.calculations import normalize_a_to_b
from pyatoa.utils.visuals.plot_utils import align_yaxis, pretty_grids, \
    format_axis

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def setup_plot(number_of, twax=True):
    """
    Dynamically set up plots according to number_of given

    :type number_of: int
    :param number_of: number of subplots t ogenerate using gridspec
    :type twax: bool
    :param twax: whether or not twin x axes are required
    :rtype (tw)axes: matplotlib axes
    :return (tw)axes: axis objects
    """
    nrows, ncols = number_of, 1
    height_ratios = [1] * number_of
    gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=height_ratios,
                               hspace=0)
    axes, twaxes = [], []
    for i in range(number_of):
        if i == 0:
            ax = plt.subplot(gs[i])
        else:
            ax = plt.subplot(gs[i],sharex=axes[0])
        if twax:
            twinax = ax.twinx()
            twaxes.append(twinax)
        else:
            twax = None
        pretty_grids(ax)
        axes.append(ax)

    # remove x-tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    return axes, twaxes


def window_maker(st_obs, st_syn, config, time_offset_sec=0., windows=None,
                 staltas=None, adj_srcs=None, length_sec=None,
                 append_title=None, dpi=100, figsize=(11.69, 8.27), show=False,
                 save=None):
    """
    Plot streams and windows. assumes you have N observation traces and
    N synthetic traces for a 2N length stream object

    NOTE: real hacky way of putting sta/lta and adjoint source on the
    same axis object, normalize them both -1 to 1 and remove the
    mean of the adjoint source to set everything onto 0

    :type st_obs: obspy.stream.Stream
    :param st_obs: observation stream object to plot
    :type st_syn: obspy.stream.Stream
    :param st_syn: synthetic stream object to plot
    :type config: pyatoa.core.config.Config
    :param config: Configuration that must contain:
        required: component_list, unit_output, event_id, min_period, max_period
        optional: pyflex_config
    :type dpi: int
    :param dpi: dots per inch of figure
    :type figsize: tuple of floats
    :param figsize: size of the waveform plot, defaults to A4 paper in landscape
    :type time_offset_sec: float
    :param time_offset_sec: Offset from origintime. Origintime is set as time=0
    :type windows: dict of pyflex.Window objects
    :param windows: If given, will plot all the misfit windows in the dictionary
    :type staltas: dict of pyflex.sta_lta objects
    :param staltas: Output of stalta analysis to plot in background, if given
    :type adj_srcs: dict of pyadjoint.Adjoint_Source objects
    :param adj_srcs: if given, plots the adjoint sources in the background
    :type length_sec: int
    :param length_sec: sets the seismogram length, useful for short waveforms
        that have lots of zeros in the back
    :type append_title: str
    :param append_title: appends any extra information after the stock title
    :type show: bool or str
    :param show: show the plot after it is created. if 'hold', will just return
        the figure object, so the User can show it later
    :type save: str
    :param save: pathname to save the figure, if not given, will not save
    """
    # Set some parameters necessary for flexible plotting
    middle_trace = len(st_obs)//2
    unit_dict = {"DISP": "displacement [m]", 
                 "VEL": "velocity [m/s]",
                 "ACC": "acceleration [m/s^2]"}
    if staltas:
        stalta_wl = config.pyflex_config[1].stalta_waterlevel

    # Instantiate plotting instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    axes, twaxes = setup_plot(number_of=len(st_obs), twax=True)
    t = np.linspace(time_offset_sec,
                    st_obs[0].stats.endtime-st_obs[0].stats.starttime,
                    len(st_obs[0].data))
    
    # Plot per component
    z = 5
    for i, comp in enumerate(config.component_list):
        obs = st_obs.select(component=comp)
        syn = st_syn.select(component=comp)
        a1, = axes[i].plot(t, obs[0].data, 'k', zorder=z,
                           label="{} (OBS)".format(obs[0].get_id()))
        a2, = axes[i].plot(t, syn[0].data, 'r', zorder=z,
                           label="{} (SYN)".format(syn[0].get_id()))
        lines_for_legend = [a1, a2]
        
        # If STA/LTA information from Pyflex given, plot behind waveforms
        if staltas is not None:
            stalta_ = normalize_a_to_b(staltas[comp], -1, 1)
            twaxes[i].plot(t, stalta_, 'gray', alpha=0.4, zorder=z-1)
            twaxes[i].axhline(y=stalta_wl-1, xmin=t[0], xmax=t[-1],
                              alpha=0.2, zorder=z-2, linewidth=1.5, c='k',
                              linestyle='--')

        # If misfit windows given, plot each window and annotate information
        if windows is not None:
            ymin, ymax = axes[i].get_ylim()
            window_anno = "maxCC:{mcc:.4f}\nccShift:{ccs:.2f}s\ndlnA:{dln:.4f}"
            try:
                for window in windows[comp]:
                    xwindow = np.arange(window.left, window.right, 1)
                    twindow = t[xwindow]
                    axes[i].fill_between(twindow, ymin, ymax,
                                         facecolor='orange', edgecolor='k',
                                         linewidth=0.5,
                                         alpha=(window.max_cc_value ** 2) * 0.25
                                         )
                    anno_str = window_anno.format(
                        mcc=window.max_cc_value,
                        ccs=window.cc_shift * st_obs[0].stats.delta,
                        dln=window.dlnA)
                    axes[i].annotate(s=anno_str, xy=(twindow[10], ymax * 0.5),
                                     zorder=z-1, fontsize=7)
                
                # If Adjoint Sources given, plot and provide legend information
                if adj_srcs is not None:
                    _adj_src = normalize_a_to_b(adj_srcs[comp].adjoint_source,
                                                a=-1, b=1)
                    _adj_src = detrend(_adj_src, type="constant")

                    # Plot adjoint source time reversed, line up with waveforms
                    b1, = twaxes[i].plot(t, _adj_src[::-1], 'g', alpha=0.5,
                                         linestyle='-.',
                                         label="Adjoint Source, Misfit={:.4f}".
                                         format(adj_srcs[comp].misfit)
                                         )
                    lines_for_legend += [b1]
            except KeyError:
                pass

        # The y-label of the middle trace contains common info from Pyflex
        if i == middle_trace:
            twaxes[i].set_ylabel("adjoint source (normalized) &\n"
                                 "STA/LTA (waterlevel = {})".format(stalta_wl),
                                 rotation=270, labelpad=20)
            comp = "{}\n{}".format(unit_dict[config.unit_output], comp)

        labels = [l.get_label() for l in lines_for_legend]
        axes[i].legend(lines_for_legend, labels, prop={"size": 9},
                       loc="upper right")
        axes[i].set_ylabel(comp)

        # Dynamic seismogram length
        if not length_sec:
            length_sec = t[-1]
        axes[i].set_xlim([t[0], length_sec])

        # Format axes and align
        twaxes[i].set_yticklabels([])
        for AX in [axes[i], twaxes[i]]:
            format_axis(AX)
        align_yaxis(axes[i], twaxes[i])

    # Set plot title with relevant information
    title = "{net}.{sta} {evt} [{t0}s, {t1}s] ".format(
        net=st_obs[0].stats.network, sta=st_obs[0].stats.station,
        evt=config.event_id, t0=config.min_period, t1=config.max_period
        )
    # Append extra information to title
    if append_title is not None:
        title = " ".join([title, append_title])

    axes[0].set_title(title)
    axes[-1].set_xlabel("time [s]")
    
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show:
        if show == "hold":
            return f
        else:
            plt.show()
    plt.close("all")
    return f
