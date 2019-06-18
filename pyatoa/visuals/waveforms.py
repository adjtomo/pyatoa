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
from pyatoa.visuals.plot_utils import align_yaxis, pretty_grids, format_axis

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


def window_maker(st_obs, st_syn, config, **kwargs):
    """
    Plot streams and windows. assumes you have N observation traces and
    N synthetic traces for a 2N length stream object

    NOTE: real hacky way of putting sta/lta and adjoint source on the
    same axis object, normalize them both -1 to 1 and remove the
    mean of the adjoint source to set everything onto 0

    :type st_*: obspy.stream.Stream
    :param st_*: stream object to plot
    :type config: pyatoa.core.config.Config
    """
    # Get keyword arguments
    dpi = kwargs.get("dpi", 100)
    figsize = kwargs.get("figsize", (11.69, 8.27))  #A4
    time_offset = kwargs.get("time_offset", 0)
    windows = kwargs.get("windows", None)
    staltas = kwargs.get("staltas", None)
    adj_srcs = kwargs.get("adj_srcs", None)
    length_s = kwargs.get("length_s", None)
    stalta_wl = kwargs.get("stalta_wl", None)
    unit_output = kwargs.get("unit_output", None)
    total_misfit = kwargs.get("total_misfit", None)
    append_title = kwargs.get("append_title", None)
    show_fig = kwargs.get("show", False)
    save = kwargs.get("save", None)

    # Set some parameters necessary for flexible plotting
    middle_trace = len(st_obs)//2
    unit_dict = {"DISP": "displacement [m]", 
                 "VEL": "velocity [m/s]",
                 "ACC": "acceleration [m/s^2]"}

    # Instantiate plotting instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    axes, twaxes = setup_plot(number_of=len(st_obs), twax=True)
    t = np.linspace(time_offset,
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
                
                # If adjoint sources given, plot and provide legend information
                if adj_srcs is not None:
                    _adj_src = normalize_a_to_b(adj_srcs[comp].adjoint_source,
                                                a=-1, b=1)
                    _adj_src = detrend(_adj_src, type="constant")
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
            comp = "{}\n{}".format(unit_dict[unit_output], comp)

        labels = [l.get_label() for l in lines_for_legend]
        axes[i].legend(lines_for_legend, labels, prop={"size": 9},
                       loc="upper right")
        axes[i].set_ylabel(comp)

        # Dynamic seismogram length
        if not length_s:
            length_s = t[-1]
        axes[i].set_xlim([t[0], length_s])
        twaxes[i].set_yticklabels([])
        for AX in [axes[i], twaxes[i]]:
            format_axis(AX)
        align_yaxis(axes[i], twaxes[i])

    # Set plot title with relevant information
    title = "{net}.{sta} {evt} [{t0}s, {t1}s] ".format(
        net=st_obs[0].stats.network, sta=st_obs[0].stats.station,
        evt=config.event_id, t0=config.min_period, t1=config.max_period
        )
    # append extra information to the title
    if total_misfit:
        title = " ".join([title, "misfit={:.3f}".format(total_misfit)])
    if append_title:
        title = " ".join([title, append_title])
    axes[0].set_title(title)
    
    axes[-1].set_xlabel("time [s]")
    
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show_fig:
        if show_fig == "hold":
            return f
        else:
            plt.show()
    plt.close("all")
    return f
