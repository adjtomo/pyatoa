"""
Waveform plotting functionality

Produces waveform plots from stream objects and plots misfit windows
outputted by pyflex as well as adjoint sources from pyadjoint
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle

from pyatoa.utils.tools.calculate import normalize_a_to_b, abs_max
from pyatoa.utils.visuals.plot_tools import align_yaxis, pretty_grids, \
    format_axis

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.6
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
            ax = plt.subplot(gs[i], sharex=axes[0])
        if twax:
            twinax = ax.twinx()
            pretty_grids(twinax, twax=True)
            twaxes.append(twinax)
        else:
            twaxes = None
        pretty_grids(ax)
        axes.append(ax)

    # remove x-tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    return axes, twaxes


def window_maker(st_obs, st_syn, config, time_offset_sec=0., windows=None,
                 staltas=None, adj_srcs=None, length_sec=None,
                 append_title=None, show=False, save=None, **kwargs):
    """
    Plot streams and windows. assumes you have N observation traces and
    N synthetic traces for a 2N length stream object

    :type st_obs: obspy.stream.Stream
    :param st_obs: observation stream object to plot
    :type st_syn: obspy.stream.Stream
    :param st_syn: synthetic stream object to plot
    :type config: pyatoa.core.config.Config
    :param config: Configuration that must contain:
        required: component_list, unit_output, event_id, min_period, max_period
        optional: pyflex_config
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
    figsize = kwargs.get("figsize", (11.689, 8.27))  # A4
    dpi = kwargs.get("dpi", 100)

    # Set some parameters necessary for flexible plotting
    middle_trace = len(st_obs)//2
    unit_dict = {"DISP": "displacement [m]", 
                 "VEL": "velocity [m/s]",
                 "ACC": "acceleration [m/s^2]"}
    adj_dict = {"DISP": "[m^-1]",
                "VEL": "[m^-1 s]",
                "ACC": "[m^-s s^-2]"}

    # Legend tagging to remove confusion during synthetic-synthetic case
    if config.synthetics_only:
        obs_tag = 'TRUE'
    else:
        obs_tag = 'OBS'

    if staltas is not None:
        stalta_wl = config.pyflex_config.stalta_waterlevel

    # Instantiate plotting instances
    f = plt.figure(figsize=figsize, dpi=dpi)
    axes, twaxes = setup_plot(number_of=len(st_obs), twax=True)
    t = np.linspace(
        time_offset_sec,
        st_obs[0].stats.endtime-st_obs[0].stats.starttime+time_offset_sec,
        len(st_obs[0].data)
        )
    
    z = 5
    # Set the annotation for misfit windows
    window_anno_template = ("max_cc={mcc:.2f}\n"
                            "cc_shift={ccs:.2f}s\n"
                            "dlnA={dln:.3f}\n"
                            "left={lft:.1f}s\n"
                            "length={lgt:.1f}s\n"
                            )
    # Plot per component
    for i, comp in enumerate(config.component_list):
        # Plot waveforms following the convention of black for obs, red for syn
        obs = st_obs.select(component=comp)
        syn = st_syn.select(component=comp)

        a1, = axes[i].plot(t, obs[0].data, 'k', zorder=z,
                           label="{} ({})".format(obs[0].get_id(), obs_tag))
        a2, = axes[i].plot(t, syn[0].data, 'r', zorder=z,
                           label="{} (SYN)".format(syn[0].get_id()))
        lines_for_legend = [a1, a2]

        # Set the seismogram length, min starttime of -10 seconds before origin
        if not length_sec:
            length_sec = t[-1]
        axes[i].set_xlim([np.maximum(time_offset_sec, -10),
                          np.minimum(length_sec, t[-1])
                          ])
        # Get plot bounds for use in setting objects at correct positions
        xmin, xmax = axes[i].get_xlim()
        ymin, ymax = axes[i].get_ylim()

        # Plot the amplitude ratio criteria if given by User
        if config.window_amplitude_ratio > 0:
            threshold_amp = abs(
                config.window_amplitude_ratio * abs_max(obs[0].data)
            )
            # Plot both negative and positive bounds
            for sign in [-1, 1]:
                axes[i].axhline(y=sign * threshold_amp, xmin=t[0], xmax=t[-1],
                                alpha=0.35, zorder=z - 4, linewidth=1.25, c='k',
                                linestyle=':')
            # Annotate the window amplitude ratio specified by User
            if i == middle_trace:
                axes[i].annotate(
                    s="{:.0f}% peak amp. obs.".format(
                        config.window_amplitude_ratio * 100), alpha=0.7,
                    xy=(0.85 * (xmax-xmin) + xmin, threshold_amp), fontsize=8,
                )

        # If misfit windows given, plot each window and annotate information
        if (windows is not None) and (comp in windows.keys()):
            for j, window in enumerate(windows[comp]):
                # Rectangle to represent misfit windows; taken from Pyflex
                tleft = window.left * window.dt + time_offset_sec
                tright = window.right * window.dt + time_offset_sec
                axes[i].add_patch(Rectangle(
                    xy=(tleft, ymin), width=tright-tleft, height=ymax-ymin,
                    color='orange', alpha=(window.max_cc_value ** 2) * 0.25)
                )
                # Annotate window information into the windows
                t_anno = (tright - tleft) * 0.025 + tleft
                y_anno = ymax * 0.2
                window_anno = window_anno_template.format(
                    mcc=window.max_cc_value,
                    ccs=window.cc_shift * st_obs[0].stats.delta,
                    dln=window.dlnA, lft=tleft, lgt=tright - tleft
                )
                # If an adjoint source is given, add info from measurement
                # Adjoint source measurements should line up with the misfit win
                if adj_srcs is not None:
                    y_anno = ymax * 0.005
                    window_anno = "type={typ}\nmsft={mft:.1E}\n".format(
                        typ=adj_srcs[comp].measurement[j]["type"],
                        mft=adj_srcs[comp].measurement[j]["misfit_dt"]
                    ) + window_anno
                axes[i].annotate(s=window_anno, xy=(t_anno, y_anno),
                                 zorder=z + 1, fontsize=8,
                                 )

            # Plot the direct arrivals from the first window in each subplot
            # as vertical ticks in blue
            for phase_arrivals in windows[comp][0].phase_arrivals:
                if phase_arrivals["name"] in ["s", "p"]:
                    axes[i].axvline(x=phase_arrivals["time"], ymin=0, ymax=0.15,
                                    color='b', alpha=0.2
                                    )
                    axes[i].annotate(s=phase_arrivals["name"],
                                     xy=(0.975 * phase_arrivals["time"],
                                         0.10 * (ymax-ymin) + ymin),
                                     fontsize=8
                                     )

            # If Adjoint Sources given, plot and provide legend information
            adj_src = None
            if (adj_srcs is not None) and (comp in adj_srcs.keys()):
                adj_src = adj_srcs[comp]

                # Plot adjoint source time reversed, line up with waveforms
                b1, = twaxes[i].plot(
                    t, adj_src.adjoint_source[::-1], 'g', alpha=0.55,
                    linestyle='-.', zorder=z-1,
                    label="Adjoint Source, Misfit={:.2E}".format(adj_src.misfit)
                )
                lines_for_legend += [b1]

            # If STA/LTA information from Pyflex given, plot behind waveforms
            if (staltas is not None) and (comp in staltas.keys()):
                # Normalize to the min and max of the adjoint source
                if adj_src is not None:
                    tymin, tymax = twaxes[i].get_ylim()
                    stalta = normalize_a_to_b(staltas[comp], tymin, tymax)
                    waterlevel = (tymax - tymin) * stalta_wl + tymin
                # If no adjoint source, simply plot STA/LTA with its own bounds
                else:
                    stalta = staltas[comp]
                    waterlevel = stalta_wl
                # Plot the STA/LTA, waterlevel and annotate the waterlevel %
                twaxes[i].plot(t, stalta, 'gray', alpha=0.4, zorder=z - 1)
                twaxes[i].axhline(y=waterlevel, xmin=t[0], xmax=t[-1],
                                  alpha=0.2, zorder=z - 2, linewidth=1.5, c='k',
                                  linestyle='--')
                if i == middle_trace:
                    twaxes[i].annotate(
                        s="waterlevel = {}".format(stalta_wl), alpha=0.7,
                        xy=(0.85 * (xmax-xmin) + xmin, waterlevel), fontsize=8
                    )

        # If no adjoint source for given component, turn off twin ax label
        else:
            twaxes[i].set_yticklabels([])

        # The y-label of the middle trace contains common info
        if i == middle_trace:
            # If STA/LTA information given, create a custom label
            if staltas is not None:
                _label = "sta/lta"
                if adj_srcs is not None:
                    _label += "\nadjoint source {}".format\
                        (adj_dict[config.unit_output]
                         )
                twaxes[i].set_ylabel(_label, rotation=270, labelpad=27.5)

            # Middle trace contains units for all traces, overwrite comp var.
            comp = "{}\n{}".format(unit_dict[config.unit_output], comp)

        axes[i].set_ylabel(comp)
        # Accumulate legend labels
        labels = [l.get_label() for l in lines_for_legend]
        axes[i].legend(lines_for_legend, labels, prop={"size": 9},
                       loc="upper right")

        # Format axes and align
        for AX in [axes[i], twaxes[i]]:
            format_axis(AX)
        align_yaxis(axes[i], twaxes[i])

    # Set plot title with relevant information
    title = "{net}.{sta}".format(net=st_obs[0].stats.network,
                                 sta=st_obs[0].stats.station
                                 )

    # Account for the fact that event id may not be given
    if config.event_id is not None:
        title = " ".join([title, config.event_id])

    # Add filter bounds to plot title
    title += " [{0}s, {1}s]".format(config.min_period, config.max_period)

    # Append extra information to title
    if append_title is not None:
        title = " ".join([title, append_title])

    axes[0].set_title(title)
    axes[-1].set_xlabel("time [s]")

    # Make sure to remove the extra whitespace before saving/showing
    f.tight_layout()
    if save:
        plt.savefig(save, figsize=figsize, dpi=dpi)
    if show:
        if show == "hold":
            return f
        else:
            plt.show()
    plt.close("all")
    return f
