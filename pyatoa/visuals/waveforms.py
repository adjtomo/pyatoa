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
from pyatoa.visuals.plot_tools import align_yaxis, pretty_grids, format_axis


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

    Note:
        The unit of adjoint source is based on Eq. 57 of Tromp et al. 2005, 
        which is a traveltime adjoint source, and includes the units of the 
        volumentric delta function since it's assumed this is  happening in a 
        3D volume. 

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
    # general kwargs
    figsize = kwargs.get("figsize", (11.689, 8.27))  # A4
    dpi = kwargs.get("dpi", 100)
    z = kwargs.get("zorder", 5)
    legend = kwargs.get("legend", True)

    # kwargs to control the look of the misfit windows and linecolors
    window_color = kwargs.get("window_color", "orange")
    window_anno_fontsize = kwargs.get("window_anno_fontsize", 8)

    # NOTE: window_anno_height is roughly percentage of visible y-axis
    window_anno_height = kwargs.get("window_anno_height", 0.5)  
    window_anno_rotation = kwargs.get("window_anno_rotation", 0)
    window_anno_fontcolor = kwargs.get("window_anno_fontcolor", "k")
    window_anno_fontweight = kwargs.get("window_anno_fontweight", "normal")

    # turn off parts of the plot at will
    plot_windows = kwargs.get("plot_windows", True)
    plot_window_anno = kwargs.get("plot_window_anno", True)
    plot_stalta = kwargs.get("plot_stalta", True)
    plot_adjsrc = kwargs.get("plot_adjsrc", True)
    plot_waterlevel = kwargs.get("plot_waterlevel", True)
    plot_arrivals = kwargs.get("plot_arrivals", True)

    # Trace colors
    obs_color = kwargs.get("obs_color", "k")
    syn_color = kwargs.get("syn_color", "r")
    adj_src_color = kwargs.get("adj_src_color", "g")
    stalta_color = kwargs.get("stalta_color", "gray")

    # allow kwargs to set global rcParams (is this hacky?)
    axes_linewidth = kwargs.get("axes_linewidth", 2)
    fontsize = kwargs.get("fontsize", 12)
    linewidth = kwargs.get("linewidth", 1.6)
    mpl.rcParams['font.size'] = fontsize
    mpl.rcParams['lines.linewidth'] = linewidth
    mpl.rcParams['axes.linewidth'] = axes_linewidth

    # Set some parameters necessary for flexible plotting
    middle_trace = len(st_obs)//2
    unit_dict = {"DISP": "displacement [m]", 
                 "VEL": "velocity [m/s]",
                 "ACC": "acceleration [m/s^2]",
                 "none": ""}
    adj_unit = "m$^{-4}$ s"
    # window_anno_template = "{ccs:.1f}s"
    window_anno_template = ("max_cc={mcc:.2f}\n"
                            "cc_shift={ccs:.2f}s\n"
                            "dlnA={dln:.3f}\n"
                            "left={lft:.1f}s\n"
                            "length={lgt:.1f}s\n"
                            )

    # Legend tag for data-synthetic or synthetic-synthetic
    obs_tag = 'OBS'
    if config.synthetics_only:
        obs_tag = 'TRUE'

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

    # Plot each component in the same fashion
    for i, comp in enumerate(config.component_list):
        obs = st_obs.select(component=comp)
        syn = st_syn.select(component=comp)

        # WAVEFORMS (convention of black for obs, red for syn)
        a1, = axes[i].plot(t, obs[0].data, obs_color, zorder=z,
                           label="{} ({})".format(obs[0].get_id(), obs_tag))
        a2, = axes[i].plot(t, syn[0].data, syn_color, zorder=z,
                           label="{} (SYN)".format(syn[0].get_id()))
        lines_for_legend = [a1, a2]

        # Seismogram length; min starttime of -10s 
        if not length_sec:
            length_sec = t[-1]
        axes[i].set_xlim([np.maximum(time_offset_sec, -10),
                          np.minimum(length_sec, t[-1])
                          ])

        # Bounds for use in setting positions
        xmin, xmax = axes[i].get_xlim()
        ymin, ymax = axes[i].get_ylim()

        # Amplitude ratio criteria
        if config.window_amplitude_ratio > 0:
            threshold_amp = abs(
                config.window_amplitude_ratio * abs_max(obs[0].data)
            )
            # Both negative and positive bounds
            for sign in [-1, 1]:
                axes[i].axhline(y=sign * threshold_amp, xmin=t[0], xmax=t[-1],
                                alpha=0.35, zorder=z - 4, linewidth=1.25, c='k',
                                linestyle=':')
            # Annotate window amplitude ratio
            if i == middle_trace:
                axes[i].annotate(
                    s="{:.0f}% peak amp. obs.".format(
                        config.window_amplitude_ratio * 100), alpha=0.7,
                    xy=(0.85 * (xmax-xmin) + xmin, threshold_amp), fontsize=8,
                )

        # MISFIT WINDOWS
        if (windows is not None) and (comp in windows) and plot_windows:
            for j, window in enumerate(windows[comp]):
                # Misfit windows as rectangle; taken from Pyflex
                tleft = window.left * window.dt + time_offset_sec
                tright = window.right * window.dt + time_offset_sec
                # Ensure that the window clears the y-bounds 
                axes[i].add_patch(Rectangle(
                    xy=(tleft, ymin*2), width=tright-tleft, 
                    height=(ymax-ymin)*5, color=window_color, 
                    alpha=(window.max_cc_value ** 2) * 0.25
                    )
                )
                # Window annotation information: location and text
                t_anno = (tright - tleft) * 0.025 + tleft
                y_anno = window_anno_height * (ymax - ymin) + ymin

                window_anno = window_anno_template.format(
                    mcc=window.max_cc_value,
                    ccs=window.cc_shift * st_obs[0].stats.delta,
                    dln=window.dlnA, lft=tleft, lgt=tright - tleft
                )
                # Annotate window information
                if plot_window_anno:
                    axes[i].annotate(s=window_anno, xy=(t_anno, y_anno),
                                     zorder=z + 1, 
                                     fontsize=window_anno_fontsize,
                                     rotation=window_anno_rotation, 
                                     color=window_anno_fontcolor,
                                     fontweight=window_anno_fontweight
                                     )

            # PHASE ARRIVALS arrivals based on first window phase arrivals
            for phase_arrivals in windows[comp][0].phase_arrivals:
                if phase_arrivals["name"] in ["s", "p"] and plot_arrivals:
                    axes[i].axvline(x=phase_arrivals["time"], ymin=0, ymax=0.15,
                                    color='b', alpha=0.2
                                    )
                    axes[i].annotate(s=phase_arrivals["name"],
                                     xy=(0.975 * phase_arrivals["time"],
                                         0.10 * (ymax-ymin) + ymin),
                                     fontsize=8
                                     )

            # ADJOINT SOURCES, attempted only if misfit windows exist
            adj_src = None
            if (adj_srcs is not None) and (comp in adj_srcs) and plot_adjsrc:
                adj_src = adj_srcs[comp]
                # Time reverse adjoint source; line up with waveforms
                b1, = twaxes[i].plot(
                    t, adj_src.adjoint_source[::-1], adj_src_color, alpha=0.55,
                    linestyle='-.', zorder=z-1,
                    label=f"Adjoint Source ($\chi$={adj_src.misfit:.2f})"
                )
                lines_for_legend += [b1]

        # STA/LTAs
        if (staltas is not None) and (comp in staltas) and plot_stalta:
            # Normalize to the min and max of adjoint source if available 
            if (adj_srcs is not None) and (comp in adj_srcs):
                tymin, tymax = twaxes[i].get_ylim()
                stalta = normalize_a_to_b(staltas[comp], tymin, tymax)
                waterlevel = (tymax - tymin) * stalta_wl + tymin
            # No adjoint source, then plot STA/LTA with its own bounds
            else:
                stalta = staltas[comp]
                waterlevel = stalta_wl
            # STA/LTA, waterlevel and annotate the waterlevel 
            b2, = twaxes[i].plot(t, stalta, stalta_color, alpha=0.4, 
                                 zorder=z - 1, label=f"STA/LTA")
            lines_for_legend += [b2]
            if plot_waterlevel:
                twaxes[i].axhline(y=waterlevel, xmin=t[0], xmax=t[-1],
                                  alpha=0.4, zorder=z - 2, linewidth=1.5, 
                                  c=stalta_color, linestyle='--')
            if i == middle_trace and plot_waterlevel:
                twaxes[i].annotate(s=f"waterlevel = {stalta_wl}", alpha=0.7,
                                   xy=(0.85 * (xmax-xmin) + xmin, waterlevel),
                                   fontsize=8
                                   )

        # Turn off twin ax label if no adjoint source
        else:
            twaxes[i].set_yticklabels([])

        # Middle trace requires additional y-label information
        if i == middle_trace and plot_stalta:
            _label = ""
            if staltas is not None:
                # _label += "sta/lta"
                if adj_srcs is not None and plot_adjsrc:
                    _label += f"\nadjoint source [{adj_unit}]"
                twaxes[i].set_ylabel(_label, rotation=270, labelpad=42.5) #  27.5)
            # middle trace contains units for all traces, overwrite comp var.
            comp = f"{unit_dict[config.unit_output]}\n{comp}"
        axes[i].set_ylabel(comp)

        # LEGEND
        if legend:
            labels = [l.get_label() for l in lines_for_legend]
            axes[i].legend(lines_for_legend, labels, prop={"size": 12},
                           loc="upper right")

        # Format axes and align with waveforms, before plotting other stuff
        for AX in [axes[i], twaxes[i]]:
            format_axis(AX)
        align_yaxis(axes[i], twaxes[i])

    # TITLE with relevant information
    title = f"{st_obs[0].stats.network}.{st_obs[0].stats.station}"

    # Event id may not be given
    if config.event_id is not None:
        title = " ".join([title, config.event_id])
    # Filter bounds to plot title
    title += f" [{config.min_period}s, {config.max_period}s]"
    # User appended title information
    if append_title is not None:
        title = " ".join([title, append_title])
    axes[0].set_title(title)

    # X label
    axes[-1].set_xlabel("time [s]")

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
