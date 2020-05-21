#!/usr/bin/env python3
"""
Functions used to create standard statistical plots for the Inspector class
This is the Inspector's Gadgte ;)
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa.utils.calculate import normalize_a_to_b


class InspectorPlotter:
    """
    A class of methods for plotting statistics from an Inspector.
    Should not be called on its own, these functions will be inherited by
    the Inspector class automatically.
    """
    def map(self, show=True, save=False, **kwargs):
        """
        Plot source and receiver locations with map view. Optional arguments
        for only plotting certain stations or events.

        :type show: bool
        :param show: Show the plot
        :type save: str
        :param save: fid to save the given figure
        """
        f, ax = plt.subplots()

        if not self.sources.empty:
            sc_sources = plt.scatter(self.sources.longitude.to_numpy(),
                                     self.sources.latitude.to_numpy(),
                                     marker="o", c="w", edgecolors="g", s=10,
                                     zorder=100
                                     )
        if not self.receivers.empty:
            sc_receiver_names, sc_receiver_list = [], []

            # Color by unique network values
            networks = self.receivers.index.get_level_values(
                "network").unique().to_numpy()
            for net in networks:
                # Random color cycle for networks
                sc_receivers = plt.scatter(
                    self.receivers.loc[net].longitude.to_numpy(),
                    self.receivers.loc[net].latitude.to_numpy(),
                    marker="v", s=10, zorder=100, label=net
                )
                sc_receiver_list.append(sc_receivers)
                sc_receiver_names.append(
                    self.receivers.loc[net].index.to_numpy())

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.legend()
        plt.title(f"{len(self.events)} events; {len(self.receivers)} stations")

        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, sc_sources, self.sources.index, **kwargs)
            for rcvs, rcv_names in zip(sc_receiver_list, sc_receiver_names):
                hover_on_plot(f, ax, rcvs, rcv_names, **kwargs)
            plt.show()

        return f, ax
    
    def event_depths(self, xaxis="latitude", show=True, save=None):
        """
        Create a scatter plot of events at depth. Compresses all events onto a
        single slice, optional choice of showing the x-axis or the y-axis

        :type xaxis: str
        :param xaxis: variable to use as the x-axis on the plot
            'latitude' or 'longitude'
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        if xaxis == "latitude":
            x_vals = self.sources.latitude.to_numpy()
        elif xaxis == "longitude":
            x_vals = self.sources.longitude.to_numpy()
        else:
            raise NotImplementedError(
                "'xaxis' must be 'latitude' or 'longitude"
            )

        # Plot initializations
        f, ax = plt.subplots(figsize=(8, 6))
        depths = self.sources.depth_km.to_numpy()
        mags = self.sources.magnitude.to_numpy()
        mags = normalize_a_to_b(mags, 100, 500)
        names = self.sources.index

        # Inverted axis for positive depth values
        if depths[0] > 0:
            plt.gca().invert_yaxis()

        # Scatter plot
        sc = plt.scatter(x_vals, depths, s=mags, c="None", marker="o",
                         edgecolors="k")
        plt.xlabel(xaxis.capitalize())
        plt.ylabel("Depth (km)")
        plt.title(f"N={len(depths)}")
        plt.grid(which="both", linestyle=":", alpha=0.5)
        hover_on_plot(f, ax, sc, names, dissapear=True)

        if save:
            plt.savefig(save)
        if show:
            plt.show

    def hist(self, model, step, model_comp=None, step_comp=None,
             choice="cc_shift_in_seconds", binsize=1., show=True, save=None,
             **kwargs):
        """
        Create a histogram of misfit information for either time shift or
        amplitude differences. Option to compare against different models,
        and to look at different choices.

        Choices are: "misfit", "dlna", "cc_shift_sec", "length_s"

        :type model: str
        :param model: model to choose for misfit
        :type step: str
        :param step: step count to query, e.g. 's00'
        :type model_comp: str
        :param model_comp: model to compare with, will be plotted in front
        :type step_comp: str
        :param step_comp: step to compare with
        :type choice: str
        :param choice: choice of 'cc_shift_s' for time shift, or 'dlnA' as
            amplitude
        :type binsize: float
        :param binsize: size of the histogram bins
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        assert model in self.models, f"model must be in {self.models}"
        assert step in self.steps.loc[model], \
            f"step must be in {self.steps.loc[model]}"
        if model_comp is not None:
            assert model_comp in self.models, \
                f"model_comp must be in {self.models}"
            assert step_comp in self.steps.loc[model_comp], \
                f"step_comp must be in {self.steps.loc[model_comp]}"

        # For fine tuning figure parameters
        title = kwargs.get("title", "")
        xlim = kwargs.get("xlim", None)
        color = kwargs.get("color", "darkorange")
        color_comp = kwargs.get("color_comp", "deepskyblue")
        fontsize = kwargs.get("fontsize", 12)
        figsize = kwargs.get("figsize", (8, 6))
        legend = kwargs.get("legend", True)
        label_range = kwargs.get("label_range", False)
        xstep = kwargs.get("xstep", 2)
        ymax = kwargs.get("ymax", None)

        def get_stats(n_, bins_):
            """get stats from a histogram"""
            mids = 0.5 * (bins_[1:] + bins_[:-1])
            mean = np.average(mids, weights=n_)
            var = np.average((mids - mean)**2, weights=n_)
            std = np.sqrt(var)

            return mean, var, std

        def get_values(m, s):
            """short hand to get the data, and the maximum value in DataFrame"""
            df_a = self.isolate(model=m, step=s)
            try:
                val_ = df_a.loc[:, choice].to_numpy()
            except KeyError as e:
                raise KeyError(f"Inspector.windows has no key {choice}") from e
            lim_ = max(abs(np.floor(min(val_))), abs(np.ceil(max(val_))))
            return val_, lim_

        # For cleaner formatting of axis labels
        label_dict = {"cc_shift_in_seconds": "Time Shift (s)",
                      "dlnA": "$\Delta\ln$(A)",
                      "misfit": "misfit",
                      "length_s": "Window Length (s)"
                      }

        f, ax = plt.subplots(figsize=figsize)
        val, lim = get_values(model, step)

        if model_comp:
            val_comp, lim_comp = get_values(model_comp, step_comp)

            # Reset the limit to be the greater of the two
            lim = max(lim, lim_comp)

            # Compare models, plot original model on top 
            n, bins, patches = plt.hist(
                x=val, bins=np.arange(-1*lim, lim+.1, binsize),
                color=color,  histtype="bar",  edgecolor="black", linewidth=3,
                label=f"{model}{step}; N={len(val)}", zorder=11, alpha=1.
            )
            mu1, var1, std1 = get_stats(n, bins)

            # Plot comparison below
            n2, bins2, patches2 = plt.hist(
                x=val_comp, bins=np.arange(-1*lim, lim+.1, binsize),
                color=color_comp,  histtype="bar", edgecolor="black",
                linewidth=2.5,
                label=f"{model_comp}{step_comp}; N={len(val_comp)}", zorder=10,
            )
            mu2, var2, std2 = get_stats(n2, bins2)

            # Plot edges of comparison over top
            plt.hist(x=val_comp, bins=np.arange(-1*lim, lim+.1, binsize),
                     color="k", histtype="step", edgecolor=color_comp,
                     linewidth=6., zorder=12,
                     )
        else:
            # No comparison model, plot single histogram
            n, bins, patches = plt.hist(
                x=val, bins=len(np.arange(-1*lim, lim, binsize)), color=color,
                histtype="bar", edgecolor="black", linewidth=2.5,
                label=f"{model}; N={len(val)}", zorder=10,
            )
            mu1, var1, std1 = get_stats(n, bins)

        # Set xlimits of the plot
        if xlim:
            plt.xlim(xlim)
        else:
            if choice == "dlna":
                plt.xlim([-1.75, 1.75])
        if ymax:
            plt.ylim([0, ymax])

        # Stats in the title by default
        if not title:
            title = f"$\mu({model})\pm\sigma({model})$={mu1:.2f}$\pm${std1:.2f}"
            if model_comp:
                title += (f"\n$\mu({model_comp})\pm\sigma({model_comp})$="
                          f"{mu2:.2f}$\pm${std2:.2f}"
                          )

        # Finalize plot details
        plt.xlabel(label_dict[choice], fontsize=fontsize)
        plt.ylabel("Count", fontsize=fontsize)
        plt.title(title, fontsize=fontsize)
        plt.tick_params(which='both', direction='in', top=True, right=True,
                        labelsize=fontsize, width=2.)
        if label_range:
            plt.xticks(np.arange(-1*label_range, label_range+.1, step=xstep))
        plt.axvline(x=0, ymin=0, ymax=1, linewidth=2., c="k", zorder=2, 
                    alpha=0.75, linestyle=':')
        if legend:
            plt.legend(fontsize=fontsize/1.25)

        plt.tight_layout() 

        if save:
            plt.savefig(save)
        if show:
            plt.show()
    
        return f, ax

    def plot_windows(self, model, step, event=None, network=None, station=None,
                     component=None):
        """
        Show lengths of windows chosen based on source-receiver distance, akin
        to Tape's Thesis or to the LASIF plots. These are useful for showing
        which phases are chosen, and window choosing behavior as distance
        increases and (probably) misfit increases.

        :type model: str
        :param model: model to analyze
        :type step: str
        :param step: step count to query, e.g. 's00'
        :type component: str
        :param component: choose a specific component to analyze
        """
        comp_dict = {"Z": "orangered",
                     "N": "forestgreen", "R": "forestgreen",
                     "E": "royalblue", "T": "royalblue"
                     }

        srcrcv = self.get_dist_baz()
        windows = self.isolate(model=model, step=step, event=event,
                               network=network, station=station
                               )
        # Only get information up to component, and times
        df = windows.loc[:, ["event", "network", "station", "component",
                             "relative_starttime", "relative_endtime"]
                         ]
        # Merge in information about source-receiver distances
        df = df.merge(srcrcv.drop("backazimuth", axis=1),
                      on=["event", "network", "station"]
                      )
        # Optional filter to look at specific event, or station up to component
        df = df.loc[
            (df["event"] == (event or df["event"].to_numpy())) &
            (df["station"] == (station or df["station"].to_numpy())) &
            (df["network"] == (network or df["network"].to_numpy())) &
            (df["component"] == (component or df["component"].to_numpy()))

            ]

        # Drop unnecessary information from dataframe
        df.drop(["event", "network", "station"], axis=1, inplace=True)
        f, ax = plt.subplots(figsize=(8, 6))
        for window in df.to_numpy():
            comp, start, end, dist = window
            # short time windows show on top
            ax.hlines(y=dist, xmin=start, xmax=end, colors=comp_dict[comp],
                      zorder=10 + 2 * (1 / end - start), alpha=0.4
                      )
            ax.hlines(y=dist, xmin=0, xmax=300, colors="k", alpha=0.1,
                      linewidth=0.1
                      )

        # Empty lines for legend
        for comp in windows.component.unique():
            ax.hlines(y=0, xmin=0, xmax=0.1, color=comp_dict[comp], label=comp)

        plt.title(f"{len(df)} misfit windows")
        plt.xlabel("Time [s]")
        plt.ylabel("Distance [km]")
        plt.legend()
        plt.xlim([min(df.relative_starttime) - 10,
                  max(df.relative_endtime) + 10])
        plt.ylim([min(df.distance_km) - 10, max(df.distance_km) + 10])

        plt.show()

    def convergence(self, by, windows_by="length_s", fontsize=15,
                    show=True, save=None):
        """
        Plot the convergence rate over the course of an inversion.
        Scatter plot of total misfit against model number, or by step count

        :type choice: str
        :param choice: choice to plot convergence through 'model' or 'iter'
        :type windows_by: str
        :param windows_by: parameter to use for Inspector.measurements() to
            determine how to illustrate measurement number, either by
            cum_win_len: cumulative window length in seconds
            num_windows: number of misfit windows
        :type fontsize: int
        :param fontsize: fontsize of all labels
        :type show: bool
        :param show: show the plot after making it
        :type save: str
        :param save: file id to save the figure to
        """
        assert(windows_by in ["n_win", "length_s"]), \
            "windows_by must be: 'n_win; or 'length_s'"

        ydict = {"length_s": "Cumulative Window Length [s]",
                 "n_win": "Numer of Measurements"}

        # Get misfit information and window lengths together
        df = self.misfits()
        df = df.merge(self.nwin(), on=["model", "step"])
        df.drop(["n_event", "summed_misfit"], axis=1, inplace=True)
        models = df.index.get_level_values("model").unique().to_numpy()

        f, ax1 = plt.subplots(figsize=(8, 6))
        ax2 = ax1.twinx()
        # Plot each iteration and every step count, not as intuitive but shows
        # the behavior of the inversion better
        if by == "iteration":
            # Color by unique model names
            start = 0
            for model in models:
                misfits = df.loc[model].misfit.to_numpy()
                windows = df.loc[model, windows_by].to_numpy()
                end = start + len(misfits)
                ax1.plot(np.arange(start, end), misfits, "o-",
                         linewidth=3, markersize=10, zorder=100,
                         )
                ax2.plot(np.arange(start, end), windows, "v--",
                         linewidth=2, markersize=8, zorder=100
                         )
                start = end
            ax1.set_xlabel("Iteration", fontsize=fontsize)
            ax1.set_xticks(np.arange(0, end, 5))
            ax1.set_xticks(np.arange(0, end, 1), minor=True)

        # Plot by the final accepted misfit per model
        elif by == "model":
            xlabels, misfits, windows = [], [], []
            for m, model in enumerate(models):
                try:
                    # Try to access model as the first step
                    df_temp = df.loc[model, "s00"]
                    # xlbl = f"{model}s00"
                except KeyError:
                    # If first step not available, search last step last model
                    df_temp = df.loc[models[m - 1]].iloc[-1]
                    # xlbl = f"{models[m - 1]}{df_temp.name}"
                xlbl = model
                misfit, window = df_temp.loc[["misfit", windows_by]].to_numpy()

                xlabels.append(xlbl)
                misfits.append(misfit)
                windows.append(window)
            ax1.plot(models, misfits, 'o-', linewidth=3, markersize=10, c="k")
            ax2.plot(models, windows, 'v--', linewidth=3, markersize=10, c="k")

            ax1.set_xlabel("Model number", fontsize=fontsize)
            ax1.xaxis.set_ticks(np.arange(0, len(models)), 1)
            ax1.set_xticklabels(xlabels, rotation=45, ha="right")
        else:
            raise ValueError("'by' must be 'model' or 'iteration")

        # Shared formatting
        ax1.set_ylabel("Total Normalized Misfit (solid)", fontsize=fontsize)
        ax2.set_ylabel(f"{ydict[windows_by]} (dashed)", rotation=270,
                       labelpad=15., fontsize=fontsize)
        ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

        # Only set ticks on the x-axis
        ax1.xaxis.grid(True, which="minor", linestyle=":")
        ax1.xaxis.grid(True, which="major", linestyle="-")

        f.tight_layout()

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax1


def colormap_colorbar(cmap, **kwargs):
    """
    Create a custom colormap and colorbar

    :type cmap: matplotlib.colors.ListedColormap
    :param cmap: colormap to use, called like plt.cm.viridis
    :type vmin: float
    :param vmin: min value for colormap
    :type vmax: float
    :param vmax: max value for colormap
    :type dv: float
    :param dv: colormap boundary separations, if None, continuous colorbar
    :type label: str
    :param label: label for colorbar
    :rtype:
    :return:
    """
    vmin = kwargs.get("vmin", 0.)
    vmax = kwargs.get("vmax", 1.)
    dv = kwargs.get("dv", None)
    label = kwargs.get("cbar_label", "")

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    sm.set_clim(vmin, vmax)
    if dv:
        boundaries = np.arange(vmin, vmax, dv)
    else:
        boundaries = None
    cbar = plt.colorbar(sm, boundaries=boundaries, shrink=0.9, pad=0.025)
    if label:
        cbar.ax.set_ylabel(label, rotation=270, labelpad=15)

    return cmap, norm, cbar


def hover_on_plot(f, ax, obj, values, dissapear=True, **kwargs):
    """
    Allow for hover on a plot for custom annotated information

    From Stackoverflow:
        https://stackoverflow.com/questions/7908636/possible-to-make-labels-
        appear-when-hovering-over-a-point-in-matplotlib

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
    """
    # Make some objects to be used for hover-over capabilities
    anno = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                       textcoords="offset points",
                       bbox=dict(boxstyle="round", fc="w"),
                       arrowprops=dict(arrowstyle="->"),
                       zorder=5000
                       )
    anno.set_visible(False)

    def update_anno(ind):
        """Functionality for getting info when hovering over a point
        during an interacting mpl session
        """
        # Choice between a 2D line and a scatter plot
        if isinstance(obj, mpl.lines.Line2D):
            x, y = obj.get_data()
            anno.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
        elif isinstance(obj, mpl.collections.PathCollection):
            pos = obj.get_offsets()[ind["ind"][0]]
            anno.xy = pos
        text = "{}".format("\n".join([values[n] for n in ind["ind"]]))
        anno.set_text(text)
        anno.get_bbox_patch().set_facecolor("w")
        anno.get_bbox_patch().set_alpha(0.5)

    def hover(event):
        """Functionality for getting info when hovering over a point
        during an interacting mpl session
        """
        vis = anno.get_visible()
        if event.inaxes == ax:
            cont, ind = obj.contains(event)
            if cont:
                update_anno(ind)
                anno.set_visible(True)
                f.canvas.draw_idle()
            # This code snippet will make the annotation dissapear when
            # the mouse moves away
            else:
                if vis and dissapear:
                    anno.set_visible(False)
                    f.canvas.draw_idle()

    f.canvas.mpl_connect("motion_notify_event", hover)
    return hover


def annotate_txt(ax, txt, anno_location="lower-right", **kwargs):
    """
    Convenience function to annotate some information

    :type ax: matplot.axes._subplots.AxesSubplot
    :param ax: axis to annotate onto
    :type txt: str
    :param txt: text to annotate
    :type anno_location: str
    :param anno_location: location on the figure to annotate
        available: bottom-right
    :return:
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    if anno_location == "lower-right":
        x = xmin + (xmax - xmin) * 0.675
        y = ymin + (ymax - ymin) * 0.025
        multialignment = "right"
    elif anno_location == "upper-right":
        x = xmin + (xmax - xmin) * 0.675
        y = ymin + (ymax - ymin) * 0.745
        multialignment = "right"
    elif anno_location == "lower-left":
        x = xmin + (xmax - xmin) * 0.050
        y = ymin + (ymax - ymin) * 0.025
        multialignment = "left"
    elif anno_location == "upper-left":
        x = xmin + (xmax - xmin) * 0.050
        y = ymin + (ymax - ymin) * 0.745
        multialignment = "left"

    ax.annotate(s=txt, xy=(x, y), multialignment=multialignment, **kwargs)

