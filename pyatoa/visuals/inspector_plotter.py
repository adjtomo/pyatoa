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

    def map(self, event=None, network=None, station=None, show=True, save=False,
            **kwargs):
        """
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
        """
        # For isolating parameters, ensure they are lists. quack
        if isinstance(event, str):
            event = [event]
        if isinstance(network, str):
            network = [network]
        if isinstance(station, str):
            station = [station]

        markersize = kwargs.get("markersize", 10)
        f, ax = plt.subplots()

        if not self.sources.empty:
            # Allow for isolation of particular events
            if event is not None:
                src_lat = self.sources.loc[event].latitude
                src_lon = self.sources.loc[event].longitude
            else:
                src_lat = self.sources.latitude
                src_lon = self.sources.longitude
            sc_sources = plt.scatter(src_lat, src_lon, marker="o", c="None",
                                     edgecolors="k", s=markersize, zorder=100
                                     )
        if not self.receivers.empty:
            sc_receiver_names, sc_receiver_list = [], []
            # Allow for isolation of networks
            if network is not None:
                networks = self.receivers.loc[network]
            else:
                networks = self.receivers
            networks = networks.index.get_level_values("network").unique()
            for net in networks:
                # Allow for isolation of stations
                if station is not None:
                    try:
                        rcv_lat = self.receivers.loc[net].loc[station].latitude
                        rcv_lon = self.receivers.loc[net].loc[station].longitude
                        rcv_nam = station
                    except KeyError:
                        continue
                else:
                    # else just plot all stations in a given network
                    rcv_lat = self.receivers.loc[net].latitude
                    rcv_lon = self.receivers.loc[net].longitude
                    rcv_nam = self.receivers.loc[net].index.to_numpy()

                # Random color cycle for networks
                sc_receivers = plt.scatter(rcv_lat, rcv_lon, marker="v",
                                           s=markersize, zorder=100, label=net)
                sc_receiver_list.append(sc_receivers)
                sc_receiver_names.append(rcv_nam)

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
        else:
            plt.close()

    def raypaths(self, iteration, step_count, show=True, save=False, **kwargs):
        """
        Plot rays connecting sources and receivers based on the availability
        of measurements. Useful for getting an approximation of resolution.

        :type iteration: int
        :param iteration: iteration to retrieve data from
        :type step_count: int
        :param step_count: step count to retrieve data from
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        ray_color = kwargs.get("ray_color", "k")
        station_color = kwargs.get("station_color", "c")
        event_color = kwargs.get("event_color", "orange")

        f, ax = plt.subplots(figsize=(8, 8))

        df = self.misfit(level="station").loc[iteration, step_count]

        # Get lat/lon information from sources and receivers
        stations = self.receivers.droplevel(0)  # remove network index
        events = self.sources.drop(["time", "magnitude", "depth_km"], axis=1)

        plotted, names = [], []
        for event, sta in df.index.to_numpy():
            elon, elat = events.loc[event].longitude, events.loc[event].latitude
            slon, slat = stations.loc[sta].longitude, stations.loc[sta].latitude
            # Plot a marker for each event and station
            if event not in plotted:
                plt.scatter(elon, elat, marker="o", c=event_color,
                            edgecolors="k", s=25, zorder=100)
                plotted.append(event)
            if sta not in plotted:
                plt.scatter(slon, slat, marker="v", c=station_color,
                            edgecolors="k", s=25, zorder=100)
                plotted.append(event)

            # Connect source and receiver with a line
            plt.plot([elon, slon], [elat, slat], color=ray_color, linestyle="-",
                     alpha=0.1, zorder=50)

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title(f"Raypaths ({len(events)} events, {len(stations)} stations)")

        if save:
            plt.savefig(save)
        if show:
            plt.show()
        else:
            plt.close()

    def measurement_hist(self, iteration, step_count, choice="event", show=True,
                         save=False):
        """
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
        """
        arr = self.nwin(
            level=choice).loc[iteration, step_count].nwin.to_numpy()

        n, bins, patches = plt.hist(x=arr, color="orange", histtype="bar",
                                    edgecolor="black", linewidth=4.,
                                    label=choice, alpha=1., zorder=20
                                    )
        mu, var, std = get_histogram_stats(n, bins)
        plt.axvline(x=mu, ymin=0, ymax=1, linewidth=2, c="k",
                    linestyle="--", zorder=15, alpha=0.5)
        for sign in [-1, 1]:
            plt.axvline(x=mu + sign * std, ymin=0, ymax=1, linewidth=2,
                        c="k", linestyle=":", zorder=15, alpha=0.5)

        if save:
            plt.savefig(save)
        if show:
            plt.show()
        else:
            plt.close()

    def station_misfit_map(self, station, iteration, step_count, choice,
                           show=True, save=False, **kwargs):
        """
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
        """
        assert (station in self.stations), "station name not found"
        cmap = kwargs.get("cmap", "viridis")

        sta = self.receivers.droplevel(0).loc[station]

        # Get misfit on a per-station basis 
        df = self.misfit(level="station").loc[
            iteration, step_count].swaplevel(0, 1)
        df = df.sort_values(by="station").loc[station]

        # Get source lat/lon values as a single dataframe with same index name
        src = self.sources.drop(["time", "magnitude", "depth_km"], axis=1)
        src.index.names = ["event"]

        # This is a dataframe of events corresponding to a single station
        df = df.merge(src, on="event")

        f, ax = plt.subplots()
        src = plt.scatter(sta.longitude, sta.latitude, marker="v", c="orange",
                          edgecolors="k", s=25, zorder=100)
        plt.scatter(df.longitude.to_numpy(), df.latitude.to_numpy(),
                    c=df[choice].to_numpy(), marker="o", s=25, zorder=99,
                    cmap=cmap)

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title(f"{station} {iteration}{step_count}; {len(df)} events")

        colormap_colorbar(cmap, vmin=df[choice].to_numpy().min(),
                          vmax=df[choice].to_numpy().max(), cbar_label=choice,
                          )

        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, src, df.index.to_numpy(), **kwargs)
            plt.show()

        return f, ax

    def event_misfit_map(self, event, iteration, step_count, choice, show=True,
                         save=False, **kwargs):
        """
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
        """
        assert (event in self.sources.index), "event name not found"
        cmap = kwargs.get("cmap", "viridis")

        f, ax = plt.subplots()
        source = self.sources.loc[event]
        src = plt.scatter(source.longitude, source.latitude, marker="o", c="r",
                          edgecolors="k", s=20, zorder=100)

        # Go through each of the stations corresponding to this source
        df = self.misfit(level="station").loc[iteration, step_count, event]
        assert (choice in df.columns), f"choice must be in {df.columns}"

        # Get lat lon values for receivers
        df = df.merge(self.receivers, on="station")
        misfit_values = df[choice].to_numpy()
        rcvs = plt.scatter(df.longitude.to_numpy(), df.latitude.to_numpy(),
                           c=misfit_values, marker="v", s=15, zorder=100,
                           cmap=cmap
                           )

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title(f"{event} {iteration}{step_count}; {len(df)} stations")

        colormap_colorbar(cmap, vmin=misfit_values.min(),
                          vmax=misfit_values.max(), cbar_label=choice,
                          )

        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, rcvs, df.index.to_numpy(), **kwargs)
            plt.show()

        return f, ax

    def hist(self, iteration=None, step_count=None, iteration_comp=None,
             step_count_comp=None, f=None, ax=None, event=None, station=None,
             choice="cc_shift_in_seconds", binsize=None, show=True, save=None,
             **kwargs):
        """
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
        """
        # Optional kwargs for fine tuning figure parameters
        title = kwargs.get("title", "")
        xlim = kwargs.get("xlim", None)
        color = kwargs.get("color", "darkorange")
        color_comp = kwargs.get("color_comp", "deepskyblue")
        fontsize = kwargs.get("fontsize", 12)
        figsize = kwargs.get("figsize", (8, 6))
        legend = kwargs.get("legend", True)
        legend_loc = kwargs.get("legend_loc", "best")
        label_range = kwargs.get("label_range", False)
        xstep = kwargs.get("xstep", 2)
        ymax = kwargs.get("ymax", None)
        xlabel = kwargs.get("xlabel", None)
        zeroline = kwargs.get("zeroline", False)
        meanline = kwargs.get("meanline", False)
        stdline = kwargs.get("stdline", False)
        linewidth = kwargs.get("linewidth", 2.5)
        label = kwargs.get("label", None)
        label_comp = kwargs.get("label_comp", None)

        # If no arguments are given, default to first and last evaluations
        if iteration is None and iteration_comp is None:
            iteration, step_count = self.initial_model
            iteration_comp, step_count_comp = self.final_model

        # Check that the provided values are available in the Inspector
        assert iteration in self.iterations, \
            f"iteration must be in {self.iterations}"
        if step_count is None:
            assert step_count in self.steps.loc[iteration], \
                f"step must be in {self.steps.loc[iteration]}"
        if iteration_comp is not None:
            assert iteration_comp in self.iterations, \
                f"iteration_comp must be in {self.iterations}"
            assert step_count_comp in self.steps.loc[iteration_comp], \
                f"step_comp must be in {self.steps.loc[iteration_comp]}"

        # Try to set a default binsize that may or may not work 
        if binsize is None:
            try:
                binsize = {"cc_shift_in_seconds": 1,
                           "dlnA": 0.25,
                           "max_cc_value": 0.05,
                           "misfit": 10,
                           "relative_starttime": 15,
                           "relative_endtime": 15}[choice]
            except KeyError:
                binsize = 1

        def get_values(m, s, e, sta):
            """short hand to get the data, and the maximum value in DataFrame"""
            df_a = self.isolate(iteration=m, step_count=s, event=e, station=sta)
            try:
                val_ = df_a.loc[:, choice].to_numpy()
            except KeyError as e:
                raise KeyError(f"Inspector.windows has no key {choice}") from e
            lim_ = max(abs(np.floor(min(val_))), abs(np.ceil(max(val_))))
            return val_, lim_

        # Instantiate the plot objects and 'goforyourlifemate'
        if f is None:
            f, ax = plt.subplots(figsize=figsize)
        if ax is None:
            ax = plt.gca()

        val, lim = get_values(iteration, step_count, event, station)

        if iteration_comp:
            val_comp, lim_comp = get_values(iteration_comp, step_count_comp,
                                            event, station)

            # Reset the limit to be the greater of the two
            lim = max(lim, lim_comp)

            # Compare iterations, plot original iteration on top 
            n, bins, patches = plt.hist(
                x=val, bins=np.arange(-1 * lim, lim + .1, binsize),
                color=color, histtype="bar", edgecolor="black",
                linewidth=linewidth,
                label=(label or f"{iteration}{step_count}") + f"; N={len(val)}",
                zorder=11, alpha=1.
            )

            # mu1, var1, std1 = get_histogram_stats(n, bins)
            mean = np.mean(val)
            std = np.std(val)
            med = np.median(val)

            # Plot comparison below
            n2, bins2, patches2 = plt.hist(
                x=val_comp, bins=np.arange(-1 * lim, lim + .1, binsize),
                color=color_comp, histtype="bar", edgecolor="black",
                linewidth=linewidth + 1, zorder=10,
                label=(label_comp or f"{iteration_comp}{step_count_comp}") +
                      f"; N={len(val_comp)}",
            )
            # mu2, var2, std2 = get_histogram_stats(n2, bins2)
            mean_comp = np.mean(val_comp)
            std_comp = np.std(val_comp)
            med_comp = np.median(val_comp)

            # Plot edges of comparison over top
            plt.hist(x=val_comp, bins=np.arange(-1 * lim, lim + .1, binsize),
                     color="k", histtype="step", edgecolor=color_comp,
                     linewidth=linewidth, zorder=12,
                     )
        else:
            # No comparison iteration, plot single histogram
            n, bins, patches = plt.hist(
                x=val, bins=len(np.arange(-1 * lim, lim, binsize)), color=color,
                histtype="bar", edgecolor="black", linewidth=linewidth,
                label=f"N={len(val)}", zorder=10,
            )
            # mu1, var1, std1 = get_histogram_stats(n, bins)
            mean = np.mean(val)
            std = np.std(val)
            med = np.median(val)

        # Plot reference lines
        if zeroline:
            plt.axvline(x=0, ymin=0, ymax=1, linewidth=2.5, c="k",
                        zorder=15, alpha=0.75, linestyle="-")

        # Plot the mean of the histogram
        if meanline:
            plt.axvline(x=mean, ymin=0, ymax=1, linewidth=2.5, c="k",
                        linestyle="--", zorder=15, alpha=0.75)
        # Plot one standard deviation
        if stdline:
            for sign in [-1, 1]:
                plt.axvline(x=mean + sign * std, ymin=0, ymax=1, linewidth=2.5,
                            c="k", linestyle=(0, (1, 1)), zorder=15, alpha=0.75)

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
            tit_fmt = "mean: {mean:.2f} / std: {std:.2f} / med: {med:.2f}"
            title = tit_fmt.format(mean=mean, std=std, med=med)
            if iteration_comp:
                tit_comp = tit_fmt.format(mean=mean_comp, std=std_comp,
                                          med=med_comp)
                title = " ".join([f"[{label or iteration}]", title, "\n",
                                  f"[{label_comp or iteration_comp}]", tit_comp
                                  ])

        # Finalize plot details
        if xlabel:
            xlab_ = xlabel
        else:
            try:
                # For cleaner formatting of x-axis label
                xlab_ = {"cc_shift_in_seconds": "Time Shift (s)",
                         "dlnA": "$\Delta\ln$(A)",
                         "misfit": "Misfit",
                         "length_s": "Window Length (s)",
                         "max_cc_value": "Peak Cross Correlation",
                         "relative_starttime": "Relative Start Time (s)",
                         "relative_endtime": "Relative End Time (s)",
                         }[choice]
            except KeyError:
                xlab_ = choice

        plt.xlabel(xlab_, fontsize=fontsize)
        plt.ylabel("Count", fontsize=fontsize)
        plt.title(title, fontsize=fontsize)
        plt.tick_params(which='both', direction='in', top=True, right=True,
                        labelsize=fontsize, width=2.)
        if label_range:
            plt.xticks(np.arange(-1 * label_range, label_range + .1, step=xstep))

        if legend:
            leg = plt.legend(fontsize=fontsize / 1.25, loc=legend_loc)
            # Thin border around legend objects, unnecessarily thick bois
            for leg_ in leg.legendHandles:
                leg_.set_linewidth(1.5)

        plt.tight_layout()

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax

    def plot_windows(self, iteration, step_count, event=None, network=None,
                     station=None, component=None):
        """
        Show lengths of windows chosen based on source-receiver distance, akin
        to Tape's Thesis or to the LASIF plots. These are useful for showing
        which phases are chosen, and window choosing behavior as distance
        increases and (probably) misfit increases.

        :type iteration: str
        :param iteration: iteration to analyze
        :type step_count: str
        :param step_count: step count to query, e.g. 's00'
        :type event: str
        :param event: filter for measurements for a given event
        :type network: str
        :param network: filter for measurements for a given network
        :type station: str
        :param station: filter for measurements for a given station
        :type component: str
        :param component: choose a specific component to analyze
        """
        assert(iteration in self.iterations and
               step_count in self.steps[iteration]), \
            f"{iteration}{step_count} does not exist in Inspector"

        comp_dict = {"Z": "orangered",
                     "N": "forestgreen", "R": "forestgreen",
                     "E": "royalblue", "T": "royalblue"
                     }

        srcrcv = self.calculate_srcrcv()
        windows = self.isolate(iteration=iteration, step_count=step_count,
                               event=event, network=network, station=station
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
        for comp in sorted(windows.component.unique()):
            ax.hlines(y=0, xmin=0, xmax=0.1, color=comp_dict[comp], label=comp)

        plt.title(f"{len(df)} misfit windows")
        plt.xlabel("Time [s]")
        plt.ylabel("Distance [km]")
        plt.legend()
        plt.xlim([min(df.relative_starttime) - 10,
                  max(df.relative_endtime) + 10])
        plt.ylim([min(df.distance_km) - 10, max(df.distance_km) + 10])

        plt.show()

    def plot_window_differences(self, iteration_a, step_a, iteration_b, step_b):
        """
        """
        srcrcv = self.calculate_srcrcv()
        windows_a = self.isolate(iteration=iteration_a, step_count=step_a)
        windows_b = self.isolate(iteration=iteration_b, step_count=step_b)
        df_a = windows_a.loc[:, ["event", "network", "station", "component",
                                 "relative_starttime", "relative_endtime"]]
        df_a.rename(columns={"relative_starttime": "start_a",
                             "relative_endtime": "end_a"}, inplace=True)
        df_b = windows_b.loc[:, ["event", "network", "station", "component",
                                 "relative_starttime", "relative_endtime"]]
        df_b.rename(columns={"relative_starttime": "start_b",
                             "relative_endtime": "end_b"}, inplace=True)

        # Comparison dataframe
        df_c = df_a.merge(df_b)
        df_c = df_c.merge(srcrcv.drop("backazimuth", axis=1),
                          on=["event", "network", "station"])
        df_c.drop(["event", "network", "station"], axis=1, inplace=True)

        f, ax = plt.subplots(figsize=(8, 6))
        for window in df_c.to_numpy()[:100]:
            comp, start_a, end_a, start_b, end_b, dist = window
            # short time windows show on top
            if start_a - start_b:
                if start_a < start_b:
                    c = "orangered"
                else:
                    c = "forestgreen"
                ax.hlines(y=dist, xmin=start_a, xmax=start_b, colors=c,
                          alpha=0.4)
            if end_a - end_b:
                if end_a < end_b:
                    c = "orangered"
                else:
                    c = "forestgreen"
                ax.hlines(y=dist, xmin=end_a, xmax=end_b, colors=c, alpha=0.4)
            ax.hlines(y=dist, xmin=0, xmax=300, colors="k", alpha=0.1,
                      linewidth=0.1
                      )

    def convergence(self, windows="length_s", trials=False, show=True,
                    save=None, normalize=False, float_tolerance=1E-3,
                    annotate=False, restarts=None, **kwargs):
        """
        Plot the convergence rate over the course of an inversion.
        Scatter plot of total misfit against iteration number, or by step count

        .. note:: 
            Because misfits are floats, they wont be exactly equal, so we need 
            to set some small tolerance in which they can differ

        :type windows: str or bool
        :param windows: parameter to use for Inspector.measurements() to
            determine how to illustrate measurement number, either by
            length_s: cumulative window length in seconds
            nwin: number of misfit windows
            None: will not plot window information
        :type trials: bool
        :param trials: plot the discarded trial step function evaluations from
            the line searches. Useful for understanding how efficient the
            optimization algorithm as
        :type normalize: bool
        :param normalize: normalize the objective function values between [0, 1]
        :type float_tolerance: float
        :param float_tolerance: acceptable floating point difference between
            adjacent misfit values to say that they are equal
        :type restarts: list of int
        :param restarts: If the inversion was restarted, e.g. for parameter
            changes, then the convergence figure should separate two line plots.
            This list allows the User to tell the function where to separate
            the convergence plot. The integers should correspond to indices of
            the Inspector.models attribute.
        :type annotate: bool
        :param annotate: annotate misfit values next to markers
        :type show: bool
        :param show: show the plot after making it
        :type save: str
        :param save: file id to save the figure to
        """
        f = kwargs.get("f", None)
        ax = kwargs.get("ax", None)
        fontsize = kwargs.get("fontsize", 15)
        figsize = kwargs.get("figsize", (8, 6))
        legend = kwargs.get("legend", True)
        misfit_label = kwargs.get("misfit_label", "misfit")
        trial_label = kwargs.get("trial_label", "trials")
        window_label = kwargs.get("window_label", "windows")
        trial_color = kwargs.get("trial_color", "r")
        window_color = kwargs.get("window_color", "orange")
        legend_loc = kwargs.get("legend_loc", "best")

        # Set some default parameters based on user choices
        if windows:
            assert (windows in ["nwin", "length_s"]), \
                "plot_windows must be: 'nwin; or 'length_s'"

        # It may take a while to calculate models so do it once here
        models = self.models
        nwin = self.nwin()

        # Set up the figure
        if f is None:
            f, ax = plt.subplots(figsize=figsize)
        if ax is None:
            ax = plt.gca()

        # First, we will sort the model values by accepted models, initial
        # evaluations, and discarded trials steps. Also need to check if
        # accepted models and initial evaluations are equal to one another.
        x = 0  # the x-position on the axis
        lines, xvals, yvals, xlabs = [], [], [], []  # main plot
        xdiscards, ydiscards, ywindows, xrestarts = [], [], [], []  # secondary
        for j in range(len(models)):
            i = j - 1  # we always need to compare to the previous misfit value

            # Status 0 means initial evaluation of iteration
            if models.state[j] == 0:
                # Ignore very first function evaluation
                if j == 0:
                    pass
                # If initial eval matches line search final, treat equally
                elif abs(models.misfit[i] - models.misfit[j]) < float_tolerance:
                    continue
                # If they differ, treat them as different points
                else:
                    x += 1
                xlab = f"{models.model[j]}_0"
            # Status 1 means final evaluation in line search
            elif models.state[j] == 1:
                x += 1
                xlab = f"{models.model[j]}_1"
            # Status -1 means discarded trial step, plot on the same X value
            elif models.state[j] == -1:
                xdiscards.append(x + 1)  # discards are related to next model
                ydiscards.append(models.misfit[j])
                continue

            xvals.append(x)
            yvals.append(models.misfit[j])
            xlabs.append(xlab)

            # Convert restart values from Inspector.models indices
            if restarts and j in restarts:
                xrestarts.append(x)

            # Get the corresponding window number based on iter/step count
            if windows:
                i_ = models.iteration[j]
                s_ = models.step_count[j]
                ywindows.append(nwin.loc[i_].loc[s_][windows])

                # !!! HACK REMOVE THIS
                # try:
                #     ywindows.append(nwin.loc[i_].loc[s_][windows])
                # except KeyError:
                #     ywindows.append(nwin.loc[i_].loc["s03"][windows])

        # Define a re-usable plotting function that takes arguments from main fx
        def plot_vals(x_, y_, idx=None, c="k", label=misfit_label):
            """
            Re-used plotting commands plot a scatter plot with a certain
            color and label. Normalizes y-values, annotates text, if required

            :type x_: np.array
            :param x_: x values to plot
            :type y_: np.array
            :param y_: y values to plot
            :type idx: int
            :param idx: index of the inversion leg for color and label, if None
                defaults to `c` and `label` for color and label
            :type c: str
            :param c: color for marker and line color
            :type label: str
            :param label: label for legend, defaults to `misfit_label` from
                kwargs of main function
            """
            # Overwrite default values 
            if idx is not None:
                c = f"C{idx}"
                label = f"{misfit_label} (leg {idx})"

            if normalize:
                y_ = [_ / max(y_) for _ in y_]

            line = ax.plot(x_, y_, "o-", linewidth=3,  markersize=10, 
                           c=c, label=label, zorder=10)
            if annotate:
                for x_anno, y_anno in zip(x_, y_):
                    plt.text(x_anno, y_anno, f"{y_anno:.3f}", zorder=11)

            return line

        # Primary: Two methods of plotting:
        if xrestarts:
            # 1) with user-defined restarts separating legs of the inversion
            first = 0  # first iteration in the current leg
            for i, last in enumerate(xrestarts):
                j = i + 1  # Leg counting should start at 1
                lines += plot_vals(xvals[first:last], yvals[first:last], j)
                first = last
            # Plot the final leg
            lines += plot_vals(xvals[last:], yvals[last:], j + 1)
        else:
            # 2) plot the entire convergence in one line
            lines += plot_vals(xvals, yvals, idx=None)

        # Secondary: Plot number of windows/ window length in a separate axis
        if windows:
            ax2 = ax.twinx()
            # Set ax2 below ax1
            ax.set_zorder(ax2.get_zorder() + 1)
            ax.patch.set_visible(False)
            lines += ax2.plot(xvals, ywindows, "v:", linewidth=2, markersize=8,
                              c=window_color, label=window_label, zorder=5
                              )
            ydict = {"length_s": "Cumulative Window Length [s]",
                     "nwin": "Numer of Measurements"}
            ax2.set_ylabel(f"{ydict[windows]} (dashed)", rotation=270,
                           labelpad=15., fontsize=fontsize
                           )
            ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

        # Secondary: Plot the discarded trial steps as x's
        if trials:
            sc = ax.scatter(xdiscards, ydiscards, c=trial_color,
                            marker="x", s=10, zorder=9, label=trial_label
                            )
            lines.append(sc)

        # Format the axes
        ax.set_xlabel("Model Number", fontsize=fontsize)
        ax.xaxis.set_ticks(xvals)
        ax.set_xticklabels(xlabs, rotation=45, ha="right")
        ax.set_ylabel("Total Normalized Misfit", fontsize=fontsize)

        # Only set ticks on the x-axis
        ax.xaxis.grid(True, which="minor", linestyle=":")
        ax.xaxis.grid(True, which="major", linestyle="-")

        if legend:
            labels = [line.get_label() for line in lines]
            ax.legend(lines, labels, prop={"size": 12}, loc=legend_loc)

        f.tight_layout()

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax


def colormap_colorbar(cmap, vmin=0., vmax=1., dv=None, cbar_label=""):
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
    :type cbar_label: str
    :param cbar_label: label for colorbar
    :rtype:
    :return:
    """
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    sm.set_clim(vmin, vmax)
    if dv:
        boundaries = np.arange(vmin, vmax, dv)
    else:
        boundaries = None
    cbar = plt.colorbar(sm, boundaries=boundaries, shrink=0.9, pad=0.025)
    if cbar_label:
        cbar.ax.set_ylabel(cbar_label, rotation=270, labelpad=15)

    return cmap, norm, cbar


def hover_on_plot(f, ax, obj, values, dissapear=True):
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


def get_histogram_stats(n, bins):
    """
    Get mean, variance and standard deviation from a histogram
    
    :type n: array or list of arrays
    :param n: values of histogram bins
    :type bins: array
    :param bins: edges of the bins
    """
    mids = 0.5 * (bins[1:] + bins[:-1])
    mean = np.average(mids, weights=n)
    var = np.average((mids - mean) ** 2, weights=n)
    std = np.sqrt(var)

    return mean, var, std


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
    acceptable_locations = ["lower-right", "upper-right",
                            "lower-left", "upper-left"]
    assert(anno_location in acceptable_locations), \
        f"anno_location must be in {acceptable_locations}"

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
