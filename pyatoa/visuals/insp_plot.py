#!/usr/bin/env python3
"""
The plotting functionality of the Inspector class. Used to generate statistics
and basemap like plots from the Inspector DataFrame objects.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa import logger
from matplotlib.patches import Rectangle
from pyatoa.utils.calculate import normalize_a_to_b


# A map from the Pyflex parameter names into cleaner looking label strings
common_labels = {"cc_shift_in_seconds": "Time Shift (s)",
                 "dlnA": "$\Delta\ln$(A)",
                 "misfit": "Misfit",
                 "length_s": "Window Length (s)",
                 "max_cc_value": "Peak Cross Correlation",
                 "relative_starttime": "Relative Start Time (s)",
                 "relative_endtime": "Relative End Time (s)",
                 }


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
                srcs = self.sources.loc[event]
            else:
                srcs = self.sources

            src_lat = srcs.latitude
            src_lon = srcs.longitude
            src_names = srcs.index
            sc_sources = plt.scatter(src_lon, src_lat, marker="o", c="None",
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
                sc_receivers = plt.scatter(rcv_lon, rcv_lat, marker="v",
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
            hover_on_plot(f, ax, sc_sources, src_names)
            for rcvs, rcv_names in zip(sc_receiver_list, sc_receiver_names):
                hover_on_plot(f, ax, rcvs, rcv_names)
            plt.show()

        return f, ax

    def scatter(self, x, y, iteration=None, step_count=None, save=None,
                show=True, **kwargs):
        """
        Create a scatter plot between two chosen keys in the windows attribute

        :type x: str
        :param x: key to choose for the x axis of the plot
        :type y: str
        :param y: key to chooose for the y axis of the plot
        :type iteration: str
        :param iteration: the chosen iteration to plot for, if None will default
            to the latest iteration available
        :type step_coutn: str
        :param step_count: chosen step count. If None, defaults to latest
        """
        if iteration is None:
            iteration, _ = self.initial_model
        if step_count is None:
            step_count = self.steps.loc[iteration][-1]

        # Ensure we have distance and backazimuth values in the dataframe
        df = self.isolate(iteration=iteration, step_count=step_count, **kwargs)
        df = df.merge(self.srcrcv, on=["event", "network", "station"])

        assert(x in df.keys()), f"X value {x} does not match keys {df.keys()}"
        assert(y in df.keys()), f"Y value {y} does not match keys {df.keys()}"

        f, ax = plt.subplots(figsize=(8, 6))
        plt.scatter(df[x].to_numpy(), df[y].to_numpy(), **kwargs)
        default_axes(ax, **kwargs)

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax

    def travel_times(self, iteration=None, step_count=None, component=None,
                     constants=None, t_offset=0, hist=False, hist_max=None, 
                     save=None, show=True, **kwargs):
        """
        Plot relative window starttime (proxy for phase arrival) against 
        source-receiver distance, to try to convey which phases are included
        in the measurement. 
        
        Similar to Figure 4.18 in Shearer's Intro to Seismology.

        :type iteration: str
        :param iteration: the chosen iteration to plot for, if None will default
            to the latest iteration available
        :type step_count: str
        :param step_count: chosen step count. If None, defaults to latest
        :type component: str
        :param component: optional specify a measurement component to isolate 
            only e.g., 'Z' components to look at Rayleigh waves
        :type constants: list of floats
        :param constants: plot lines of constant velocity to estimate the 
            average wavespeed that leads to some of the linear trends
        :type t_offset: float
        :param t_offset: if the synthetic offset time in SPECFEM is set then 
            the constant lines will need to be offset by the same amount to 
            match the measurements.
        :type hist: bool
        :param hist: create a histogram binning the approximate seismic 
            velocities
        """
        hist_color = kwargs.get("hist_color", "deepskyblue")
        title_plot = kwargs.get("title_plot", None)
        title_hist = kwargs.get("title_hist", None)

        if iteration is None:
            iteration, _ = self.final_model
        if step_count is None:
            step_count = self.steps.loc[iteration][-1]

        # Ensure we have distance and backazimuth values in the dataframe
        df = self.isolate(iteration=iteration, step_count=step_count,
                          component=component)
        df = df.merge(self.srcrcv, on=["event", "network", "station"])

        # Assuming that isolate has only picked values from a single iterstep
        iterstep = f"{df.iteration[0]}{df.step[0]}"

        dist, start, length = df[["distance_km", "relative_starttime", 
                                  "length_s"]].to_numpy().T

        # Shift relative starttimes by the user-defined offset
        start -= t_offset

        # size of the markers based on the length of the window
        length = normalize_a_to_b(length, .5, .5)

        f, ax = plt.subplots(figsize=(8, 6))

        plt.scatter(dist, start, c="k", s=.25, marker="x", zorder=5)
        if title_plot is not None:
            plt.title(title_plot)
        else:
            plt.title(f"Apparent travel times ({iterstep} N={len(dist)})")
        plt.xlabel("Source-receiver distance [km]")
        plt.ylabel("Relative start time [s]")

        if constants is not None:
            x = np.linspace(0, dist.max(), len(dist))
            for i, c in enumerate(constants):
                y = x / c
                plt.plot(x, y, c=f"C{i}", lw=1, zorder=1, label=f"{c} km/s")
            plt.legend()

        plt.xlim([0, dist.max()])
        plt.ylim([0, start.max()])
        f.tight_layout()
        default_axes(ax, **kwargs)

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        plt.close()

        # Now make a separate histogram showing the apparent velocities
        if hist:
            velocities = dist / start
            # Max velocity based on PREM highest (ish) Vp
            n, bins, patches = plt.hist(x=velocities, 
                                        bins=np.arange(0, 12, .5), 
                                        color=hist_color, histtype="bar", 
                                        edgecolor="black", linewidth=2, 
                                        zorder=11, alpha=1.
                                        )
            if hist_max:
                plt.ylim([0, hist_max])

            if title_hist is not None:
                plt.title(title_hist)
            else:
                plt.title(f"Apparent velocities "
                          f"({iterstep} N={len(velocities)})")
            plt.xlabel("Velocity [km/s]")
            plt.ylabel("Count")
            plt.gcf().tight_layout()
            default_axes(plt.gca(), **kwargs)

            if save:
                plt.savefig(f"hist_{save}")
            if show:
                plt.show()


    def event_depths(self, xaxis="longitude", show=True, save=None, **kwargs):
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

        default_axes(ax, **kwargs)

        if save:
            plt.savefig(save)
        if show:
            plt.show

        return f, ax

    def raypaths(self, iteration, step_count, color_by=None, show=True, 
                 save=False, vmin=None, vmax=None, **kwargs):
        """
        Plot rays connecting sources and receivers based on the availability
        of measurements. Useful for getting an approximation of resolution.

        :type iteration: int
        :param iteration: iteration to retrieve data from
        :type step_count: int
        :param step_count: step count to retrieve data from
        :type color_by: str
        :param color_by: allow rays to be colored based on a normalized value.
            nwin: color rays by the number of windows available for that path
            misfit: color rays by total misfit
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        cmap = kwargs.get("cmap", "viridis")
        ray_color = kwargs.get("ray_color", "k")
        ray_linewidth = kwargs.get("ray_linewidth", 1)
        ray_alpha = kwargs.get("ray_alpha", 0.1)
        station_color = kwargs.get("station_color", "c")
        event_color = kwargs.get("event_color", "orange")
        figsize = kwargs.get("figsize", (8, 8))
        markersize = kwargs.get("markersize", 25)

        f, ax = plt.subplots(figsize=figsize)

        df = self.misfit(level="station").loc[iteration, step_count]

        # Get lat/lon information from sources and receivers
        stations = self.receivers.droplevel(0)  # remove network index
        events = self.sources.drop(["time", "magnitude", "depth_km"], axis=1)

        # Set up the normalized colorbar 
        cbar, extend = None, None
        if color_by is not None:
            assert(color_by in df.keys()), f"{color_by} must be in {df.keys()}"
            if vmin is None:
                vmin = df[color_by].min()
            elif vmin > df[color_by].min():
                extend = "min"
            if vmax is None:
                vmax = df[color_by].max()
            elif vmax < df[color_by].max():
                if extend == "min":
                    extend = "both"
                else:
                    extend = "max"

            sm, norm, cbar = colormap_colorbar(cmap, vmin=vmin, vmax=vmax,
                                               cbar_label=color_by.capitalize(),
                                               extend=extend
                                               )

        plotted, names = [], []
        for event, sta in df.index.to_numpy():
            elon, elat = events.loc[event].longitude, events.loc[event].latitude
            slon, slat = stations.loc[sta].longitude, stations.loc[sta].latitude
            # Plot a marker for each event and station
            if event not in plotted:
                plt.scatter(elon, elat, marker="o", c=event_color,
                            edgecolors="k", s=markersize, zorder=100)
                plotted.append(event)
            if sta not in plotted:
                plt.scatter(slon, slat, marker="v", c=station_color,
                            edgecolors="k", s=markersize, zorder=100)
                plotted.append(event)

            if color_by is not None:
                ray_color = sm.cmap(norm(df.loc[event].loc[sta][color_by]))

            # Connect source and receiver with a line
            plt.plot([elon, slon], [elat, slat], color=ray_color, linestyle="-",
                     alpha=0.1, zorder=50, linewidth=ray_linewidth)

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title(f"{len(df)} raypaths")
        # plt.title(f"{len(df)} raypaths ({len(events)} events, "
        #           f"{len(stations)} stations)")

        default_axes(ax, cbar, **kwargs)

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax

    def raypath_density(self, iteration, step_count, point_spacing_km=.5,
                        bin_spacing_km=8, cmap="viridis", show=True, save=False, 
                        **kwargs):
        """
        Create a raypath density plot to provide a more deatiled illustration of 
        raypath gradients, which may be interpreted alongside tomographic 
        inversion results as a preliminary resolution test.

        The idea behind this is to partition each individual raypath line into
        discrete points and then create a 2D histogram with all points

        :type point_spacing_km: float
        :param point_spacing_km: approximate discretization interval for each
            raypath line. Smaller numbers will lead to higher resolution but
            also longer computation time.
        :type bin_spacing_km: float
        :param bin_spacing_km: the bin size in km of the 2d histogram. If
            the same as 'point_spacing_km' then you'll probably just see the
            lines. Should be larger than 'point_spacing_km' for a more
            contour plot looking feel.
        """
        figsize = kwargs.get("figsize", (8, 8))

        f, ax = plt.subplots(figsize=figsize)
        df = self.misfit(level="station").loc[iteration, step_count]

        # Get lat/lon information from sources and receivers
        stations = self.receivers.droplevel(0)  # remove network index
        events = self.sources.drop(["time", "magnitude", "depth_km"], axis=1)

        # Determine grid bounds and required number of bins for histograms
        x_min = min(stations.longitude.min(), events.longitude.min())
        x_max = max(stations.longitude.max(), events.longitude.max())
        y_min = min(stations.latitude.min(), events.latitude.min())
        y_max = max(stations.latitude.max(), events.latitude.max())

        # 111.11 VERY roughly converts degrees to km, not really geographically
        # correct though. Should be okay for this low-res application
        x_bins = int(abs(x_max - x_min)  * 111.11 / bin_spacing_km)
        y_bins = int(abs(y_max - y_min)  * 111.11 / bin_spacing_km)

        # Convert station names and event ids into coordinates
        dx = point_spacing_km / 111.11  # grid spacing in degrees

        # Initiate empty arrays to be filled
        x = np.array([])
        y = np.array([])
        plotted = []
        for event, sta in df.index.to_numpy():
            elon, elat = events.loc[event].longitude, events.loc[event].latitude
            slon, slat = stations.loc[sta].longitude, stations.loc[sta].latitude

            # Plot a marker for each event and station
            if event not in plotted:
                plt.scatter(elon, elat, marker="o", c=event_color,
                            edgecolors="k", s=markersize, zorder=100)
                plotted.append(event)
            if sta not in plotted:
                plt.scatter(slon, slat, marker="v", c=station_color,
                            edgecolors="k", s=markersize, zorder=100)
                plotted.append(event)
           
            # Calculate the necessary number of discrete points to create line
            nlon = int(abs(elon - slon) * 111.11 / point_spacing_km)
            nlat = int(abs(elat - slat) * 111.11 / point_spacing_km)
            nvals = max(nlon, nlat)

            x_ = np.linspace(elon, slon, nvals)
            y_ = np.linspace(elat, slat, nvals)

            x = np.concatenate((x, x_))
            y = np.concatenate((y, y_))

        # Create the 2D histogram of raypath density
        plt.hist2d(x, y, bins=(x_bins, y_bins), cmap=plt.get_cmap(cmap), 
                   zorder=5)
        cbar = plt.colorbar(label="counts", shrink=0.9, pad=0.025)

        plt.scatter(coast[:,1], coast[:,0], s=.05, c="w", zorder=20)  # DELETE
        plt.title(f"Raypath Density (N={len(df)} src-rcv pairs)")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")

        # Calculate aspect ratio based on latitude
        w = 1 / np.cos(np.radians(y[0]))
        plt.gca().set_aspect(w)

        default_axes(plt.gca(), cbar)

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        plt.close()





    def event_hist(self, choice):
        """
        Make a histogram of event information
        :return:
        """
        assert choice in self.sources.keys(), \
            f"Choice must be in {self.sources.keys()}"

        f, ax = plt.subplots()
        arr = self.sources[choice].to_numpy()

        # Compare iterations, plot original iteration on top
        n, bins, patches = plt.hist(x=arr, color="w", histtype="bar",
                                    bins=list(np.arange(4.5, 6.1, .1)),
                                    edgecolor="black", linewidth=2.,
                                    label=choice, alpha=1., zorder=20
                                    )
        mu, var, std = get_histogram_stats(n, bins)
        default_axes(ax)


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

    def station_event_misfit_map(self, station, iteration, step_count, choice,
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

        _, _, cbar = colormap_colorbar(cmap, vmin=df[choice].to_numpy().min(),
                                       vmax=df[choice].to_numpy().max(), 
                                       cbar_label=choice)

        default_axes(ax, cbar, **kwargs)

        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, src, df.index.to_numpy())
            plt.show()

        return f, ax

    def event_station_misfit_map(self, event, iteration, step_count, choice,
                                 show=True, save=False, **kwargs):
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

        _, _, cbar = colormap_colorbar(cmap, vmin=misfit_values.min(),
                                       vmax=misfit_values.max(), 
                                       cbar_label=choice,)

        default_axes(ax, cbar, **kwargs)

        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, rcvs, df.index.to_numpy())
            plt.show()

        return f, ax

    def event_misfit_map(self, choice=None, iteration=None, step_count=None,
                         show=True, save=False, **kwargs):
        """
        Plot all events on a map and their corresponding scaled misfit value

        :type iteration: str
        :param iteration: iteration number e.g. 'i00'
        :type step_count: str
        :param step_count: step count e.g. 's00'
        :type choice: str
        :param choice: choice of misfit value, either 'misfit' or 'nwin' or
            'unscaled_misfit'
        :type show: bool
        :param show: Show the plot
        :type save: str
        :param save: fid to save the given figure
        """
        cmap = kwargs.get("cmap", "viridis")
        markersize = kwargs.get("markersize", 20)
        marker = kwargs.get("marker", "o")

        if iteration is None:
            iteration, step_count = self.final_model

        if choice is None:
            choice = "misfit"

        f, ax = plt.subplots()
        sources = self.sources.drop(["time", "magnitude", "depth_km"], axis=1)
        # Rename from event_id to event to match the naming of the dataframe
        sources.rename_axis("event", inplace=True)
        df = self.misfit(level="event").loc[iteration, step_count]
        df = df.merge(sources, on="event")

        srcs = plt.scatter(df.longitude.to_numpy(), df.latitude.to_numpy(),
                           c=df[choice].to_numpy(), marker=marker, s=markersize,
                           zorder=100, cmap=cmap)

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title(f"{iteration}{step_count}; {choice} {len(df)} events")

        _, _, cbar = colormap_colorbar(cmap, vmin=df[choice].to_numpy().min(),
                                       vmax=df[choice].to_numpy().max(),
                                       cbar_label=choice,)

        default_axes(ax, cbar, **kwargs)

        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, srcs, df.index.to_numpy())
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
                xlab_ = common_labels[choice]
            except KeyError:
                xlab_ = choice

        plt.xlabel(xlab_, fontsize=fontsize)
        plt.ylabel("Count", fontsize=fontsize)
        plt.title(title)
        if label_range:
            plt.xticks(np.arange(-1 * label_range, label_range + .1, 
                                 step=xstep))

        if legend:
            leg = plt.legend(fontsize=fontsize / 1.25, loc=legend_loc)
            # Thin border around legend objects, unnecessarily thick bois
            for leg_ in leg.legendHandles:
                leg_.set_linewidth(1.5)

        default_axes(ax, **kwargs)

        plt.tight_layout()

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax

    def plot_windows(self, iteration, step, iteration_comp=None,
                     step_comp=None, choice="cc_shift_in_seconds",
                     event=None, network=None, station=None, component=None,
                     no_overlap=True, distances=False, annotate=False,
                     bounds=False, show=True, save=False, **kwargs):
        """
        Show lengths of windows chosen based on source-receiver distance, akin
        to Tape's Thesis or to the LASIF plots. These are useful for showing
        which phases are chosen, and window choosing behavior as distance
        increases and (probably) misfit increases.

        :type iteration: str
        :param iteration: iteration to analyze
        :type step: str
        :param step: step count to query, e.g. 's00'
        :type iteration_comp: str
        :param iteration_comp: Optional, if provided, difference the 'choice'
            values with the chosen 'iteration/step'. Useful for easily checking
            for improvement. Only works if the windows are the same.
        :type step_comp: str
        :param step_comp: associated step count for 'iteration_comp'
        :type event: str
        :param event: filter for measurements for a given event
        :type network: str
        :param network: filter for measurements for a given network
        :type station: str
        :param station: filter for measurements for a given station
        :type component: str
        :param component: choose a specific component to analyze
        :type choice: str
        :param choice: choice of value to define the colorscale by. These relate
            to the keys of Inspector.windows. Default is 'cc_shift_in_seconds'
        :type no_overlap: bool
        :param no_overlap: If real distances are used, many src-rcv pairs are
            at the same or very similar distances, leading to overlapping
            rectangles. If this is set to True, to minimize overlap, the
            function will try to shift the distance to a value that hasn't yet
            been plotted. It will alternate larger positive and negative values
            until something is found. Will lead to non-real distances.
        :type distances: bool
        :param distances: If set False, just plot one window atop the other,
            which makes for more concise, easier to view plots, but
            then real distance information is lost, only relative distance
            kept.
        :type annotate: bool
        :param annotate: If True, will annotate event and station information
            for each window. May get messy if `distances == True` and
            `no_overlap == False` because you will get many overlapping
            annotations. Works ideally if `distances == False`.
        :type bounds: bool or list of float
        :param bounds:
            * (bool) False: set default bounds based on the min and max of data
            * (bool) True: set default bounds equal, based on abs max of data
            * (list) Manually set the bounds of the colorbar
        :type show: bool
        :param show: show the plot after generating
        :type save: str
        :param save: save the plot to the given filename

        Keyword Arguments
        ::
            float alpha:
                The opacity of the rectangles, defaults to 0.25
            str cmap:
                The colormap used to plot the values of `choice`
            str cbar_label:
                The label for the colorbar
            float rectangle_height:
                The vertical size of the rectangles, defaults to 1.
            float anno_shift:
                The distance in seconds to shift the plot to accomodate
                annotations. This needs to be played as its based on the length
                of the strings that are used in the annotations.
        """
        alpha = kwargs.get("alpha", 0.6)
        cmap = kwargs.get("cmap", "viridis")
        cbar_label = kwargs.get("cbar_label", None)
        rectangle_height = kwargs.get("rectangle_height", 1.0)
        anno_shift = kwargs.get("anno_shift", 50)

        assert(iteration in self.iterations and
               step in self.steps[iteration]), \
            f"{iteration}{step} does not exist in Inspector"

        assert(choice in self.windows.keys()), (f"Color by choice {choice} not "
                                                f"in list of available keys")

        # Filter out the specific windows that we're interested in
        df = self.isolate(iteration=iteration, step_count=step,
                          event=event, network=network, station=station,
                          component=component)
    
        # If a comparison iteration is given, isolate the 'choice' key, and
        # subtract it from the main dataframe. The new plotted values are diffs!
        if iteration_comp:
            df_comp = self.isolate(iteration=iteration_comp,
                                   step_count=step_comp, event=event,
                                   network=network, station=station,
                                   component=component)
            # This is enough unique info to identify a specific window
            merge_keys = ["event", "network", "station", "channel",
                          "relative_starttime", choice]
            df_comp = df_comp.loc[:, merge_keys]
            df_comp.rename({choice: f"{choice}_comp"}, axis=1, inplace=True)

            # Crude check to see if the number of windows is comparable
            assert(len(df) == len(df_comp)), (f"Number of windows does not "
                                              f"match between "
                                              f"{iteration}{step} and "
                                              f"{iteration_comp}{step_comp}")

            df = df.merge(df_comp, on=merge_keys[:-1])
            # Subtract the comparison iteration from the initial check
            df[choice] = df[choice] - df[f"{choice}_comp"]

        # Merge window information with source-receiver distances, not BAz
        df = df.merge(self.srcrcv.drop("backazimuth", axis=1),
                      on=["event", "network", "station"]
                      )

        # Drop unnecessary information except that needed to plot
        # IMPORTANT: Sort by distance so that when the dataframe is iterated on
        #            it starts from the smallest distance and goes up
        df = df.loc[:, ["event", "station", "component", "relative_starttime",
                        "relative_endtime", "distance_km", choice]
                    ].sort_values(by="distance_km")
        if df.empty:
            logger.warning("Filtered dataframe is empty, no windows to plot")
            return

        # Plotting begins here
        f, ax = plt.subplots(figsize=(8, 6))

        # Create a custom color scale based on the min and max values of choice
        if cbar_label is None:
            try:
                # For cleaner formatting of colorbar label
                cbar_label = common_labels[choice]
            except KeyError:
                cbar_label = choice
        if iteration_comp:
            cbar_label = f"DIFF {cbar_label}"

        # Set the bounds of the colorbar
        if isinstance(bounds, list):
            vmin, vmax = bounds
        else:
            if bounds:
                vmax = max(abs(df[choice].min()), abs(df[choice].max()))
                vmin = -1 * vmax
            else:
                vmin = df[choice].min()
                vmax = df[choice].max()
        sm, norm, _ = colormap_colorbar(cmap, vmin=vmin, vmax=vmax,
                                        cbar_label=cbar_label)

        # Determine the global xmin and xmax which will be used more than once
        xmin = df.relative_starttime.min()
        if annotate:
            # Shift to accomodate annotations
            xmin -= anno_shift
        xmax = df.relative_endtime.max()

        dist_values, y_value = [], 0  # keep track of what y-values are used
        for window in df.to_numpy():
            ev, sta, comp, start, end, dist, value = window

            if not distances:
                # Ignore distances and simply plot linearly
                dist_ = y_value
                y_value += rectangle_height
            else:
                # Try not to overlap windows that are very close in distance
                dist_ = int(dist)
                if no_overlap:
                    if dist_ in dist_values:
                        shift, sign = 1, -1
                        while dist_ in dist_values:
                            dist_ += shift
                            # Alternate shift so that we search
                            # 1, -1, 2, -2, 3, -3, etc...
                            shift = sign * (abs(shift) + rectangle_height)
                            sign *= -1
                        dist_values.append(dist_)
                        logger.warning(f"Shifted {ev} {sta}: {dist - dist_}km")
                    else:
                        dist_values.append(dist_)

            # Plot the windows as rectangles to sort of match waveform plots
            ax.add_patch(Rectangle(xy=(start, dist_ - rectangle_height / 2),
                                   width=end - start, ec="k", alpha=alpha,
                                   height=rectangle_height, 
                                   fc=sm.cmap(norm(value)),
                                   zorder=12)
                         )
            # Black background line for frame of reference / gridding
            ax.hlines(y=dist_, xmin=xmin, xmax=xmax, colors="k",
                      alpha=0.3, linewidth=0.3, zorder=10
                      )
            # Annotate event, station, component, distance and value for
            # easier identification. Can be messy with a lot of windows
            if annotate:
                plt.text(xmin, dist_,
                         f"{ev} {sta} {comp} {dist:.2f}km {value:.2f}",
                         fontsize=4.5, zorder=11)

        # Finalize the look of the plot
        plt.title(f"Window Plot: N = {len(df)} "
                  f"[{iteration}{step}] [{iteration_comp}{step_comp}]\n"
                  f"Event: {event} / Station: {station} / Network: {network} / "
                  f"Component: {component}")
        plt.xlabel("Time [s]")
        plt.xlim([xmin, xmax])

        if distances:
            plt.ylabel("Distance [km]")
            plt.ylim([df.distance_km.min() - 10, df.distance_km.max() + 10])
        else:
            # Relative distances means the y-axis values are useless
            plt.ylabel("Relative Distance")
            plt.ylim([-rectangle_height, dist_ + rectangle_height])
            ax.yaxis.set_ticks([])

        if save:
            plt.savefig(save)
        if show:
            plt.show()
        else:
            plt.close()

    def convergence(self, windows="length_s", trials=False, show=True,
                    save=None, normalize=False, float_precision=3,
                    annotate=False, restarts="default", restart_annos=None, 
                    xvalues="model", **kwargs):
        """
        TO DO:
        Separate the sorting functionality from the plotting functionality,
        this function is too confusing.

        Plot the convergence rate over the course of an inversion.
        Scatter plot of total misfit against iteration number, or by step count

        .. note:: 
            Because misfits are floats, they wont be exactly equal, so we need 
            to set some small tolerance in which they can differ

        :type windows: str or bool
        :param windows: parameter to use for Inspector.measurements() to
            determine how to illustrate measurement number, either by:

            * length_s: cumulative window length in seconds
            * nwin: number of misfit windows
            * None: will not plot window information
        :type trials: str
        :param trials: plot the discarded trial step function evaluations from
            the line searches. Useful for understanding optimization efficiency

            * marker: plot trial steps as red x's at their respective misfit val
            * text: annotate the number of trial steps but not their misfit val
        :type normalize: bool
        :param normalize: normalize the objective function values between [0, 1]
        :type float_precision: int
        :param float_precision: acceptable floating point precision for 
            comparisons of misfits. Defaults to 3 values after decimal
        :type restarts: list of int
        :param restarts: If the inversion was restarted, e.g. for parameter
            changes, then the convergence figure should separate two line plots.
            This list allows the User to tell the function where to separate
            the convergence plot. The integers should correspond to indices of
            the Inspector.models attribute.
        :type annotate: bool
        :param annotate: annotate misfit values next to markers
        :type restart_annos: list of str
        :param restart_annos: if restarts is not None, allow annotating text 
            next to each restart. Useful for annotating e.g. parameter changes
            that accompany each restart
        :type xvalues: str
        :param xvalues: How the x-axis should be labelled, available:

            * model: plot the model number under each point
            * eval: number sequentially from 1
        :type show: bool
        :param show: show the plot after making it
        :type save: str
        :param save: file id to save the figure to
        """
        f = kwargs.get("f", None)
        ax = kwargs.get("ax", None)
        dpi = kwargs.get("dpi", 100)
        fontsize = kwargs.get("fontsize", 15)
        anno_fontsize = kwargs.get("anno_fontsize", 15)
        figsize = kwargs.get("figsize", (8, 6))
        legend = kwargs.get("legend", True)
        title = kwargs.get("title", None)
        misfit_label = kwargs.get("misfit_label", "misfit")
        trial_label = kwargs.get("trial_label", "trials")
        window_label = kwargs.get("window_label", "windows")
        trial_color = kwargs.get("trial_color", "r")
        window_color = kwargs.get("window_color", "orange")
        legend_loc = kwargs.get("legend_loc", "best")
        axis_linewidth = kwargs.get("axis_linewidth", 2.)

        # Set some default parameters based on user choices, check parameters
        if windows:
            assert (windows in ["nwin", "length_s"]), \
                "plot_windows must be: 'nwin; or 'length_s'"
        # Default restart values are chosen automatically by the Inspector
        if restarts == "default":
            restarts = self.restarts.index.values

        if restarts is not None and restart_annos is not None:
            assert(len(restarts) + 1 == len(restart_annos)), \
                    "Length of restart anno must match length of `restarts` + 1"

        assert(xvalues in ["model", "eval"]), \
                "xvalues must be 'model' or 'eval'"

        # It may take a while to calculate models so do it once here
        models = self.models
        misfit = models.misfit.round(decimals=float_precision)
        nwin = self.nwin()

        # Set up the figure
        if f is None:
            f, ax = plt.subplots(figsize=figsize, dpi=dpi)
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
                    xlab = "m00"
                    pass
                # If initial eval matches line search final, treat equally
                elif misfit[i] == misfit[j]:
                    continue
                # If they differ, treat them as different points
                else:
                    x += 1
                    xlab = ""
                # xlab = f"{models.model[j]}_r"
            # Status 1 means final evaluation in line search
            elif models.state[j] == 1:
                x += 1
                xlab = f"{models.model[j]}"
            # Status -1 means discarded trial step, plot on the same X value
            elif models.state[j] == -1:
                xdiscards.append(x + 1)  # discards are related to next model
                ydiscards.append(misfit[j])
                continue

            xvals.append(x)
            yvals.append(misfit[j])
            xlabs.append(xlab)

            # Convert restart values from Inspector.models indices
            if restarts is not None and j in restarts:
                xrestarts.append(x)

            # Get the corresponding window number based on iter/step count
            if windows:
                i_ = models.iteration[j]
                s_ = models.step_count[j]
                ywindows.append(nwin.loc[i_].loc[s_][windows])

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
                label = f"{misfit_label}"

            if normalize:
                y_ = [_ / max(y_) for _ in y_]

            line = ax.plot(x_, y_, "o-", linewidth=3,  markersize=10, 
                           c=c, label=label, zorder=10, markeredgecolor="k",
                           markeredgewidth=1.5)
            if annotate:
                for x_anno, y_anno in zip(x_, y_):
                    ax.text(x_anno, y_anno, f"{y_anno:.3f}", zorder=11,
                             fontsize=anno_fontsize)

            if restart_annos:
                ax.text(x_[0], y_[0], restart_annos[idx - 1], zorder=12,
                         fontsize=anno_fontsize, verticalalignment="bottom")

            return line

        # Primary: Two methods of plotting:
        if xrestarts:
            # 1) with user-defined restarts separating legs of the inversion
            first = 0  # first iteration in the current leg
            for i, last in enumerate(xrestarts):
                j = i + 1  # Leg counting should start at 1
                plot_vals(xvals[first:last], yvals[first:last], j)
                first = last
            # Plot the final leg
            lines += plot_vals(xvals[last:], yvals[last:], j + 1)
            if restart_annos:
                ax.text(xvals[last], yvals[last], restart_annos[j])
        else:
            # 2) plot the entire convergence in one line
            lines += plot_vals(xvals, yvals, idx=None)

        # Secondary: Plot number of windows/ window length in a separate axis
        if windows:
            ax2 = ax.twinx()
            # Set ax2 below ax1
            ax.set_zorder(ax2.get_zorder() + 1)
            ax.patch.set_visible(False)
            lines += ax2.plot(xvals, ywindows, "d--", linewidth=2, markersize=8,
                              c=window_color, label=window_label, zorder=5,
                              markeredgecolor="k", markeredgewidth=2
                              )
            ydict = {"length_s": "Cumulative Window Length [s]",
                     "nwin": "Number of Measurements"}
            ax2.set_ylabel(f"{ydict[windows]} (dashed)", rotation=270,
                           labelpad=15., fontsize=fontsize
                           )
            ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
            ax2.yaxis.get_offset_text().set_fontsize(fontsize)

            ax2.tick_params(labelsize=fontsize)

        # Secondary: Plot the discarded trial steps
        if trials == "marker":
            # Scatterplot as red X's to show the misfit value. Not the best
            # because it throws off the scaling of the normal misfit values
            sc = ax.scatter(xdiscards, ydiscards, c=trial_color,
                            marker="x", s=10, zorder=9, label=trial_label
                            )
            lines.append(sc)
        elif trials == "text":
            # Annotate the number of trial steps next to the corresponding value
            for xdiscard in set(xdiscards):
                # Since yvalues are normalized elsewhere, just plot the text
                # near the bottom of the visible axis
                ymin, ymax = ax.get_ylim()
                yval = 0.25 * (ymax - ymin) + ymin

                num_discards = xdiscards.count(xdiscard)
                ax.text(xdiscard, yval, f"{num_discards} trial(s)")

        # Format the axes
        if xvalues.lower() == "model":
            xlabel_ = "Model Number"
            ax.set_xticklabels(xlabs, rotation=60, ha="center")
        elif xvalues.lower() == "eval":
            xlabel_ = "Function Evaluation"
            ax.set_xticklabels(np.arange(1, len(xvals) + 1, 1))
        else:
            xlabel = "Iteration"
    
        ax.set_xlabel(xlabel_, fontsize=fontsize)
        ax.xaxis.set_ticks(xvals)
        ax.set_ylabel("Total Normalized Misfit", fontsize=fontsize)
        ax.tick_params(axis="both", which="major", labelsize=fontsize)

        # Only set ticks on the x-axis
        ax.xaxis.grid(True, which="minor", linestyle=":")
        ax.xaxis.grid(True, which="major", linestyle="-")

        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(axis_linewidth)

        if title is None:
            ax.set_title(f"{self.tag.title()} Convergence\n"
                         f"{len(self.events)} Events / "
                         f"{len(self.stations)} Stations")
        else:
            ax.set_title(title)
                        
        if legend:
            labels = [line.get_label() for line in lines]
            ax.legend(lines, labels, prop={"size": 12}, loc=legend_loc)


        f.tight_layout()

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, ax


def default_axes(ax, cbar=None, **kwargs):
    """
    Ensure that all plots have the same default look. Should be more flexible
    than setting rcParams or having a style sheet. Also allows the same kwargs 
    to be thrown by all functions so that the function calls have the same 
    format.

    Keyword Arguments
    ::
    """
    tick_fontsize = kwargs.get("tick_fontsize", 10)
    tick_linewidth = kwargs.get("tick_linewidth", 1.5)
    tick_length = kwargs.get("tick_length", 5)
    tick_direction = kwargs.get("tick_direction", "in")
    label_fontsize = kwargs.get("label_fontsize", 12)
    axis_linewidth = kwargs.get("axis_linewidth", 2.)
    title_fontsize = kwargs.get("title_fontsize", 14)
    cbar_tick_fontsize = kwargs.get("cbar_tick_fontsize", 10)
    cbar_label_fontsize = kwargs.get("cbar_label_fontsize", 12)
    cbar_outline_color = kwargs.get("cbar_outline_color", "k")
    cbar_linewidth = kwargs.get("cbar_linewdith", 2.)

    # Re-set font sizes for labels already created
    ax.title.set_fontsize(title_fontsize)
    ax.xaxis.label.set_fontsize(label_fontsize)
    ax.yaxis.label.set_fontsize(label_fontsize)
    ax.tick_params(axis="both", which="both", width=tick_linewidth, 
                   direction=tick_direction, labelsize=tick_fontsize, 
                   length=tick_length)

    # Thicken up the bounding axis lines
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(axis_linewidth)

    # Adjust font and bounding bar of colorbar if available
    if cbar is not None:
        cbar.ax.tick_params(labelsize=cbar_tick_fontsize)
        cbar.ax.yaxis.label.set_fontsize(cbar_label_fontsize)
        cbar.outline.set_edgecolor(cbar_outline_color)
        cbar.outline.set_linewidth(cbar_linewidth)


def colormap_colorbar(cmap, vmin=0., vmax=1., dv=None, cbar_label="", 
                      extend="neither"):
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
    cbar = plt.colorbar(sm, boundaries=boundaries, shrink=0.9, pad=0.025,
                        extend=extend)
    if cbar_label:
        cbar.ax.set_ylabel(cbar_label, rotation=270, labelpad=15)

    return sm, norm, cbar


def hover_on_plot(f, ax, obj, values, dissapear=True):
    """
    Allow for hover on a plot for custom annotated information

    .. note::
        This functionality is copied from StackOverflow:
        https://stackoverflow.com/questions/7908636/possible-to-make-labels-\
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
