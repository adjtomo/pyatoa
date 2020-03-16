"""
Functions used to create standard statistical plots for the Inspector
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa.utils.calculate import normalize_a_to_b


class Artist:
    """
    A class of methods for plotting statistics from an Inspector
    """
    def map(self, sta_codes=None, event_ids=None, show=True, save=False,
            **kwargs):
        """
        Plot source and receiver locations

        :type sta_codes: list
        :param sta_codes: unique stations to plot
        :type event_ids: list
        :param event_ids: unique events to plot
        :type show: bool
        :param show: Show the plot
        :type save: str
        :param save: fid to save the given figure
        """
        f, ax = plt.subplots()

        # Plot all sources as scatter
        event_x, event_y, event_s = [], [], []
        for event in self.event_ids:
            if event_ids and event not in event_ids:
                continue
            event_s.append(f"{event}\n"
                           f"M{self.srcrcv[event]['mag']:.2f}\n"
                           f"{self.srcrcv[event]['time'][:10]}\n"
                           f"{self.srcrcv[event]['depth_m']*1E-3:.2f}km")
            event_x.append(self.srcrcv[event]["utm_x"])
            event_y.append(self.srcrcv[event]["utm_y"])
        sc_events = plt.scatter(event_x, event_y, marker="o", c="w",
                                edgecolors="r", s=10, zorder=100)

        # Plot all receivers as scatter
        station_x, station_y, station_s = [], [], []
        for sta in self.stations:
            if sta_codes and sta not in sta_codes:
                continue
            station_s.append(sta)
            station_x.append(self.srcrcv[sta]["utm_x"])
            station_y.append(self.srcrcv[sta]["utm_y"])
        sc_stations = plt.scatter(station_x, station_y, marker="v",
                                  edgecolors="g", c="w", s=10, zorder=100)

        plt.grid(which="both", linestyle=":")
        plt.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
        plt.xlabel("Easting (m)")
        plt.ylabel("Northing (m)")
        plt.title(f"Source-Receiver Map;\n"
                  f"N_sta={len(station_x)}; N_ev={len(event_x)}")

        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, sc_events, event_s, **kwargs)
            hover_on_plot(f, ax, sc_stations, station_s, **kwargs)
            plt.show()

        return f, ax

    def windows_by_distance(self, model, choice="cc_shift_sec",
                            show=True, save=False, event_id=None,
                            sta_code=None, cha_code=None, color_by=None):
        """
        Make a plot of window attributes versus source-receiver distance

        :type model: str
        :param model: model to query, e.g. "m00"
        :type choice: str
        :param choice: variable to query for windows
            available: "cc_shift_sec", "dlna", "length_s", "max_cc", "weight"
        :type event_id: str
        :param event_id: only plot for a given event id
        :type sta_code: str
        :param sta_code: only plot for a given station
        :type cha_code: str
        :param cha_code: only plot for a given channel
        :type color_by: str
        :param color_by: color the markers based on similar groupings, None
            defaults to no color. Available: 'event', 'channel'
        :type show: bool
        :param show: Show the plot
        :type save: str
        :param save: fid to save the given figure
        :type event_id: str
        :param event_id: only plot for a given event
        """
        choices = ["cc_shift_sec", "dlna", "length_s", "max_cc", "weight"]
        assert (choice in choices), f"Choice must be in {choices}"
        assert self.srcrcv, "No distance information"
        assert (model in self.models), (f"Model must be in "
                                             f"{self.models}")
        if event_id:
            assert (event_id in self.event_ids), \
                f"event_id must be in {self.event_ids}"
        if sta_code:
            assert (sta_code in self.stations), \
                f"sta_code must be in {self.stations}"
    
        # Collect misfit and distance information per event
        f, ax = plt.subplots()
        stations, distances, values, cidx = [], [], [], []
        for i, event in enumerate(self.windows):
            if event_id and event != event_id:
                continue
            for j, sta in enumerate(self.windows[event][model]):
                if sta_code and sta != sta_code:
                    continue
                for k, cha in enumerate(self.windows[event][model][sta]):
                    if cha_code and cha != cha_code:
                        continue
                    window = self.windows[event][model][sta][cha][choice]
                    # For each window, assign a point
                    for value in window:
                        stations.append(f"{sta}.{cha}\n{event}\n{value:.2f}")
                        distances.append(
                            self.srcrcv[event][sta]["dist_km"])
                        values.append(value)
                        if color_by == "event":
                            cidx.append(i)
                        # Wouldn't work because station list isn't the same
                        # length each time
                        # elif color_by == "station":
                        #     cidx.append(j)
                        elif color_by == "channel":
                            cidx.append(k)
                        else:
                            cidx.append("w")
    
        # Quick plot so each event gets a different color
        sc = plt.scatter(distances, values, c=cidx, s=7, label=event,
                         marker="o", edgecolors="k", cmap=plt.cm.jet)
    
        plt.title(f"{choice} vs. event-station distance; N={len(values)}\n"
                  f"model: {model} / event: {event_id} / station: {sta_code}")
        plt.xlabel("Distances (km)")
        plt.ylabel(choice)
        # plt.ylim([0, np.ceil(max(values))])
        plt.grid(which="both", linestyle="--", alpha=0.5, linewidth=.5)
    
        hover = hover_on_plot(f, ax, sc, stations)
    
        if save:
            plt.savefig(save)
        if show:
            f.canvas.mpl_connect("motion_notify_event", hover)
            plt.show()
    
        return f, ax

    def misfit_by_distance(self, model, show=True, save=False,
                           event_id=None, sta_code=None, color_by=None):
        """
        Make a plot of misfit versus source-receiver distance

        :type model: str
        :param model: model to query, e.g. "m00"
        :type event_id: str
        :param event_id: only plot for a given event id
        :type sta_code: str
        :param sta_code: only plot for a given station
        :type show: bool
        :param show: Show the plot
        :type save: str
        :param save: fid to save the given figure
        :type event_id: str
        :param event_id: only plot for a given event
        :type color_by: str
        :param color_by: color the markers based on similar groupings, None
            defaults to no color. Available: 'event'
        """
        assert self.srcrcv, "No distance information"
        assert(model in self.models), f"Model must be in {self.models}"
        if event_id:
            assert(event_id in self.event_ids), \
                f"event_id must be in {self.event_ids}"
        if sta_code:
            assert(sta_code in self.stations), \
                f"sta_code must be in {self.stations}"
    
        # Collect misfit and distance information per event
        f, ax = plt.subplots()
        stations, distances, misfits, cidx = [], [], [], []
        for i, event in enumerate(self.misfits):
            if event_id and event != event_id:
                continue
            for j, sta in enumerate(self.misfits[event][model]):
                if sta_code and sta != sta_code:
                    continue
                misfit = self.misfits[event][model][sta]["msft"]
                stations.append(f"{sta}\n{event}\n{misfit:.2f}")
                distances.append(self.srcrcv[event][sta]["dist_km"])
                misfits.append(misfit)
                if color_by == "event":
                    cidx.append(i)
                else:
                    cidx.append("w")
    
        # Quick plot so each event gets a different color
        sc = plt.scatter(distances, misfits, c=cidx, s=7, label=event,
                         marker="o", edgecolors="k")
    
        plt.title(f"misfit vs. event-station distance; N={len(values)}\n"
                  f"model: {model} / event: {event_id} / station: {sta_code}")
        plt.xlabel("Distances (km)")
        plt.ylabel("Misfit")
        plt.ylim([0, np.ceil(max(misfits))])
        plt.grid(which="both", linestyle="--", alpha=0.5, linewidth=.5)
    
        hover = hover_on_plot(f, ax, sc, stations)
    
        if save:
            plt.savefig(save)
        if show:
            f.canvas.mpl_connect("motion_notify_event", hover)
            plt.show()
    
        return f, ax

    def plot_paths(self, model, choice, event_id=None, sta_code=None,
                   threshold=None, hover_on_lines=False, annotate=True,
                   colormap=plt.cm.Spectral_r, color="k",
                   normalize_to_model=False,
                   show=True, save=None, **kwargs):
        """
        Plot misfit by source-receiver path to try to highlight portions of
        the model that may give rise to larger misfit

        :type model: str
        :param model: model to choose for misfit
        :type choice: str
        :param choice: choice between `misfit` and `cc_shift_sec` time shift
        :type event_id: str
        :param event_id: only plot for a given event id
        :type sta_code: str
        :param sta_code: only plot for a given station
        :type threshold: float
        :param threshold: normalized misfit value below which, paths will not be
            plotted. Good for looking at only high misfit values. Values must
            be between 0 and 1
        :type hover_on_lines: bool
        :param hover_on_lines: for interactive plots, show misfit values when
            hovering over the source-receiver raypath lines. This can get a bit
            messy with a lot of lines
        :type colormap: matplotlib.colors.ListedColormap
        :param colormap: colormap for coloring lines
        :type normalize_to_model: bool
        :param normalize_to_model: normalize the misfit value to the entire
            model. Only relevant is event_id or sta_code is specified. Defaults
            to False, colors will be normalized to the misfits shown
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        choices = ["misfit", "cc_shift_sec"]
        assert(choice in choices), f"choice must be in {choices}"

        # Pick which value is going to be plotted
        if choice == "misfit":
            misfits = self.sort_misfits_by_model()
            values = self.misfit_values(model)
        elif choice == "cc_shift_sec":
            misfits = self.sort_windows_by_model()
            values = self.window_values(model, choice)
        assert(model in misfits), f"model must be in {misfits.keys()}"

        # Ensure that optional arguments are acceptable
        if event_id:
            assert(event_id in self.event_ids), \
                f"event_id must be in {self.event_ids}"
        if sta_code:
            assert(sta_code in self.stations), \
                f"sta_code must be in {self.stations}"
        if threshold:
            assert(0 <= threshold <= 1), "Threshold must be between 0 and 1"
    
        # Empty values to be filled when looping through data
        x_values, y_values, alphas, misfit_list = [], [], [], []
        coverage = 0  # to count number of srcs or rcvs covered

        for event in misfits[model]:
            # Skip event if requested
            if event_id and event != event_id:
                continue
            # Loop through full station list rather than subset that have misfits
            for sta in self.stations:
                # Skip station if requested
                if sta_code and sta != sta_code:
                    continue
                # Here we determine the value that is to be used as the color
                try:
                    # For misfit, just append misfit value, normalize if wanted
                    if choice == "misfit":
                        misfit = misfits[model][event][sta]["msft"]
                        if normalize_to_model:
                            misfit /= max(values)

                    # Take the mean of all time shifts for all measurements
                    elif choice == "cc_shift_sec":
                        misfit = n_windows = 0
                        for i, cha in enumerate(misfits[model][event][sta]):
                            cs = misfits[model][event][sta][cha]["cc_shift_sec"]
                            n_windows += len(cs)
                            misfit += sum([abs(_) for _ in cs])
                        misfit /= n_windows

                    # Ignore any misfit values below a certain threshold
                    if threshold and misfit < threshold:
                        continue
                except KeyError:
                    continue
    
                # Append the relevant values to plot srcrcv lines later
                x_values.append([self.srcrcv[event]["utm_x"],
                                 self.srcrcv[sta]["utm_x"]]
                                )
                y_values.append([self.srcrcv[event]["utm_y"],
                                 self.srcrcv[sta]["utm_y"]]
                                )
                misfit_list.append(misfit)
                coverage += 1

        # Instantiate plot
        f, ax = plt.subplots(figsize=(8, 7))

        # Set the colorbar based on the misfit list retrieved
        if colormap:
            try:
                vmin, vmax = min(misfit_list), max(misfit_list)
            except ValueError:
                vmin = vmax = 0
            if normalize_to_model:
                vmin = min([abs(_) for _ in values])
                vmax = max([abs(_) for _ in values])
            cbar_label = choice
            if normalize_to_model:
                cbar_label = f"normalized {choice}"
            cmap, norm, cbar = colormap_colorbar(colormap, vmin=vmin, vmax=vmax,
                                                 cbar_label=cbar_label, 
                                                 **kwargs)

        # Plot the misfit lines
        for x, y, m in zip(x_values, y_values, misfit_list):
            # Change the alpha for aggregate plots to remove visual clutter
            if colormap:
                alpha = misfit
                if sta_code or event_id:
                    alpha = None
                    c = cmap(norm(m))
            # If no colormap is given, just plot source receiver coverage
            else:
                c = color
                alpha = 0.05
            line, = plt.plot(x, y, c=c, alpha=alpha, zorder=10 + m)
            if hover_on_lines:
                hover_on_plot(f, ax, line, [f"{m:.2f}"], dissapear=True)

        # Plot all sources as scatter
        event_x, event_y, event_s = [], [], []
        for event in misfits[model]:
            event_s.append(event)
            event_x.append(self.srcrcv[event]["utm_x"])
            event_y.append(self.srcrcv[event]["utm_y"])
        sc_events = plt.scatter(event_x, event_y, marker="o", c="w",
                                edgecolors="r", s=10, zorder=100)

        # Plot all receivers as scatter
        station_x, station_y, station_s = [], [], []
        for sta in self.stations:
            station_s.append(sta)
            station_x.append(self.srcrcv[sta]["utm_x"])
            station_y.append(self.srcrcv[sta]["utm_y"])
        sc_stations = plt.scatter(station_x, station_y, marker="v",
                                  edgecolors="g", c="w", s=10, zorder=100)

        # Annotation information
        append_anno = ""
        if sta_code and not event_id:
            append_anno += (
                f"{coverage}/{len(self.event_ids)} events\n"
                f"{1E2*coverage/len(self.event_ids):.2f}% coverage")
        elif event_id and not sta_code:
            append_anno += (
                f"{coverage}/{len(self.stations)} stations\n"
                f"{1E2*coverage/len(self.stations):.2f}% coverage\n")

        # Determine min, max and mean misfit for annotation
        min_misfit = max_misfit = mean_misfit = "NaN"
        if misfit_list:
            min_misfit = f"{np.min(misfit_list):.2f}"
            max_misfit = f"{np.max(misfit_list):.2f}"
            mean_misfit = f"{np.mean(misfit_list):.2f}"

        # Fill in the rest of the plot attributes
        plt.grid(which="both", linestyle=":")
        plt.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
        plt.xlabel("Easting (m)")
        plt.ylabel("Northing (m)")
        plt.title(f"Source-Receiver misfit; N={len(misfit_list)}")
        
        anno = (f"model: {model}\n"
                f"event: {event_id}\n"
                f"station: {sta_code}\n"
                f"min: {min_misfit}\n"
                f"max: {max_misfit}\n"
                f"mean: {mean_misfit}\n")
        if append_anno:
            anno += append_anno

        if annotate:
            annotate_txt(ax, anno, **kwargs)
    
        # Make source and receiver markers interactive
        if save:
            plt.savefig(save)
        if show:
            hover_on_plot(f, ax, sc_events, event_s, **kwargs)
            hover_on_plot(f, ax, sc_stations, station_s, **kwargs)
            plt.show()
    
        return f, ax
    
    def event_depths(self, xaxis="x", show=True, save=None):
        """
        Create a scatter plot of events at depth

        :type xaxis: str
        :param xaxis: variable to use as the x-axis on the plot, can either be
            'x', 'y', 'lat' or 'lon'
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        choices = ["x", "y", "lat", "lon"]
        assert(xaxis in choices), f"xaxis must be in {choices}"
    
        # Plot initializations
        f, ax = plt.subplots(figsize=(8, 6))
        label_dict = {"x": "Easting (m)", "y": "Northing (m)",
                      "lat": "Latitude (deg)", "lon": "Longitude (deg)"
                      }
    
        depth_vals = self.depths
        mag_vals = self.mags
        x_vals, depths, mags, events, s = [], [], [], [], []
        for event in self.event_ids:
            # Differentiate between lat/lon and utm, allow User to plot for both
            # x or y on the x axis
            if xaxis in ["x", "y"]:
                x = self.srcrcv[event]["utm_x"]
                y = self.srcrcv[event]["utm_y"]
            else:
                x = self.srcrcv[event]["lon"]
                y = self.srcrcv[event]["lat"]
            if xaxis in ["x", "lon"]:
                x_vals.append(x)
            elif xaxis in ["y", "lat"]:
                x_vals.append(y)
    
            # Get some other information to use in the plot
            depths.append(depth_vals[event])
            mags.append(mag_vals[event])
            s.append(f"{event}\nM{mag_vals[event]:.2f}; "
                     f"{depth_vals[event]*1E-3:.2f}km")
    
        # Adjust units before plotting
        depths = [_ * 1E-3 for _ in depths]
        mags = normalize_a_to_b(mags, 100, 500)
    
        # Inverted axis for positive depth values
        if depths[0] > 0:
            plt.gca().invert_yaxis()
    
        # Scatter plot
        sc = plt.scatter(x_vals, depths, s=mags, c="None", marker="o",
                         edgecolors="k")
        plt.xlabel(label_dict[xaxis])
        plt.ylabel("Depth (m)")
        plt.title(f"Event Depth Cross Section; N={len(depths)}")
        plt.grid(which="both", linestyle=":", alpha=0.5)
        hover_on_plot(f, ax, sc, s, dissapear=True)
    
        if save:
            plt.savefig(save)
        if show:
            plt.show
    
    def misfit_histogram(self, model, model_comp=None, choice="cc_shift_sec",
                         binsize=1., show=True, save=None, title=None, 
                         anno=True):
        """
        Create a histogram of misfit information for either time shift or
        amplitude differences

        :type model: str
        :param model: model to choose for misfit
        :type model_comp: str
        :param model_comp: model to compare with, will be plotted in front
        :type choice: str
        :param choice: choice of 'cc_shift_s' for time shift, or 'dlna' as
            amplitude
        :type binsize: float
        :param binsize: size of the histogram bins
        :type show: bool
        :param show: show the plot
        :type save: str
        :param save: fid to save the figure
        """
        def get_stats(n, bins):
            """get stats from a histogram"""
            mids = 0.5 * (bins[1:] + bins[:-1])
            mean = np.average(mids, weights=n)
            var = np.average((mids - mean)**2, weights=n)
            std = np.sqrt(var)

            return mean, var, std

        def retrieve_val(model):
            """retrieve some model information"""
            if choice == "misfit":
                val = self.misfit_values(model)
            else:
                val = self.window_values(model, choice)

            # Determine bound for histogram
            min_value = np.floor(min(val))
            max_value = np.ceil(max(val))
            bound = max(abs(min_value), abs(max_value))

            return val, bound

        choices = ["misfit", "dlna", "cc_shift_sec", "length_s"]
        assert(choice in choices), f"choice must be in {choices}"
        assert(model in self.models), f"model {model} not in Inspector"
        if model_comp:
            pass
            # assert(model != model_comp), "model_comp must be different"
    
        label_dict = {"cc_shift_sec": "Time Shift (s)",
                      "dlna": "dlnA=ln(A_obs/A_syn)",
                      "misfit": "misfit",
                      "length_s": "Window Length (s)"}
   
        # Get the values to plot on histogram, and the bounds of the bins
        val, bound = retrieve_val(model)
        if model_comp:
            val_comp, bound_comp = retrieve_val(model_comp)
            
        # Make the main histogram
        f, ax = plt.subplots(figsize=(8, 6))
        n, bins, patches = plt.hist(
                 x=val, bins=len(np.arange(-1*bound, bound, binsize)),
                 color="darkorange",  histtype="bar", edgecolor="black",
                 linewidth=2.5, label=f"{model}; N={len(val)}", zorder=10
                 )
        mu1, var1, std1 = get_stats(n, bins)

        # If a model is to be compared to the initial model, plot infront,
        # transparent, with overlying lines, and behind as full
        if model_comp:
            n, bins, pathces = plt.hist(
                     x=val_comp, 
                     label=f"{model_comp}; N={len(val_comp)}",
                     bins=len(np.arange(-1 * bound_comp, bound_comp, binsize)),
                     histtype="bar", edgecolor="black", linewidth=2.,
                     alpha=.8, zorder=9)
            plt.hist(x=val_comp,
                     bins=len(np.arange(-1 * bound_comp, bound_comp, binsize)),
                     histtype="step", edgecolor="deepskyblue", linewidth=2.,
                     alpha=0.7, zorder=11)
            mu2, var2, std2 = get_stats(n, bins)
    
        # dln(A) information only relevant in these bounds
        if choice == "dlna":
            plt.xlim([-1.75, 1.75])

        # Annotate stats information
        if anno:
            anno = f"$\mu({model})\pm\sigma({model})$={mu1:.2f}$\pm${std1:.2f}"
            if model_comp:
                anno += (f"\n$\mu({model_comp})\pm\sigma({model_comp})$="
                         f"{mu2:.2f}$\pm${std2:.2f}")
            # annotate_txt(ax, anno, "upper-left", fontsize=12, zorder=12) 
    
        # Finalize plot details
        plt.xlabel(label_dict[choice])
        plt.ylabel("Count")
        if not title:
            title = f"{choice} histogram"
        plt.title(anno)
        plt.tick_params(which='both', direction='in', top=True, right=True)
        plt.grid(linewidth=1, linestyle=":", which="both", zorder=1)
        plt.axvline(x=0, ymin=0, ymax=1, linewidth=1.5, c="k", zorder=2, 
                    alpha=0.75, linestyle='--')
        plt.legend()
    
        if save:
            plt.savefig(save)
        if show:
            plt.show()
    
        return f, ax

    def plot_windows(self, model, event_id=None, sta_choice=None,
                     component=None, color_by_comp=False,
                     velocities=[2, 4, 7], alpha=0.2):
        """
        Show lengths of windows chosen based on source-receiver distance, akin
        to Tape's Thesis or to the LASIF plots. These are useful for showing
        which phases are chosen, and window choosing behavior as distance
        increases and (probably) misfit increases.

        :type model: str
        :param model: model to analyze
        :type component: str
        :param component: choose a specific component to analyze
        :type color_by_comp: bool
        :param color_by_comp: individually color each component
        :type velocities: list of floats
        :param velocities: for lines showing given wavespeeds for comparison
            against certain seismic phases, if none, no lines plotted
        :return:
        """
        windows = self.sort_windows_by_model()[model]
        comp_dict = {'Z': 'r', 'N': 'g', 'R': 'g', 'E': 'b', 'T': 'b'}
        if component:
            assert(component.upper() in comp_dict), \
                f"component not in {comp_dict.keys()}"

        starts, ends, dists, labels, colors = [], [], [], [], []
        for event in windows:
            if event_id and event != event_id:
                continue
            for sta in windows[event]:
                if sta_choice and sta != sta_choice:
                    continue
                distance = self.srcrcv[event][sta]["dist_km"]
                for cha in windows[event][sta]:
                    comp = cha[-1]
                    if component and component != comp:
                        continue
                    for i, (start, end) in enumerate(zip(
                            windows[event][sta][cha]["rel_start"],
                            windows[event][sta][cha]["rel_end"])):
                        starts.append(start)
                        ends.append(end)
                        dists.append(distance)
                        colors.append(comp_dict[comp.upper()])
                        labels.append(f"{event}, {sta}, {cha}, {i}")

        f, ax = plt.subplots(figsize=(8, 6))
        if not color_by_comp:
            colors = "k"
        
        # Plot the lines that represent the windows
        # Set alpha by the length of window, short windows darker
        alphas = 1/(np.array(ends) - np.array(starts))
        alphas = normalize_a_to_b(alphas, 0.2, .5)
        linewidths = normalize_a_to_b(alphas, 0.2, .9)
        for i in range(len(dists)):
            l = ax.hlines(y=dists[i], xmin=starts[i], xmax=ends[i], 
                          colors=colors, alpha=alphas[i], 
                          linewidth=linewidths[i])
            ax.hlines(y=dists[i], xmin=0, xmax=300, colors="k", alpha=0.1,
                      linewidth=0.1)

        # Define lines starting from the start of the data, slopes related to
        # certain wavespeeds (km/s)
        if velocities:
            def wavespeed_line(x, v):
                """Slope of the line of a given wavespeed"""
                return v * x - (v * min(x))

            def trailing_edge(y):
                """Slope of the line of the trailing edge of the plot"""
                m = (max(dists) - min(dists)) / (max(ends) - min(ends))
                return (y - min(dists)) / m + min(ends)

            x_vals = np.linspace(0, 300, 100)
            for vel in velocities:
                # Conditions for stopping the line plotting
                for i, (x_, y_) in enumerate(
                        zip(x_vals, wavespeed_line(x_vals, vel))):
                    if y_ >= max(dists):
                        break
                    elif x_ > trailing_edge(y_):
                        break
                plt.plot(x_vals[:i], wavespeed_line(x_vals[:i], vel),
                         linewidth=1.25, label=f"{vel}km/s", alpha=0.75)

        plt.title(f"{len(dists)} misfit windows")
        plt.xlabel("Time [s]")
        plt.ylabel("Distance [km]")
        if velocities:
            plt.legend(fontsize=8)
        plt.grid(linestyle='--', linewidth=.5, alpha=.5)
        plt.xlim([min(starts)-10, max(ends)+10])
        plt.ylim([min(dists)-10, max(dists)+10])

        plt.show()

    def convergence(self, choice="iter", linewidth=2., markersize=8., c="k",
                    show=True, save=None):
        """
        Plot the convergence rate over the course of an inversion

        :type choice: str
        :param choice: choice between plotting "step" lengths or "iter" ations
        :type linewidth: float
        :param linewidth: line width for the connecting lines
        :type markersize: float
        :param markersize: marker size of the points
        :type show: bool
        :param show: show the plot after making it
        :type save: str
        :param save: file id to save the figure to
        """
        # Set up the values to plot
        misfits = self.sum_misfits()
        windows = self.cum_win_len()
        models = list(misfits.keys())
        xvalues = np.arange(0, 10, 1)
        misfits = list(misfits.values())
        windows = list(windows.values())

        f, ax1 = plt.subplots(figsize=(8, 6))
        ax2 = ax1.twinx()
        ax1.plot(xvalues, misfits, 'o-', c=c, linewidth=linewidth,
                 markersize=markersize)
        ax2.plot(xvalues, windows, 'v--', c=c, linewidth=linewidth,
                 markersize=markersize)

        ax1.set_xlabel("Model Number")
        ax1.set_ylabel("Total Normalized Misfiti (solid)")
        ax2.set_ylabel("Cumulative Window Length [s] (dashed)", rotation=270, 
                       labelpad=15.)
        ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
        ax1.grid(True, alpha=0.5, linestyle='--', linewidth=1.)

        # Change the labels to the model numbers
        labels = [item.get_text() for item in ax1.get_xticklabels()]
        for i, model in enumerate(models):
            labels[i] = model
        ax1.set_xticklabels(labels)

        if save:
            plt.savefig(save)
        if show:
            plt.show()

        plt.close()


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


def hover_on_plot(f, ax, obj, values, dissapear=False, **kwargs):
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

