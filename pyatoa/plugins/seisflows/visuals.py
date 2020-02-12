"""
Functions used to create standard statistical plots for the Inspector
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa.utils.tools.srcrcv import lonlat_utm
from pyatoa.utils.tools.calculate import normalize_a_to_b


def windows_by_distance(insp, model, choice="cc_shift_sec",
                        show=True, save=False, event_id=None,
                        sta_code=None, cha_code=None, color_by=None):
    """
    Make a plot of window attributes versus source-receiver distance

    :type insp: pyatoa.plugins.seisflows.inspector.Inspector
    :param insp: Inspector object containing the relevant information
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
    assert insp.srcrcv, "No distance information"
    assert (model in insp.models), f"Model must be in {insp.models}"
    if event_id:
        assert (event_id in insp.event_ids), \
            f"event_id must be in {insp.event_ids}"
    if sta_code:
        assert (sta_code in insp.stations), \
            f"sta_code must be in {insp.stations}"

    # Collect misfit and distance information per event
    f, ax = plt.subplots()
    stations, distances, values, cidx = [], [], [], []
    for i, event in enumerate(insp.windows):
        if event_id and event != event_id:
            continue
        for j, sta in enumerate(insp.windows[event][model]):
            if sta_code and sta != sta_code:
                continue
            for k, cha in enumerate(insp.windows[event][model][sta]):
                if cha_code and cha != cha_code:
                    continue
                window = insp.windows[event][model][sta][cha][choice]
                # For each window, assign a point
                for value in window:
                    stations.append(f"{sta}.{cha}\n{event}\n{value:.2f}")
                    distances.append(insp.srcrcv[event][sta]["dist_km"])
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


def misfit_by_distance(insp, model, show=True, save=False,
                       event_id=None, sta_code=None, color_by=None):
    """
    Make a plot of misfit versus source-receiver distance

    :type insp: pyatoa.plugins.seisflows.inspector.Inspector
    :param insp: Inspector object containing the relevant information
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
    assert insp.srcrcv, "No distance information"
    assert(model in insp.models), f"Model must be in {insp.models}"
    if event_id:
        assert(event_id in insp.event_ids), \
            f"event_id must be in {insp.event_ids}"
    if sta_code:
        assert(sta_code in insp.stations), \
            f"sta_code must be in {insp.stations}"

    # Collect misfit and distance information per event
    f, ax = plt.subplots()
    stations, distances, misfits, cidx = [], [], [], []
    for i, event in enumerate(insp.misfits):
        if event_id and event != event_id:
            continue
        for j, sta in enumerate(insp.misfits[event][model]):
            if sta_code and sta != sta_code:
                continue
            misfit = insp.misfits[event][model][sta]["msft"]
            stations.append(f"{sta}\n{event}\n{misfit:.2f}")
            distances.append(insp.srcrcv[event][sta]["dist_km"])
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


def misfit_by_path(insp, model, event_id=None, sta_code=None,
                   threshold=None, hover_on_lines=False, show_all=True,
                   colormap=plt.cm.Spectral_r, show=True, save=None,
                   **kwargs):
    """
    Plot misfit by source-receiver path to try to highlight portions of
    the model that may give rise to larger misfit

    :type insp: pyatoa.plugins.seisflows.inspector.Inspector
    :param insp: Inspector object containing the relevant information
    :type model: str
    :param model: model to choose for misfit
    :type event_id: str
    :param event_id: only plot for a given event id
    :type sta_code: str
    :param sta_code: only plot for a given station
    :type threshold: float
    :param threshold: normalized misfit value below which, paths will not be
        plotted. Good for looking at only high misfit values. Values must
        be between 0 and 1
    :type show_all: bool
    :param show_all: Only relevant is threshold is not None. Show the other
        stations or events if plotting specific station or event
    :type hover_on_lines: bool
    :param hover_on_lines: for interactive plots, show misfit values when
        hovering over the source-receiver raypath lines. This can get a bit
        messy with a lot of lines
    :type colormap: matplotlib.colors.ListedColormap
    :param colormap: colormap for coloring lines
    :type show: bool
    :param show: show the plot
    :type save: str
    :param save: fid to save the figure
    """
    misfits = insp.sort_misfits_by_model()
    assert(model in misfits)
    assert(model in insp.models), f"Model must be in {insp.models}"
    if event_id:
        assert(event_id in insp.event_ids), \
            f"event_id must be in {insp.event_ids}"
    if sta_code:
        assert(sta_code in insp.stations), \
            f"sta_code must be in {insp.stations}"
    if threshold:
        assert(0 <= threshold <= 1), "Threshold must be between 0 and 1"

    # Instantiate plot
    f, ax = plt.subplots(figsize=(8, 7))
    cmap, norm, cbar = colormap_colorbar(colormap, **kwargs)

    # Empty lists to be filled when looping through data
    stations_plotted = []
    event_x, event_y, event_s = [], [], []
    station_x, station_y, station_s = [], [], []
    msftlist = []

    # Used to count the number of events that a station recorded
    # or the number of stations that recorded an event
    coverage = 0

    for event in misfits[model]:
        # Skip event if requested
        if event_id and event != event_id:
            continue

        ev_x, ev_y = lonlat_utm(lon_or_x=insp.srcrcv[event]["lon"],
                                lat_or_y=insp.srcrcv[event]["lat"],
                                utm_zone=-60, inverse=False
                                )

        # Append event information to list for plotting
        # Skip adding stations if thresholding and show all is turned off
        event_s.append(event)
        event_x.append(ev_x)
        event_y.append(ev_y)

        # Loop through full station list rather than subset that have misfits
        coords = insp.coords
        for sta in insp.stations:
            # Convert station coordinates and append to list
            sta_x, sta_y = lonlat_utm(
                lon_or_x=coords[sta]["lon"], lat_or_y=coords[sta]["lat"],
                utm_zone=-60, inverse=False
            )
            # Ensure we only plot each station once
            if sta not in stations_plotted:
                stations_plotted.append(sta)
                station_x.append(sta_x)
                station_y.append(sta_y)
                station_s.append(sta)

            # Skip station if requested
            if sta_code and sta != sta_code:
                continue

            # Normalize misfit by the largest value, plot srcrcv as line
            try:
                misfit = (misfits[model][event][sta]["msft"] /
                          max(insp.misfit_values(model))
                          )
                msftlist.append(misfit)
                coverage += 1
            except KeyError:
                continue

            # Ignore misfit values below a certain threshold if requested
            if threshold and misfit < threshold:
                continue

            # Change the alpha for aggregate plots to remove visual clutter
            alpha = misfit
            if sta_code or event_id:
                alpha = None

            # Plot the line between src and rcv colored by misfit
            line, = plt.plot([ev_x, sta_x], [ev_y, sta_y],
                             c=cmap(norm(misfit)), alpha=alpha,
                             zorder=10 + misfit)
            # Set hover on for the misfit lines, can get messy
            if hover_on_lines:
                hover_on_plot(f, ax, line, [f"{misfit:.2f}"],
                              dissapear=True)

    # Determine the coverage if a single station or event was chosen
    if sta_code and not event_id:
        cov = f"{coverage}/{len(insp.event_ids)} events\n" \
              f"{1E2*coverage/len(insp.event_ids):.2f}% coverage"
    elif event_id and not sta_code:
        cov = f"{coverage}/{len(insp.stations)} stations\n" \
              f"{1E2*coverage/len(insp.stations):.2f}% coverage\n"

    # Plot sources and receivers as scatterplots
    sc_events = plt.scatter(event_x, event_y, marker="o", c="w",
                            edgecolors="r", s=10, zorder=100)
    sc_stations = plt.scatter(station_x, station_y, marker="v",
                              edgecolors="g", c="w", s=10, zorder=100)

    # Make the plot a bit prettier
    plt.grid(which="both", linestyle=":")
    plt.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
    plt.xlabel("Easting (m)")
    plt.ylabel("Northing (m)")
    plt.title("Source-Receiver misfit")

    anno = (f"model: {model}\n"
            f"event: {event_id}\n"
            f"station: {sta_code}\n"
            f"min misfit: {np.min(msftlist):.2f}\n"
            f"max misfit: {np.max(msftlist):.2f}\n"
            f"mean misfit: {np.mean(msftlist):.2f}\n")
    if cov:
        anno += cov

    annotate_txt(ax, anno, **kwargs)

    # Make source and receiver markers interactive
    if save:
        plt.savefig(save)
    if show:
        hover_on_plot(f, ax, sc_events, event_s)
        hover_on_plot(f, ax, sc_stations, station_s)
        plt.show()

    return f, ax


def event_depths(insp, xaxis="x", show=True, save=None):
    """
    Create a scatter plot of events at depth

    :type insp: pyatoa.plugins.seisflows.inspector.Inspector
    :param insp: Inspector object containing the relevant information
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

    depth_vals = insp.depths
    mag_vals = insp.mags
    x_vals, depths, mags, events, s = [], [], [], [], []
    for event in insp.event_ids:
        # Differentiate between lat/lon and utm, allow User to plot for both
        # x or y on the x axis
        if xaxis in ["x", "y"]:
            x, y = lonlat_utm(lon_or_x=insp.srcrcv[event]["lon"],
                              lat_or_y=insp.srcrcv[event]["lat"],
                              utm_zone=-60, inverse=False
                              )
        else:
            x = insp.srcrcv[event]["lon"]
            y = insp.srcrcv[event]["lat"]
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


def misfit_histogram(insp, model, model_comp=None, choice="cc_shift_sec",
                     binsize=1., show=True, save=None):
    """
    Create a histogram of misfit information for either time shift or
    amplitude differences

    :type insp: pyatoa.plugins.seisflows.inspector.Inspector
    :param insp: Inspector object containing the relevant information
    :type model: str
    :param model: model to choose for misfit
    :type model_comp: str
    :param model_comp: model to compare with, will be plotted in front
    :type choice: str
    :param choice: choice of 'cc_shift_s' for time shift, or 'dlna' as amplitude
    :type binsize: float
    :param binsize: size of the histogram bins
    :type show: bool
    :param show: show the plot
    :type save: str
    :param save: fid to save the figure
    """
    choices = ["misfit", "dlna", "cc_shift_sec", "length_s"]
    assert(choice in choices), f"choice must be in {choices}"
    assert(model in insp.models), f"model {model} not in Inspector"
    if model_comp:
        assert(model != model_comp), "model_comp must be different"

    label_dict = {"cc_shift_sec": "Time Shift (s)",
                  "dlna": "dlnA=ln(A_obs/A_syn)",
                  "misfit": "misfit",
                  "length_s": "Window Length (s)"}

    # Retrive list of relevant data
    if choice == "misfit":
        val = insp.misfit_values(model)
    else:
        val = insp.window_values(model, choice)

    # Determine bound for histogram
    min_value = np.floor(min(val))
    max_value = np.ceil(max(val))
    bound = max(abs(min_value), abs(max_value))

    # Make the main histogram
    f, ax = plt.subplots(figsize=(8, 6))
    plt.hist(x=val, bins=len(np.arange(-1*bound, bound, binsize)),
             color="darkorange",  histtype="bar", edgecolor="black",
             linewidth=2.5, label=model, zorder=10
             )
    # If a model is to be compared to the initial model, plot infront,
    # transparent, with overlying lines
    if model_comp:
        val_comp = insp.misfit_values(model_comp)
        plt.hist(x=val_comp, label=model_comp,
                 bins=len(np.arange(-1 * bound, bound, binsize)),
                 histtype="step", edgecolor="deepskyblue", linewidth=2.,
                 alpha=0.7, zorder=100)

    # dln(A) information only relevant in these bounds
    if choice == "dlna":
        plt.xlim([-1.75, 1.75])

    # Finalize plot details
    plt.xlabel(label_dict[choice])
    plt.ylabel("Count")
    plt.title(f"{choice} histogram")
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.grid(linewidth=1, linestyle=":", which="both", zorder=1)
    plt.axvline(x=0, ymin=0, ymax=1, linewidth=1.5, c="k", zorder=2, alpha=0.75,
                linestyle='--')
    plt.legend()

    if save:
        plt.savefig(save)
    if show:
        plt.show()

    return f, ax


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
    vmax = kwargs.get("vmax", 1.1)
    dv = kwargs.get("dv", None)
    label = kwargs.get("label", "")

    norm = mpl.colors.Normalize(vmin=0, vmax=1)
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


def hover_on_plot(f, ax, obj, values, dissapear=False):
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


def annotate_txt(ax, txt, anno_location="lower-right"):
    """
    Convenience function to annotate some information

    :type ax: matplot.axes._subplots.AxesSubplot
    :param ax: axis to annotate onto
    :type txt: str
    :param txt: text to annotate
    :type location: str
    :param location: location on the figure to annotate
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

    ax.annotate(s=txt, xy=(x, y), multialignment=multialignment, fontsize=10)


if __name__ == "__main__":
    pass
