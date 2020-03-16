"""
Create waveform plots that show the changes in synthetic waveforms with
progressive model updates. Each individual model gets its on row in the plot
"""
import sys
import os
import glob
import pyasdf
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
  
from obspy.signal.cross_correlation import correlate, xcorr_max
from pyatoa import Config, Manager
from pyatoa.utils.asdf.extractions import windows_from_ds
from pyatoa.utils.tools.srcrcv import gcd_and_baz, seismogram_length

mpl.rcParams['lines.linewidth'] = 1.8
mpl.rcParams['axes.linewidth'] = 2.
mpl.rcParams['font.size'] = 14

warnings.simplefilter("ignore")


def setup_plot(nrows, ncols, label_y_units=False, label_x_units=None):
    """
    Dynamically set up plots according to number_of given
    Returns a list of lists of axes objects
    e.g. axes[i][j] gives the ith column and the jth row

    :type nrows: int
    :param nrows: number of rows in the gridspec
    :type ncols: int
    :param ncols: number of columns in the gridspec
    :rtype axes: matplotlib axes
    :return axes: axis objects
    """
    gs = mpl.gridspec.GridSpec(nrows, ncols, hspace=0, wspace=0.05, 
                               width_ratios=[2] * ncols,
                               height_ratios=[1] * nrows
                               )

    axes = []
    for row in range(0, gs.get_geometry()[0]):
        components = []
        for col in range(0, gs.get_geometry()[1]):
            if row == 0:
                ax = plt.subplot(gs[row, col])
            else:
                # Share the x axis with the same component in the column
                ax = plt.subplot(gs[row, col], sharex=axes[0][col],
                                 sharey=axes[0][col]
                                 )
            ax.set_axisbelow(True)
            ax.tick_params(which='both', direction='in', top=True,
                                 right=True)
            # Set the grids on
            ax.minorticks_on()
            for axis_ in ['major', 'minor']:
                ax.grid(which=axis_, linestyle=':', linewidth='0.5',
                              color='k', alpha=0.75)
            components.append(ax)
        axes.append(components)

    # remove x-tick labels except for last axis, middle row
    for row in axes[:-1]:
        for col in row:
            plt.setp(col.get_xticklabels(), visible=False)
    if label_x_units:
        for i, col in enumerate(axes[-1]):
            if i not in label_x_units:
                plt.setp(col.get_xticklabels(), visible=False)

    # remove y-tick labels except for first axis
    if label_y_units:
        for i, row in enumerate(axes):
            for j, col in enumerate(row):
                if j != 0: 
                    plt.setp(col.get_yticklabels(), visible=False)
                else:    
                    col.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    else:
        for i, row in enumerate(axes):
            for col in row:
                plt.setp(col.get_yticklabels(), visible=False)

    return axes


def center_on_peak_energy(st, thresh_value=0.1):
    """
    An attempt to determine the extent of a peak energy in a seismogram so that
    plots do not show unnecessary tails, otherwise for small plots, a lot of 
    information gets lost
    """
    # Set the starting indices at values that will get easily overwritten
    start_index = len(st[0].data)
    end_index = 0
    for tr in st:
        threshold = tr.data.max() * thresh_value
        indices = np.where(tr.data > threshold)[0]
        if indices[0] < start_index:
            start_index = indices[0]
        if indices[-1] > end_index:
            end_index = indices[-1]

    return start_index, end_index


def plot_iterative_waveforms(datasets_path, output_dir, min_period, max_period,
                             event_id="*.h5", select_models=[], 
                             components=["Z", "N", "E"], figsize=(11.69, 8.27),
                             select_stations=[], synthetics_only=False,
                             trace_length=[], cross_corr=False, 
                             time_shift=False, label_y_units=False, 
                             label_x_units=None, show=True):
    """
    Main function to plot waveforms iterative based on model updates

    :type datasets_path: str
    :param datasets_path: path to directory containing pyasdf datasets
    :type output_dir: str
    :param output_dir: path to save the figure to
    :type min_period: float
    :param min_period: minimum filter period for waveforms
    :type max_period: float
    :param max_period: maximum filter period for waveforms
    :type event_id: str
    :param event_id: event id to identify the asdf dataset, can be wildcard
    :type select_models: list of str
    :param select_models: allows User to select which models are included in the
        waveform plots
    :type select_stations: list of str
    :param select_stations: allows User to select which stations queried
        for the waveform plots
    :type synthetics_only: bool
    :param synthetics_only: synthetics only parameter to be passed to the
        Config. if False, instrument response will be removed from obs.
    :type trace_length: list of floats
    :param trace_length: [trace_start, trace_end] will be used to set the x
        limit on the waveform data. If none, no xlim will be set
    :type cross_corr: bool
    :param cross_corr: if True, cross correlation will be taken between traces
        at each model and max correlation will be annotated to each plot.
    :type label_units: bool
    :param label_units: if True, units will be placed on the y-axis
    :type show: bool
    :param show: Show the plot or do not
    """
    # Plotting parameters
    dpi = 125
    z = 5
    comp_colors = ["r", "b", "g"]
    window_colors = ["orange", "yellow", "lawngreen", "cyan",]
    unit_dict = {"DISP": "displacement [m]",
                 "VEL": "velocity [m/s]",
                 "ACC": "acceleration [m/s^2]"}

    # Set up a dictionary to collect peak amplitude information
    max_amp = {}
    for comp in components:
        max_amp[comp] = 0

    # Read in each of the datasets
    for dataset in glob.glob(os.path.join(datasets_path, event_id)):
        with pyasdf.ASDFDataSet(dataset) as ds:
            event_id = ds.events[0].resource_id.id.split('/')[-1]

            # User defined filtering parameters
            config = Config(
                event_id=event_id,
                model_number=0,
                min_period=min_period,
                max_period=max_period,
                filter_corners=4,
                component_list=components,
                rotate_to_rtz=False,
                unit_output="DISP",
                synthetics_only=synthetics_only
            )

            # Loop through the available stations
            for sta in ds.waveforms.list():
                if select_stations and (sta not in select_stations):
                    continue 
                print(sta)
                mgmt = Manager(config=config, empty=True)
                mgmt.inv = ds.waveforms[sta].StationXML
                mgmt.event = ds.events[0]
            
                # Hacky way to preprocess all the synthetic traces using Pyatoa
                synthetics, windows = {}, {}
                for syn_tag in ds.waveforms[sta].get_waveform_tags():
                    if syn_tag == "observed":
                        continue
                    else:
                        if select_models and (syn_tag not in select_models):
                            continue
                        mgmt.st_obs = ds.waveforms[sta].observed.copy()
                        mgmt.st_syn = ds.waveforms[sta][syn_tag].copy()
                        mgmt.standardize()
                        mgmt.preprocess()
                        synthetics[syn_tag] = mgmt.st_syn.copy()
                        windows[syn_tag], _ = windows_from_ds(
                                               ds, model=syn_tag.split('_')[-1],
                                               net=sta.split('.')[0],
                                               sta=sta.split('.')[1])

                # Collect the observed trace last
                st_obs = mgmt.st_obs.copy()

                # Determine a rough estimate for the length of the seismogram
                length_sec = seismogram_length(
                    distance_km=gcd_and_baz(event=ds.events[0],
                                            sta=mgmt.inv[0][0])[0],
                    slow_wavespeed_km_s=1, binsize=50, minimum_length=100
                )

                # Instantiate plotting instances
                f = plt.figure(figsize=figsize)
                axes = setup_plot(nrows=len(synthetics.keys()), 
                                  ncols=len(components), 
                                  label_y_units=label_y_units,
                                  label_x_units=label_x_units
                                  )
                middle_column = len(components) // 2

                # Create time axis based on data statistics
                t = st_obs[0].times(
                       reftime=st_obs[0].stats.starttime - mgmt.time_offset_sec)

                # Make sure the models are in order
                synthetic_keys = list(synthetics.keys())
                synthetic_keys.sort()
                synthetic_init = synthetics[synthetic_keys[0]]

                # Plot each model on a different row
                for row, syn_key in enumerate(synthetic_keys):
                    ylab = ""
                    if len(select_models) == 2:
                        if syn_key == "synthetic_m00":
                            ylab += "initial model"
                        else:
                            ylab += "updated model"
                    # this is the model number e.g. "m00"
                    else:
                        ylab += syn_key.split('_')[-1]
                    if label_y_units:
                        ylab += f"\n{unit_dict[config.unit_output]}"
                
                    # Plot each component in a different column
                    for col, comp in enumerate(config.component_list):
                        obs = st_obs.select(component=comp)
                        syn = synthetics[syn_key].select(component=comp)

                        # Plot waveform
                        a1, = axes[row][col].plot(t, obs[0].data, 'k', zorder=z,
                                                  label="True")
                        a2, = axes[row][col].plot(t, syn[0].data, 
                                                  comp_colors[col], zorder=z,
                                                  label="Syn")

                        # Use the max amplitude to set the ylimits
                        max_amp_ = max([abs(obs.max()[0]),abs(syn.max()[0])])
                        if max_amp_ > max_amp[comp]:
                            max_amp[comp] = max_amp_
                            ylim = max_amp[comp] + 0.1 * max_amp[comp]
                            axes[row][col].set_ylim([-1 * ylim, ylim])

                        # Min and max are now set based on waveform plots
                        xmin, xmax = axes[row][col].get_xlim()
                        ymin, ymax = axes[row][col].get_ylim()

                        # Plot windows, if no windows move on
                        try:
                            window_list = windows[syn_key][comp]
                            for i, win in enumerate(window_list):
                                tleft = (win.left * win.dt + 
                                                mgmt.time_offset_sec)
                                tright = (win.right * win.dt + 
                                                mgmt.time_offset_sec)

                                axes[row][col].add_patch(
                                        mpl.patches.Rectangle(
                                              xy=(tleft, -1 * max_amp[comp]), 
                                              width=tright-tleft, 
                                              height=2.5*max_amp[comp], 
                                              color=window_colors[i],
                                              edgecolor="k",
                                              alpha=(win.max_cc_value **2) / 4)
                                                        )
                                # Annotate time shift information for each win
                                if time_shift:
                                    anno_y = .75 * (ymax - ymin) + ymin
                                    axes[row][col].annotate(
                                            f"{win.cc_shift*win.dt:.2f}s",
                                            xy=(tleft, anno_y), fontsize=8,
                                            zorder=z+5)
                        except KeyError:
                            pass

                        # cross-correlate obs and syn for a sense of misfit
                        if cross_corr:
                            cc = correlate(obs[0].data, syn[0].data, 
                                           shift=len(obs[0].data), 
                                           domain="freq")
                            f_shift, corr = xcorr_max(cc)
                            axes[row][col].annotate(f"cc:{corr:.2f}", xy=(1,0),
                                                    xycoords="axes fraction",
                                                    horizontalalignment="right",
                                                    verticalalignment="bottom",
                                                    weight="bold")

                        if row == 0:
                            # determine how long the traces should be
                            # hardcode the trace length based on user params
                            if isinstance(trace_length, list):
                                if len(trace_length) == 2:
                                    axes[row][col].set_xlim(trace_length)
                                else:
                                    start, end, step = trace_length
                                    axes[row][col].set_xticks(
                                            np.arange(start, end, step)
                                            )
                                    axes[row][col].set_xlim(trace_length[:2])
                            # center trace on the peak energy
                            elif trace_length == "center_on_peak":
                                st_idx, end_idx = center_on_peak_energy(
                                                                      obs + syn)
                                t_start = max(t[st_idx] - 10, 0)
                                t_end = min(t[end_idx] + 10, t[-1])

                                axes[row][col].set_xlim([t_start, t_end])                            
                            # determine the length of the seismogram based on
                            # arrival of latest energy
                            else:
                                if not length_sec:
                                    length_sec = t[-1]
                                axes[row][col].set_xlim(
                                    [np.maximum(mgmt.time_offset_sec, -10),
                                     np.minimum(length_sec, t[-1])
                                     ])

                            # Set title for the first row
                            if col == middle_column:
                                axes[row][col].set_title(
                                        f"{st_obs[0].stats.network}."
                                        f"{st_obs[0].stats.station} "
                                        f"{event_id}")

                            # Annotate component to the top corner
                            axes[row][col].annotate(f"{comp[-1]}", xy=(1,0),
                                                    xycoords="axes fraction",
                                                    horizontalalignment="right",
                                                    verticalalignment="bottom")

                    
                    # y-label after all the processing has occurred
                    axes[row][0].set_ylabel(ylab)

                # Label the time axis on the bottom row
                axes[-1][middle_column].set_xlabel("time [sec]")

                # Save the generated figure
                if output_dir:
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)

                    final_model = synthetic_keys[-1].split('_')[-1]
                    fid_out = os.path.join(
                        output_dir, f"{event_id}_{sta}_{final_model}.png")
                    plt.savefig(fid_out, figsize=figsize, dpi=dpi)

                # Show the plot
                if show:
                    plt.show()
                
                plt.close()


if __name__ == "__main__":
    try:
        # Choose
        i = 2

        # Set parameters here 
        datasets_path = "./data"

        # Path to save figures to, if None given, figures wil not be saved
        output_dir = "./waveforms"
        
        # If you only want to choose one event in your directory, wildcards okay
        event_ids = ["2017p059122.h5", "2013p142607.h5", "2019p304574.h5",
                     "2019p738432.h5"]  
        event_id = event_ids[i]

        # If you don't want to plot all models, can add e.g. 'synthetic_m00' 
        select_models = []

        # Pick stations, if left empty, will plot all stations in dataset
        select_stations_list = [['NZ.WAZ'], ['NZ.WEL'], ['NZ.KNZ'], ['NZ.MKAZ']]
        select_stations = select_stations_list[i]

        # list of two ints, list of three ints, "dynamic" (default) or 
        # "center_on_peak"
        trace_lengths = [[55, 200], [40, 250], [40, 160], [45, 170]]
        trace_length = trace_lengths[i]

        # Synthetic only tests need to be treated differently
        min_period = 10
        max_period = 30
        synthetics_only = True

        # User-defined figure parameters
        figsize = (7, 10)
        label_y_units = False  # label the units of the traces, otherwise blank
        label_x_units = None    # index of column to label x units, otherwise all
        cross_corr = False  # cross-correlate traces and annotate max correlation
        time_shift = True  # put time shift somewhere in the plot
        components = ["Z", "N", "E"]
        show = True


        plot_iterative_waveforms(datasets_path=datasets_path, 
                                 output_dir=output_dir, 
                                 min_period=min_period, max_period=max_period,
                                 event_id=event_id, select_models=select_models, 
                                 components=components,
                                 select_stations=select_stations, 
                                 synthetics_only=synthetics_only,
                                 trace_length=trace_length, 
                                 cross_corr=cross_corr, 
                                 time_shift=time_shift, 
                                 figsize=figsize,
                                 label_y_units=label_y_units, 
                                 label_x_units=label_x_units, show=show)
    except KeyboardInterrupt:
        sys.exit()
