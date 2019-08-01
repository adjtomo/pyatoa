"""
Create waveform plots that show the changes in synthetic waveforms with
progressive model updates. All waveforms superimposed on one another
"""
import os
import glob
import pyasdf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyatoa import Config, Manager
from pyatoa.utils.operations.source_receiver import gcd_and_baz,\
    seismogram_length

mpl.rcParams['lines.linewidth'] = 1.1


def setup_plot(nrows, ncols):
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
    gs = mpl.gridspec.GridSpec(nrows, ncols, hspace=0, wspace=0.1)

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

    # remove x-tick labels except for last axis
    for row in axes[:-1]:
        for column in row:
            plt.setp(column.get_xticklabels(), visible=False)

    # # remove y-tick labels expcet for first axis
    # middle_row = nrows // 2
    # for i, row in enumerate(axes):
    #     if i != middle_row:
    #         for column in row:
    #             plt.setp(column.get_yticklabels(), visible=False)

    # remove all y-tick labels
    for i, row in enumerate(axes):
        for column in row:
            plt.setp(column.get_yticklabels(), visible=False)

    return axes


# Path to hdf5 files
datasets_path = "./"

# Path to save figures to, if none given, figures wil not be saved
output_dir = "./waveforms"

# User-defined figure parameters
show_plots = False
dpi = 125
figsize = (11.69, 8.27)

# For use when plotting labels
z = 5
color_list = ["r", "b", "g"]
unit_dict = {"DISP": "displacement [m]",
             "VEL": "velocity [m/s]",
             "ACC": "acceleration [m/s^2]"}


# Read in each of the datasets
for dataset in glob.glob(os.path.join(datasets_path, '*')):
    with pyasdf.ASDFDataSet(dataset) as ds:
        event_id = ds.events[0].resource_id.id.split('/')[-1]

        # User defined filtering parameters
        config = Config(
            event_id=event_id,
            model_number=0,
            min_period=10,
            max_period=30,
            filter_corners=4,
            rotate_to_rtz=False,
            unit_output="DISP",
            synthetics_only=True
        )

        # Loop through the available stations
        for sta in ds.waveforms.list():
            print(sta)
            mgmt = Manager(config=config, empty=True)
            mgmt.inv = ds.waveforms[sta].StationXML
            mgmt.event = ds.events[0]

            # Hacky way to preprocess all the synthetic traces using Pyatoa
            synthetics = {}
            for syn_tag in ds.waveforms[sta].get_waveform_tags():
                if syn_tag == "observed":
                    continue
                mgmt.st_obs = ds.waveforms[sta].observed.copy()
                mgmt.st_syn = ds.waveforms[sta][syn_tag].copy()
                mgmt.preprocess()
                synthetics[syn_tag] = mgmt.st_syn.copy()

            # Collect the observed trace last
            st_obs = mgmt.st_obs.copy()

            # Determine a rough estimate for the length of the seismogram
            length_sec = seismogram_length(
                distance_km=gcd_and_baz(event=ds.events[0],
                                        sta=mgmt.inv[0][0])[0],
                slow_wavespeed_km_s=1, binsize=50, minimum_length=100
            )

            # Instantiate plotting instances
            f = plt.figure(figsize=figsize, dpi=dpi)
            axes = setup_plot(nrows=len(synthetics.keys()), ncols=len(st_obs))
            middle_column = len(st_obs) // 2
            middle_row = len(synthetics.keys()) // 2

            # Create time axis based on data statistics
            t = np.linspace(
                mgmt.time_offset_sec,
                (st_obs[0].stats.endtime - st_obs[0].stats.starttime +
                 mgmt.time_offset_sec),
                len(st_obs[0].data)
            )

            # Make sure the models are in order
            synthetic_keys = list(synthetics.keys())
            synthetic_keys.sort()
            synthetic_init = synthetics[synthetic_keys[0]]

            # Plot each model on a different row
            for row, syn_key in enumerate(synthetic_keys):
                axes[row][0].set_ylabel("{}".format(syn_key.split('_')[-1]))

                # This will put units on the ylabel
                # if row == middle_row:
                #     # Label the leftmost column by the model number and units
                #     axes[row][0].set_ylabel(
                #         "{}\n{}".format(unit_dict[config.unit_output],
                #                         syn_key.split('_')[-1]
                #                         )
                #     )
                # else:
                #     axes[row][0].set_ylabel("{}".format(syn_key.split('_')[-1]))

                # Plot each component in a different column
                for col, comp in enumerate(config.component_list):
                    obs = st_obs.select(component=comp)
                    syn = synthetics[syn_key].select(component=comp)
                    syn_init = synthetic_init.select(component=comp)

                    # Plot waveform
                    a1, = axes[row][col].plot(t, obs[0].data, 'k', zorder=z)
                    a2, = axes[row][col].plot(t, syn[0].data, color_list[col],
                                              zorder=z)
                    a3, = axes[row][col].plot(t, syn_init[0].data, color_list[col],
                            zorder=z-1, alpha=0.35, linestyle=':', linewidth=2.)

                    if row == 0:
                        # Set the seismogram length but only for the first row
                        axes[row][col].set_xlim([25,100])
                        # if not length_sec:
                        #     length_sec = t[-1]
                        # axes[row][col].set_xlim(
                        #     [np.maximum(mgmt.time_offset_sec, -10),
                        #      np.minimum(length_sec, t[-1])
                        #      ])

                        # Set titles for the first row
                        if col == middle_column:
                            title = "{net}.{sta} {eid}\n".format(
                                net=st_obs[0].stats.network,
                                sta=st_obs[0].stats.station,
                                eid=event_id
                            )
                            title += comp
                        else:
                            title = comp
                        axes[row][col].set_title(title)

            # Label the time axis on the bottom row
            axes[-1][middle_column].set_xlabel("time [sec]")

            # Save the generated figure
            if output_dir:
                final_model = synthetic_keys[-1].split('_')[-1]
                fid_out = os.path.join(
                    output_dir, "{eid}_{sta}_{mod}.png".format(eid=event_id,
                                                               sta=sta,
                                                               mod=final_model)
                )
                plt.savefig(fid_out, figsize=figsize, dpi=dpi)

            # Show the plot
            if show_plots:
                plt.show()

            plt.close()
