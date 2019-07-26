"""
Create waveform plots that show the changes in synthetic waveforms with
progressive model updates. All waveforms superimposed on one another
"""
import os
import glob
import pyasdf
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm

from pyatoa import Config, Manager
from pyatoa.utils.operations.source_receiver import gcd_and_baz,\
    seismogram_length
from pyatoa.utils.visuals.waveforms import setup_plot

# Path to hdf5 files
datasets_path = "path/to/dataset"

# Path to save figures to, if none given, figures wil not be saved
output_dir = "path/to/output"

# User-defined figure parameters
show_plots = True
colormap_choice = 'viridis'
dpi = 100
figsize = (11.69, 8.27)

# For use when plotting labels
z = 5
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

            # Set some parameters necessary for flexible plotting
            middle_trace = len(st_obs) // 2

            # Instantiate plotting instances
            f = plt.figure(figsize=figsize, dpi=dpi)
            axes, _ = setup_plot(number_of=len(st_obs), twax=False)

            # Create time axis based on data statistics
            t = np.linspace(
                mgmt.time_offset_sec,
                (st_obs[0].stats.endtime - st_obs[0].stats.starttime +
                 mgmt.time_offset_sec),
                len(st_obs[0].data)
            )

            # Plot per component
            for i, comp in enumerate(config.component_list):
                # Plot observation waveform
                obs = st_obs.select(component=comp)
                a1, = axes[i].plot(t, obs[0].data, 'k', zorder=z,
                                   label="{} (OBS)".format(obs[0].get_id()))
                lines_for_legend = [a1]

                # Plot each synthetic waveforms
                synthetic_keys = list(synthetics.keys())
                synthetic_keys.sort()
                colormap = cm.get_cmap(colormap_choice, len(synthetic_keys))
                for j, syn_key in enumerate(synthetic_keys):
                    syn = synthetics[syn_key].select(component=comp)
                    a2, = axes[i].plot(
                        t, syn[0].data, 'r', zorder=z, label=syn_key,
                        c=colormap(j)
                    )
                    lines_for_legend.append(a2)

                # Set the seismogram length
                if not length_sec:
                    length_sec = t[-1]
                axes[i].set_xlim([np.maximum(mgmt.time_offset_sec, -10),
                                  np.minimum(length_sec, t[-1])
                                  ])

                # The y-label of the middle trace contains common info
                if i == middle_trace:
                    comp = "{}\n{}".format(unit_dict[config.unit_output], comp)
                axes[i].set_ylabel(comp)

                # Accumulate legend labels and create a legend
                labels = [l.get_label() for l in lines_for_legend]
                axes[i].legend(lines_for_legend, labels, prop={"size": 9},
                               loc="upper right")

            # Set plot title with relevant information
            title = "{net}.{sta}".format(net=st_obs[0].stats.network,
                                         sta=st_obs[0].stats.station
                                         )

            # Account for the fact that event id may not be given
            if config.event_id is not None:
                title = " ".join([title, config.event_id])

            # Final plotting parameters
            axes[0].set_title(title)
            axes[-1].set_xlabel("time [s]")

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
