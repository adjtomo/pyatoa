"""
Create waveform plots that show the changes in synthetic waveforms with
progressive model updates
"""
import pyasdf
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm

from pyatoa import Config, Manager
from pyatoa.utils.operations.source_receiver import gcd_and_baz,\
    seismogram_length
from pyatoa.utils.visuals.waveforms import setup_plot



dataset = "/Users/chowbr/Documents/subduction/tomo/seisflows/checkerboard/data/2016p275188.h5"
unit_dict = {"DISP": "displacement [m]",
             "VEL": "velocity [m/s]",
             "ACC": "acceleration [m/s^2]"}
z = 5
dpi = 100
figsize = (11.69, 8.27)


# Read in the dataset
with pyasdf.ASDFDataSet(dataset) as ds:
    event_id = ds.events[0].resource_id.id.split('/')[-1]
    config = Config(
        event_id=event_id,
        model_number="m00",
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

        # Hacky way to preprocess all the synthetic traces
        synthetics = {}
        for syn_tag in ds.waveforms[sta].get_waveform_tags():
            if syn_tag == "observed":
                continue
            mgmt.st_obs = ds.waveforms[sta].observed.copy()
            mgmt.st_syn = ds.waveforms[sta][syn_tag].copy()
            mgmt.preprocess()
            synthetics[syn_tag] = mgmt.st_syn.copy()

        # Collect some more information
        st_obs = mgmt.st_obs.copy()
        length_sec = seismogram_length(
            distance_km=gcd_and_baz(event=ds.events[0], sta=mgmt.inv[0][0])[0],
            slow_wavespeed_km_s=1, binsize=50, minimum_length=100
        )

        # Set some parameters necessary for flexible plotting
        middle_trace = len(st_obs) // 2

        # Instantiate plotting instances
        f = plt.figure(figsize=figsize, dpi=dpi)
        axes, _ = setup_plot(number_of=len(st_obs), twax=False)
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

            # Plot synthetics
            lines_for_legend = [a1]
            synthetic_keys = list(synthetics.keys())
            synthetic_keys.sort()
            viridis = cm.get_cmap('Reds', len(synthetic_keys))
            for j, syn_key in enumerate(synthetic_keys):
                syn = synthetics[syn_key].select(component=comp)
                a2, = axes[i].plot(t, syn[0].data, 'r', zorder=z, label=syn_key,
                                   c=viridis(j)
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
            # Accumulate legend labels
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

        axes[0].set_title(title)
        axes[-1].set_xlabel("time [s]")
        plt.show()
        plt.close()
