"""
An example run script to populate a pyasdf dataset with event, stations,
waveforms, misfit windows and adjoint sources.
Created to be called from within the Seisflows workflow, so command line 
arguments to specify the event id
"""
import os
import sys
import glob
import pyasdf
import pyatoa
import logging
import traceback
import numpy as np
from obspy import read_inventory

from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.operations.file_generation import create_stations_adjoint, \
                                                    write_adj_src_to_ascii

# initiate logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.DEBUG)

# initiate config
model_number = "m00"
event_id_list = ["2018p130600"]
synthetics_only = True
working_dir = os.getcwd()
obs_traces = os.path.join(working_dir, "obs")
syn_traces = os.path.join(working_dir, "syn")

for event_id in event_id_list:
    # make a directory to store generated figures
    fig_dir = os.path.join(os.getcwd(), event_id)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # set the pyatoa config object for misfit quantification 
    config = pyatoa.Config(event_id=event_id, 
                           model_number=model_number, 
                           min_period=10, 
                           max_period=30, 
                           filter_corners=4, 
                           rotate_to_rtz=False, 
                           unit_output="DISP", 
                           window_amplitude_ratio=0.2,
                           pyflex_config="hikurangi_strict", 
                           adj_src_type="cc_hikurangi_strict", 
                           synthetics_only=synthetics_only,
                           cfgpaths={'synthetics':[syn_traces],
                                     'waveforms':[obs_traces],
                                     'responses':[]}
                           )

    # initiate pyasdf dataset where all data will be saved
    ds_path = os.path.join(working_dir, "{}.h5".format(config.event_id))
    with pyasdf.ASDFDataSet(ds_path) as ds:
        clean_ds(ds)
        config.write_to_asdf(ds)
        stations = np.loadtxt(os.path.join(working_dir, "STATIONS"), 
                              usecols=[0, 1], dtype=str)

        # begin the Pyatoa Workflow, loop through all stations
        mgmt = pyatoa.Manager(config=config, ds=ds)
        for station in stations: 
            sta, net = station
            try:
                mgmt.reset()
                mgmt.gather_data(station_code="{net}.{sta}.{loc}.{cha}".format(
                                 net=net, sta=sta, loc="*", cha="HH?")
                )
                mgmt.preprocess()
                mgmt.run_pyflex()
                mgmt.run_pyadjoint()
                mgmt.plot_wav(save=os.path.join(fig_dir, "wav_{sta}".format(
                    sta=sta)), show=False
                )
                mgmt.plot_map(save=os.path.join(fig_dir, "map_{sta}".format(
                    sta=sta)), show=False
                )
            except Exception as e:
                print("\n")
                traceback.print_exc()
                mgmt.reset()
                continue

        sem_path = os.path.join(working_dir, "SEM")
        write_adj_src_to_ascii(ds, model_number, sem_path)
        create_stations_adjoint(ds, model_number, 
                                os.path.join(working_dir, "STATIONS"), sem_path)


