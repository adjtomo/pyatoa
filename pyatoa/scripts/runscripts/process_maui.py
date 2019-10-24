"""
An example run script to populate a pyasdf dataset with event, stations,
waveforms, misfit windows and adjoint sources.
Created to be called from within the Seisflows workflow, so command line 
arguments to specify the event id
"""
import os

import pyasdf
import pyatoa
import logging
import traceback
import numpy as np

from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.tools.io import create_stations_adjoint, \
    write_adj_src_to_ascii

# initiate logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.DEBUG)

# initiate config
model_number = "m00"
event_id_list = ["2013p142607"]
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
                                     'responses':["/scale_wlg_persistent/filesets/home/chowbr/primer/auxiliary/seed/RESPONSE"]
                                     }
                           )

    # append some information to the title
    append_title = "pyflex={}, pyadjoint={}".format(config.pyflex_config[0],
                                                    config.pyadjoint_config[0])

    # initiate pyasdf dataset where all data will be saved
    ds_path = os.path.join(working_dir, "{}.h5".format(config.event_id))
    with pyasdf.ASDFDataSet(ds_path) as ds:
        clean_ds(ds)
        config.write_to_asdf(ds)
        
        # begin the Pyatoa Workflow, loop through all stations
        mgmt = pyatoa.Manager(config=config, ds=ds)

        # load in station information
        stations = np.loadtxt(os.path.join(working_dir, "STATIONS"), 
                              usecols=[0, 1], dtype=str)
        empties = 0
        for station in stations: 
            sta, net = station
            try:
                mgmt.reset()
                mgmt.gather_data(station_code="{net}.{sta}.{loc}.{cha}".format(
                                 net=net, sta=sta, loc="*", cha="HH?")
                )
                mgmt.preprocess()
                mgmt.run_pyflex()

                # count all stations w/ no windows
                if not mgmt.windows:
                    empties += 1

                mgmt.run_pyadjoint()
                mgmt.plot_wav(save=os.path.join(fig_dir, "wav_{sta}".format(
                    sta=sta)), append_title=append_title, show=False
                )
                mgmt.plot_map(save=os.path.join(fig_dir, "map_{sta}".format(
                    sta=sta)), show=False
                )
            except Exception as e:
                print("\n")
                traceback.print_exc()
                mgmt.reset()
                continue

        # path for the adjoint sources
        sem_path = os.path.join(working_dir, "SEM")
        if not os.path.exists(sem_path):
            os.makedirs(sem_path)

        write_adj_src_to_ascii(ds, model_number, sem_path)
        create_stations_adjoint(ds, model_number, 
                                os.path.join(working_dir, "STATIONS"), sem_path)

        print("{}/{} STATIONS W >1 WINDOWS".format(len(stations), empties))
