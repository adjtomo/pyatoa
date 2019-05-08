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
from obspy import read_inventory

from pyatoa.utils.operations.pyasdf_editing import clean_ds
from pyatoa.utils.operations.formatting import write_adj_src_to_ascii
from pyatoa.utils.operations.file_generation import create_stations_adjoint

# initiate logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.DEBUG)

# initiate config
model_number = "m00"
basepath = "/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/"
event_id_list = [os.path.basename(_) for _ in glob.glob(os.path.join(
                              os.getcwd(), "test_synthetics", model_number, "*")
                  )]
working_dir = os.getcwd()

for event_id in event_id_list:
    # make a directory to store generated figures
    fig_dir = os.path.join(os.getcwd(), event_id)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # set the pyatoa config object for misfit quantification 
    config = pyatoa.Config(
        event_id=event_id, model_number=model_number, min_period=10, 
        max_period=30, filter_corners=4, rotate_to_rtz=False, 
        unit_output="DISP", pyflex_config="UAF", 
        adj_src_type="multitaper_misfit", 
        paths_to_synthetics=[os.path.join(os.getcwd(), "test_synthetics")],
        # paths_to_waveforms=[os.path.join(basepath, "primer", "seismic"],
        # paths_to_responses=[os.path.join(basepath, "primer", "seed", "RESPONSE"]
        )

    # initiate pyasdf dataset where all data will be saved
    ds = pyasdf.ASDFDataSet(os.path.join(
                            os.getcwd(), "hdf5", "{}.h5".format(config.event_id))
                           )
    clean_ds(ds)
    config.write_to_asdf(ds)

    # begin the Pyatoa Workflow, loop through all stations
    mgmt = pyatoa.Manager(config=config, ds=ds)
    master_inventory = read_inventory(os.path.join(basepath, "primer", 
                                                   "auxiliary_data", "stationxml", 
                                                   "master_inventory_slim.xml")
                                     )
    for net in master_inventory:
        for sta in net:
            if sta.is_active(time=mgmt.event.preferred_origin().time):
                try:
                    mgmt.gather_data(
                        station_code="{net}.{sta}.{loc}.{cha}".format(
                            net=net.code, sta=sta.code, loc="*", cha="HH?")
                    )
                    mgmt.preprocess()
                    mgmt.run_pyflex()
                    mgmt.run_pyadjoint()
                    mgmt.plot_wav(save=os.path.join(fig_dir, "wav_{sta}".format(
                        sta=sta.code)), show=False
                    )
                    mgmt.plot_map(save=os.path.join(fig_dir, "map_{sta}".format(
                        sta=sta.code)), show=False
                    )
                    mgmt.reset()
                except Exception as e:
                    print("\n")
                    traceback.print_exc()
                    mgmt.reset()
                    continue

    sem_path = os.path.join(os.getcwd(), "SEM")
    write_adj_src_to_ascii(ds, model_number=model_number, filepath=sem_path)
    create_stations_adjoint(ds, model_number=model_number, filepath=sem_path)


