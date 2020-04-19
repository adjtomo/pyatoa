"""
This script is meant to be used when no data has been collected (aside from
Synthetic waveforms generated from Specfem). Pyatoa will instantiate a new
HDF5 object and fill it up if possible. The workflow will complete by creating
ascii files for the adjoint sources, and a STATIONS_ADJOINT file, which are
necessary for a Specfem3D adjoint simulation.
"""
import os
import pyatoa
import logging
import traceback
import numpy as np
from pyasdf import ASDFDataSet as asdf
from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.write import create_stations_adjoint, write_adj_src_to_ascii

# User defined event id(s) and location of STATION file from Specfem3D
event_ids = [""]
station_file = "./STATIONS"
synthetics_only = True

# initiate logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.DEBUG)

stations = np.loadtxt(station_file, usecols=[0, 1], dtype=str)

# If only one station present in file, ensure that the station loop still works
if len(stations) == 2:
    stations = [stations]

for event_id in event_ids:
    if not os.path.exists(f"./figures/{event_id}"):
        os.makedirs(f"./figures/{event_id}")

    # instantiate the Pyatoa Config object
    config = pyatoa.Config(
        event_id=event_id,
        model="m00",
        min_period=10,
        max_period=30,
        filter_corners=4,
        rotate_to_rtz=False,
        synthetics_only=synthetics_only,
        unit_output="DISP",
        pyflex_preset="default",
        adj_src_type="cc",
        cfgpaths={"waveforms": [], "synthetics": './', "responses": []}
    )

    # additional text to add to title of waveform plots
    append_title = (f"pyflex={config.pyflex_preset}; "
                    f"pyadjoint={config.adj_src_type} ")

    # initiate pyasdf dataset where all data will be saved
    complete = 0
    with asdf(f"./{config.event_id}.h5") as ds:
        # clean the dataset for a new workflow
        clean_ds(ds)
        config.write(ds)

        # begin the workflow by looping through all stations
        mgmt = pyatoa.Manager(config=config, ds=ds)
        for station in stations:
            sta, net = station
            try:
                mgmt.gather(station_code=f"{net}.{sta}.*.HH?")
                mgmt.standardize()
                mgmt.preprocess()
                mgmt.window()
                mgmt.measure()
                mgmt.plot_wav(save=f"./figures/{config.event_id}/wav_{sta}",
                              append_title=append_title,
                              show=True,
                              )
                mgmt.plot_map(save=f"./figures/{config.event_id}/map_{sta}",
                              show=True
                              )

                if mgmt.num_windows:
                    complete += 1
                mgmt.reset()
            except Exception as e:
                traceback.print_exc()
                import ipdb;ipdb.set_trace()
                mgmt.reset()
                continue

        # Create Specfem necessary files
        if complete:
            sem_path = "./SEM"
            if not os.path.exists(sem_path):
                os.makedirs(sem_path)
            write_adj_src_to_ascii(ds=ds, model=config.model, pathout=sem_path)
            create_stations_adjoint(ds=ds, model=config.model,
                                    specfem_station_file=station_file,
                                    pathout=sem_path)


