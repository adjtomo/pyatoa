"""
This script is meant to be used when no data has been collected (aside from
Synthetic waveforms generated from Specfem). Pyatoa will instantiate a new
HDF5 object and fill it up if possible. The workflow will complete by creating
ascii files for the adjoint sources, and a STATIONS_ADJOINT file, which are
necessary for a Specfem3D adjoint simulation.
"""
import os
import pyasdf
import pyatoa
import logging
import traceback
import numpy as np
from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.io import create_stations_adjoint, write_adj_src_to_ascii

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
    if not os.path.exists("./figures/{}".format(event_id)):
        os.makedirs("./figures/{}".format(event_id))

    # instantiate the Pyatoa Config object
    config = pyatoa.Config(
        event_id=event_id,
        model_number="m00",
        min_period=10,
        max_period=30,
        filter_corners=4,
        rotate_to_rtz=False,
        synthetics_only=synthetics_only,
        unit_output="DISP",
        pyflex_preset="default",
        adj_src_type="cc_traveltime_misfit",
        cfgpaths={"waveforms": [], "synthetics": './', "responses": []}
    )

    # additional text to add to title of waveform plots
    append_title = "pyflex={}; pyadjoint={} ".format(config.pyflex_preset,
                                                     config.adj_src_type)

    # initiate pyasdf dataset where all data will be saved
    complete = 0
    with pyasdf.ASDFDataSet("./{}.h5".format(config.event_id)) as ds:
        # clean the dataset for a new workflow
        clean_ds(ds)
        config.write_to_asdf(ds)

        # begin the workflow by looping through all stations
        mgmt = pyatoa.Manager(config=config, ds=ds)
        for station in stations:
            sta, net = station
            try:
                mgmt.gather_data(station_code="{net}.{sta}.{loc}.{cha}".format(
                                 net=net, sta=sta, loc="*", cha="HH?"))

                mgmt.preprocess()
                mgmt.run_pyflex()
                mgmt.run_pyadjoint()
                mgmt.plot_wav(save="./figures/{eid}/wav_{sta}".format(
                    eid=config.event_id, sta=sta), append_title=append_title,
                    show=True,
                )
                mgmt.plot_map(save="./figures/{eid}/map_{sta}".format(
                    eid=config.event_id, sta=sta), show=True
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
            write_adj_src_to_ascii(ds, config.model_number, sem_path)
            create_stations_adjoint(ds, config.model_number, station_file,
                                    sem_path)


