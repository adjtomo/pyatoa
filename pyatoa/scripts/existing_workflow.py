"""
This script is meant to be run when Pyatoa has been run before, and therefore
an HDF5 file already exists with the data. In this case data is read in from
the HDF5 file to be run in the workflow
"""
import os
import glob
import pyasdf
import pyatoa
import logging
import traceback
from pyatoa.utils.tools.io import create_stations_adjoint, \
    write_adj_src_to_ascii

# User defines HDF5 filename
ds_fid = glob.glob(os.path.join(os.getcwd(), '*.h5'))[0]
model_number = "m00"
synthetics_only = False
set_logging = False

# initiate logging
if set_logging:
    logger_pyatoa = logging.getLogger("pyatoa")
    logger_pyatoa.setLevel(logging.DEBUG)
    logger_pyflex = logging.getLogger("pyflex")
    logger_pyflex.setLevel(logging.DEBUG)


# assuming the HDF5 file is named after the event id
event_id = os.path.basename(ds_fid).split('.')[0]
if not os.path.exists("./figures/{}".format(event_id)):
    os.makedirs("./figures/{}".format(event_id))

# instantiate the Pyatoa Config object
config = pyatoa.Config(
    event_id=event_id,
    model_number=model_number,
    min_period=10,
    max_period=30,
    filter_corners=4,
    rotate_to_rtz=False,
    synthetics_only=synthetics_only,
    unit_output="DISP",
    pyflex_map="hikurangi",
    adj_src_type="cc_traveltime_misfit",
    cfgpaths={"waveforms": [], "synthetics": './', "responses": []}
)

# additional text to add to title of waveform plots
append_title = "\npyflex={}; pyadjoint={} ".format(config.pyflex_map,
                                                   config.adj_src_type)

# initiate pyasdf dataset where all data will be saved
with pyasdf.ASDFDataSet(ds_fid) as ds:
    mgmt = pyatoa.Manager(config=config, ds=ds)
    stations = ds.waveforms.list()

    # begin the workflow by looping through all stations
    for station in stations:
        net, sta = station.split('.')
        try:
            mgmt.populate(station_code=f"{net}.{sta}.??.HH?",
                          model=config.model_number
                          )
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
            mgmt.reset()
        except Exception as e:
            traceback.print_exc()
            import ipdb;ipdb.set_trace()
            mgmt.reset()
            continue

    # Create Specfem necessary files
    sem_path = "./SEM"
    if not os.path.exists(sem_path):
        os.makedirs(sem_path)
    write_adj_src_to_ascii(ds, config.model_number, sem_path)
    create_stations_adjoint(ds, config.model_number, sem_path)

