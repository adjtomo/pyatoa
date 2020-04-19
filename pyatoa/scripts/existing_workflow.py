"""
This script is meant to be run when Pyatoa has been run before, and therefore
an HDF5 file already exists with the data. In this case data is read in from
the HDF5 file to be run in the workflow
"""
import os
import glob
import pyatoa
import logging
import traceback
from pyasdf import ASDFDataSet as asdf

# User defines HDF5 filename
ds_fid = glob.glob(os.path.join(os.getcwd(), '*.h5'))[0]
model = "m00"
step = "s00"
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
if not os.path.exists(f"./figures/{event_id}"):
    os.makedirs(f"./figures/{event_id}")

# additional text to add to title of waveform plots
append_title = (f"pyflex={config.pyflex_preset}; "
                f"pyadjoint={config.adj_src_type} ")

# initiate pyasdf dataset where all data will be saved
with asdf(ds_fid) as ds:
    cfg = pyatoa.Config(ds=ds, path=f"{model}{step}")
    mgmt = pyatoa.Manager(config=cfg, ds=ds)
    stations = ds.waveforms.list()

    # begin the workflow by looping through all stations
    for station in stations:
        net, sta = station.split('.')
        try:
            mgmt.load(station)
            mgmt.standardize()
            mgmt.preprocess()
            mgmt.window()
            mgmt.measure()
            mgmt.plot_wav(save="./figures/{config.event_id}/wav_{sta}",
                          append_title=append_title, show=True,
                          )
            mgmt.plot_map(save="./figures/{config.event_id}/map_{sta}",
                          show=True
                          )
            mgmt.reset()
        except Exception as e:
            traceback.print_exc()
            import ipdb;ipdb.set_trace()
            mgmt.reset()
            continue


