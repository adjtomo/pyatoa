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
from pyatoa.utils.form import format_event_name

# User defines HDF5 filename
ds_fid = "2017p819775.h5"
model = "m00"
step = "s00"
synthetics_only = False
set_logging = True

# initiate logging
if set_logging:
    logger_pyatoa = logging.getLogger("pyatoa")
    logger_pyatoa.setLevel(logging.DEBUG)
    logger_pyflex = logging.getLogger("pyflex")
    logger_pyflex.setLevel(logging.DEBUG)


# initiate pyasdf dataset where all data will be saved
with asdf(ds_fid) as ds:
    event_id = format_event_name(ds)
    if not os.path.exists(f"./figures/{event_id}"):
        os.makedirs(f"./figures/{event_id}")

    mgmt = pyatoa.Manager(ds=ds)
    stations = ds.waveforms.list()

    # begin the workflow by looping through all stations
    for station in ds.waveforms.list():
        net, sta = station.split('.')
        try:
            mgmt.load(station, path=f"{model}/{step}")
            mgmt.standardize()
            mgmt.preprocess()
            mgmt.window()
            mgmt.measure()
            mgmt.plot(save="./figures/{event_id}/wav_{sta}",
                      append_title=(f"pyflex={mgmt.cfg.pyflex_preset}; "
                                    f"pyadjoint={mgmt.cfg.adj_src_type}"), 
                      show=True,
                      )
            mgmt.reset()
        except Exception as e:
            traceback.print_exc()
            import ipdb;ipdb.set_trace()
            mgmt.reset()
            continue


