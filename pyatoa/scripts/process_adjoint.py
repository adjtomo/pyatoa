"""
An example run script to populate a pyasdf dataset with event, stations,
waveforms, misfit windows and adjoint sources
"""
import os
import pyasdf
import pyatoa
import logging
import traceback
import numpy as np
from obspy import read_inventory

from pyatoa.utils.operations.formatting import write_adj_src_to_ascii
from pyatoa.utils.operations.file_generation import create_stations_adjoint

# initiate logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.DEBUG)

# initiate config
model_number = "m00"
config = pyatoa.Config(event_id="2018p266243", model_number=model_number,
                       min_period=10, max_period=30, filter_corners=4,
                       rotate_to_rtz=False, unit_output="DISP",
                       pyflex_config="UAF",
                       adj_src_type="multitaper_misfit",
                       paths_to_waveforms=[
                           '/seis/prj/fwi/bchow/seismic', '/geonet/seismic',
                           '/Users/chowbr/Documents/subduction/seismic'],
                       paths_to_synthetics=[
                           '/seis/prj/fwi/bchow/tomo/adjoint_test/adjoint_master',
                           '/Users/chowbr/Documents/subduction/tomo/adjoint_test/adjoint_master'],
                       paths_to_responses=[
                           '/seis/prj/fwi/bchow/seed/RESPONSE',
                           '/Users/chowbr/Documents/subduction/seed/RESPONSE'],
                       )

# initiate pyasdf dataset where all data will be saved
ds = pyasdf.ASDFDataSet("./{}_adjoint_test.h5".format(config.event_id))
config.write_to_asdf(ds)

# begin the Pyatoa Workflow
mgmt = pyatoa.Manager(config=config, ds=ds)

# loop through all stations that were interested in processing
master_inventory = read_inventory(
    "/Users/chowbr/Documents/subduction/data/STATIONXML/"
    "MASTER/master_inventory.xml")
# "/seis/prj/fwi/bchow/data/STATIONXML/"
# "MASTER/master_inventory.xml"
# )
for net in master_inventory:
    for sta in net:
        if sta.is_active(time=mgmt.event.preferred_origin().time):
            try:
                mgmt.gather_data(
                    station_code="{net}.{sta}.{loc}.{cha}".format(net=net.code,
                                                                  sta=sta.code,
                                                                  loc="*",
                                                                  cha="HH?")
                )
                mgmt.preprocess()
                mgmt.run_pyflex()
                mgmt.run_pyadjoint()
                mgmt.plot_wav(save="./figures/{eid}/wav_{sta}".format(
                    eid=config.event_id, sta=sta.code), show=False
                )
                mgmt.reset()
            except Exception as e:
                traceback.print_exc()
                import ipdb; ipdb.set_trace()
                mgmt.reset()
                continue

write_adj_src_to_ascii(ds, model_number=model_number)
create_stations_adjoint(ds, model_number=model_number)


