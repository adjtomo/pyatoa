"""
An example run script to populate a pyasdf dataset with event, stations,
waveforms, misfit windows and adjoint sources
"""
import os
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
event_ids = ["2016p355601", "2017p2922ss46"]

for event_id in event_ids:
    if not os.path.exists("./figures/{}".format(event_id)):
        os.makedirs("./figures/{}".format(event_id))

    config = pyatoa.Config(event_id=event_id, model_number=model_number,
                           min_period=10, max_period=30, filter_corners=4,
                           rotate_to_rtz=False, unit_output="DISP",
                           pyflex_config="UAF",
                           adj_src_type="multitaper_misfit",
                           paths_to_waveforms=[
                               '/seis/prj/fwi/bchow/seismic', '/geonet/seismic',
                               '/Users/chowbr/Documents/subduction/seismic'],
                           paths_to_synthetics=[
                               '/seis/prj/fwi/bchow/oficial/AGU18',
                               '/seis/prj/fwi/bchow/tomo/adjoint_test/adjoint_master',
                               '/Users/chowbr/Documents/subduction/tomo/adjoint_test/adjoint_master'],
                           paths_to_responses=[
                               '/seis/prj/fwi/bchow/seed/RESPONSE',
                               '/Users/chowbr/Documents/subduction/seed/RESPONSE'],
                           )

    # initiate pyasdf dataset where all data will be saved
    ds = pyasdf.ASDFDataSet(
        "/Users/chowbr/Documents/subduction/tomo/adjoint_test/hdf5/{}.h5".format(config.event_id))
    # clean_ds(ds)
    # config.write_to_asdf(ds)

    # begin the Pyatoa Workflow
    mgmt = pyatoa.Manager(config=config, ds=ds)
    # loop through all stations that were interested in processing
    master_inventory = read_inventory(
        "/Users/chowbr/Documents/subduction/data/STATIONXML/"
        "MASTER/master_inventory_slim.xml")
    for net in master_inventory:
        for sta in net:
            if sta.is_active(time=mgmt.event.preferred_origin().time):
                try:
                    mgmt.gather_data(
                        station_code="{net}.{sta}.{loc}.{cha}".format(
                            net=net.code, sta=sta.code, loc="*", cha="HH?")
                    )
                    if not mgmt.check_full():
                        continue
                    mgmt.preprocess()
                    mgmt.run_pyflex()
                    mgmt.run_pyadjoint()
                    # from pyatoa.visuals.plot_waveforms import window_maker
                    # window_maker(st_obs=mgmt.crate.st_obs,
                    #              st_syn=mgmt.crate.st_syn,
                    #              time_offset=mgmt.crate.time_offset,
                    #              unit_output=config.unit_output, config=config,
                    #              show=True)

                    mgmt.plot_wav(save="./figures/{eid}/wav_{sta}".format(
                        eid=config.event_id, sta=sta.code), show="hold",
                    )
                    mgmt.plot_map(save="./figures/{eid}/map_{sta}".format(
                        eid=config.event_id, sta=sta.code), show=True
                    )
                    mgmt.reset()
                except Exception as e:
                    traceback.print_exc()
                    mgmt.reset()
                    continue

    write_adj_src_to_ascii(ds, model_number=model_number, filepath="./SEM/")
    create_stations_adjoint(ds, model_number=model_number, filepath="./SEM/")


