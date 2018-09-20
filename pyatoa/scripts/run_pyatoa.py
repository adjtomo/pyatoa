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


def pathing():
    """
    simple pathname calling for my own testing purposes
    """
    basecheck = os.getcwd().split('/')[1]
    if basecheck == "Users":
        datapath = "/Users/chowbr/Documents/subduction/data/"
        plotpath = "/Users/chowbr/Documents/subduction/plots/adjtomo/{eventid}"
    else:
        datapath = "/seis/prj/fwi/bchow/data/"
        plotpath = "/seis/prj/fwi/bchow/plots/adjtomo/{eventid}"

    return datapath, plotpath


# initiate logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.DEBUG)

# initiate config
config = pyatoa.Config(event_id="2018p130600", model_number=0, min_period=10,
                       max_period=30, filter_corners=4, rotate_to_rtz=True,
                       unit_output="DISP", pyflex_config="UAF",
                       adj_src_type="multitaper_misfit", raw_sampling_rate=50.,
                       paths_to_waveforms=[
                           '/Users/chowbr/Documents/subduction/seismic',
                           '/seis/prj/fwi/bchow/seismic', '/geonet/seismic'],
                       paths_to_responses=[
                           '/seis/prj/fwi/bchow/seed/RESPONSE',
                           '/Users/chowbr/Documents/subduction/seed/RESPONSE'],
                       )

# initiate pyasdf dataset where all data will be saved
ds = pyasdf.ASDFDataSet(os.path.join(pathing()[0], "PYASDF", "{}.h5".format(
    config.event_id))
                        )
config.write_to_pyasdf(ds)

# begin the Pyatoa Workflow
proc = pyatoa.Processor(config=config, ds=ds)
proc.launch()

# loop through all stations that were interested in processing
master_inventory = read_inventory(
    os.path.join(pathing()[0], "STATIONXML", "MASTER", "master_inventory.xml")
    )
for net in master_inventory:
    for sta in net:
        if sta.is_active(time=proc.event.preferred_origin().time):
            try:
                proc.gather_data(
                    station_code="{net}.{sta}.{loc}.{cha}".format(net=net.code,
                                                                  sta=sta.code,
                                                                  loc="*",
                                                                  cha="HH?")
                )
                proc.preprocess()
                proc.run_pyflex()
                proc.run_pyadjoint()
                proc.plot_wav(save=os.path.join(pathing()[1].format(
                    eventid=config.event_id), "wav_{}".format(sta.code)
                    ), show=False
                )
                proc.plot_map(save=os.path.join(pathing()[1].format(
                    eventid=config.event_id), "map_{}".format(sta.code)
                    ), show=False, show_faults=False
                )
                proc.reset()
            except Exception as e:
                traceback.print_exc()
                import ipdb;ipdb.set_trace()
                proc.reset()
                continue


