"""
An example run script to populate a pyasdf dataset with event, stations,
waveforms, misfit windows and adjoint sources
"""
import pyasdf
import pyatoa
import logging
from obspy import read_inventory

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
                           '/seis/prj/fwi/bchow/seismic','/geonet/seismic'],
                       paths_to_responses=[
                           '/seis/prj/fwi/bchow/seed/RESPONSE',
                           '/Users/chowbr/Documents/subduction/seed/RESPONSE'],
                       )

# initiate pyasdf dataset where all data will be saved
ds = pyasdf.ASDFDataSet("/seis/prj/fwi/bchow/data/PYASDF/{}.h5".format(
                                                            config.event_id))

# begin the Pyatoa Workflow
proc = pyatoa.Processor(config=config, ds=ds)
proc.launch()
print(proc)

# loop through all stations that were interested in processing
master_inventory = read_inventory(
    '/seis/prj/fwi/bchow/data/STATIONXML/MASTER/master_inventory.xml')
station_code_template = ""
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
                proc.plot_wav()
                proc.plot_map()
                import ipdb; ipdb.set_trace()
            except Exception as e:
                print(e)
                import ipdb; ipdb.set_trace()


