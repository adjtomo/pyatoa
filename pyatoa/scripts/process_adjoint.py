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


def write_adj_src(ds, model_number):
    """
    write adj src to text file for use in specfem
    :param ds:
    :param model_number:
    :return:
    """
    shortcut = ds.auxiliary_data.AdjointSources[model_number]
    for adj_src in shortcut.list():
        template = "./{}.adj".format(shortcut[adj_src].path.replace('_', '.'))
        np.savetxt(template, shortcut[adj_src].data.value)


def write_all_adj_src(ds, model_number, comp_list=["N","E","Z"]):
    """
    write adj src to text file for use in specfem
    :param ds:
    :param model_number:
    :return:
    """
    shortcut = ds.auxiliary_data.AdjointSources[model_number]
    for adj_src in shortcut.list():
        station = shortcut[adj_src].path.replace('_', '.')
        template = "./{}.adj".format(station)
        np.savetxt(template, shortcut[adj_src].data.value)
        for comp in comp_list:
            station_check = station[:-1] + comp
            if station_check.replace('.', '_') not in shortcut.list():
                blank_adj_src = shortcut[adj_src].data.value
                blank_adj_src[:, 1] = np.zeros(len(blank_adj_src[:, 1]))
                blank_template = "./{}.adj".format(station_check)
                np.savetxt(blank_template, blank_adj_src)



def create_stations_adjoint(ds, model_number):
    """
    generate an adjoint stations file for specfem iput
    :param ds:
    :param model_number:
    :return:
    """
    master_station_list = ('/Users/chowbr/Documents/subduction/data/'
                           'STATIONXML/MASTER/master_station_list.txt'
                           )
    stas_with_adjsrcs = []
    for code in ds.auxiliary_data.AdjointSources[model_number].list():
        stas_with_adjsrcs.append(code.split('_')[1])
    stas_with_adjsrcs = set(stas_with_adjsrcs)

    with open(master_station_list, "r") as f:
        lines = f.readlines()

    with open("STATIONS_ADJOINT", "w") as f:
        for line in lines:
            if line.split()[0] in stas_with_adjsrcs:
                    f.write(line)


# initiate logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.DEBUG)

# initiate config
config = pyatoa.Config(event_id="2018p130600", model_number=0, min_period=10,
                       max_period=30, filter_corners=4, rotate_to_rtz=False,
                       unit_output="DISP", pyflex_config="UAF",
                       adj_src_type="multitaper_misfit",
                       paths_to_waveforms=[
                           '/Users/chowbr/Documents/subduction/seismic'],
                       paths_to_synthetics=[
                           '/Users/chowbr/Documents/subduction/tomo/adjoint_test/adjoint_master'],
                       paths_to_responses=[
                           '/Users/chowbr/Documents/subduction/seed/RESPONSE']
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
for net in master_inventory:
    for sta in net:
        if sta.is_active(time=mgmt.event.preferred_origin().time):
            try:
                # if os.path.exists('./figures/wav_{}'.format(sta.code)):
                #     continue
                mgmt.gather_data(
                    station_code="{net}.{sta}.{loc}.{cha}".format(net=net.code,
                                                                  sta=sta.code,
                                                                  loc="*",
                                                                  cha="HH?")
                )
                mgmt.preprocess()
                mgmt.run_pyflex()
                mgmt.run_pyadjoint()
                mgmt.plot_wav(save="./figures/wav_{}".format(sta.code),
                              show=False
                              )
                mgmt.plot_map(save="./figures/map_{}".format(sta.code),
                              show=False, show_faults=False
                              )
                mgmt.reset()
            except Exception as e:
                traceback.print_exc()
                mgmt.reset()
                continue

# write_adj_src(ds, model_number="m00")


