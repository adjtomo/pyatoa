"""
A plugin to seisflows solver.specfem3d_nz.eval_func() to use in evaluating the
misfit functional within the automated workflow, in the context of
requirements mandated by seisflows
"""
import os
import sys
import glob
import time
import json
import pyasdf
import pyatoa
import logging
import argparse
import warnings
import traceback
import numpy as np
from obspy import read_inventory

from pyatoa.utils.operations.pyasdf_editing import clean_ds, sum_misfits
from pyatoa.utils.operations.formatting import write_adj_src_to_ascii
from pyatoa.utils.operations.file_generation import create_stations_adjoint


def initialize_parser():
    """
    Seisflows calls pyatoa via subprocess, so we use an argparser to provide
    Pyatoa with information that it requires
    :return:
    """
    parser = argparse.ArgumentParser(description="Inputs for Seisflows")
    parser.add_argument("-i", "--event_id")
    parser.add_argument("-m", "--model_number")
    parser.add_argument("-w", "--working_dir")
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-c", "--current_dir")
    parser.add_argument("-s", "--suffix")
    parser.add_argument("-l", "--logging", default=True)

    return parser.parse_args()


def set_logging(set_bool=True):
    """
    Turn logging on or off for Pyatoa, lots of information is spit out so for
    fully automatic workflows, best to turn off
    :param set_bool:
    :return:
    """
    if set_bool:
        logger = logging.getLogger("pyatoa")
        logger.setLevel(logging.DEBUG)


def get_paths(working_dir):
    """
    Take advantage of the fact that seisflows stores all its path information
    in a json file. Use that json file to set all the paths in Pyatoa
    Not used
    """
    seisflows_paths = os.path.join(
                                  working_dir, "output", "seisflows_paths.json")
    with open(seisflows_paths, "r") as f:
        path_dict = json.load(f)

    return path_dict  
 

def make_directories(args):
    """
    Make the necessary output directories for Pyatoa
    """
    # pyatoa specific directories
    fig_dir = os.path.join(args.output_dir, "figures", args.model_number, 
                                                                  args.event_id)
    data_dir = os.path.join(args.output_dir, "data")
    misfit_txt = os.path.join(args.output_dir, "pyatoa.misfits")
    for d in [fig_dir, data_dir]:
        if not os.path.exists(d):
            os.makedirs(d)
   
    pyatoa_paths = {
        "FIGURES": fig_dir, 
        "PYATOA_DATA": data_dir,
        "MISFIT_FILE": misfit_txt,
        "ADJ_TRACES": os.path.join(args.current_dir, "traces", "adj"), 
        "SYN_TRACES": os.path.join(args.current_dir, "traces", "syn"),
        "EVENT_DATA": os.path.join(args.current_dir, "DATA")
                   }
    return pyatoa_paths


def write_pyatoa_misfit(misfit, misfit_file):
    """
    Seisflows specific function
    write misfit into a text file for each event for easy checks
    """
    with open(misfit_file, "a") as f:
        f.write("{}\n".format(misfit))


def process_data(args):
    """
    Main workflow for Pyatoa to process data
    :param args:
    :return:
    """
    pyatoa_paths = make_directories(args)   

    # set the pyatoa config object for misfit quantification
    config = pyatoa.Config(
        event_id=args.event_id, model_number=args.model_number, 
        min_period=10, max_period=30, filter_corners=4, 
        rotate_to_rtz=False, unit_output="DISP", pyflex_config="UAF",
        adj_src_type="multitaper_misfit",
        paths_to_synthetics= [pyatoa_paths["SYN_TRACES"]],
        paths_to_waveforms=[], paths_to_responses=[]
        )

    # save output by event id
    ds_name = os.path.join(pyatoa_paths["PYATOA_DATA"], 
                           "{}.h5".format(config.event_id)
                          ) 
    with pyasdf.ASDFDataSet(ds_name) as ds:
        clean_ds(ds)
        config.write_to_asdf(ds)

        # calculate misfit by station
        mgmt = pyatoa.Manager(config=config, ds=ds)
        master_inventory = read_inventory(os.path.join(
                                  args.working_dir, "master_inventory_slim.xml")
        )
        for net in master_inventory:
            for sta in net:
                if sta.is_active(time=mgmt.event.preferred_origin().time):
                    try:
                        mgmt.gather_data(
                            station_code="{net}.{sta}.{loc}.{cha}".format(
                                net=net.code, sta=sta.code, loc="*", cha="HH?")
                        )
                        mgmt.preprocess()
                        mgmt.run_pyflex()
                        mgmt.run_pyadjoint()
                        mgmt.plot_wav(save=os.path.join(pyatoa_paths["FIGURES"], 
                                   "wav_{sta}".format(sta=sta.code)), show=False
                        )
                        mgmt.plot_map(save=os.path.join(pyatoa_paths["FIGURES"], 
                                   "map_{sta}".format(sta=sta.code)), show=False
                        )
                        mgmt.reset()
                    except Exception as e:
                        print("\n")
                        traceback.print_exc()
                        mgmt.reset()
                        continue

        # generate output files for specfem and seisflows
        write_adj_src_to_ascii(
                ds, args.model_number, pyatoa_paths["ADJ_TRACES"])
        create_stations_adjoint(
                ds, args.model_number, pyatoa_paths["EVENT_DATA"])
        source_misfit = sum_misfits(ds, args.model_number)
        write_pyatoa_misfit(source_misfit, pyatoa_paths["MISFIT_FILE"])
    
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    print("Ignoring warnings due to H5PY deprecation warning")
    args = initialize_parser()
    set_logging(set_bool=args.logging)
    process_data(args)
    time.sleep(15)    


