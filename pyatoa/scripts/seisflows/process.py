"""
A plugin to seisflows solver.specfem3d_nz.eval_func() to use in evaluating the
misfit functional within the automated workflow, in the context of
requirements mandated by seisflows
"""
import os
import sys
import json
import pyasdf
import pyatoa
import logging
import argparse
import warnings
import traceback
import numpy as np

from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.asdf.extractions import sum_misfits
from pyatoa.utils.asdf.additions import write_stats_to_asdf
from pyatoa.utils.operations.file_generation import create_stations_adjoint, \
                                     write_adj_src_to_ascii, write_misfit_json
from pyatoa.visuals.convert_images import tile_and_combine


def initialize_parser():
    """
    Seisflows calls pyatoa via subprocess, so we use an argparser to provide
    Pyatoa with information that it requires
    :return:
    """
    parser = argparse.ArgumentParser(description="Inputs from Seisflows")
    parser.add_argument("-e", "--event_id")
    parser.add_argument("-m", "--model_number")
    parser.add_argument("-i", "--step_count")
    parser.add_argument("-w", "--working_dir")
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-c", "--current_dir")
    parser.add_argument("-s", "--suffix")
    parser.add_argument("-l", "--logging", default=False)

    return parser.parse_args()


def get_paths(args):
    """
    Make the necessary output directories for Pyatoa
    :type args: argparse.Parser
    :param args: arguments passed in from seisflows
    :rtype *_paths: dict
    :return *_paths: dictionaries for given pyatoa and seisflows pathing
    """
    # pyatoa specific directories
    fig_dir = os.path.join(args.output_dir, "figures", args.model_number, 
                                                                  args.event_id)
    data_dir = os.path.join(args.output_dir, "data")
    for d in [fig_dir, data_dir]:
        if not os.path.exists(d):
            os.makedirs(d)
   
    paths = {
        "EVENT_FIGURES": fig_dir, 
        "PYATOA_DATA": data_dir,
        "PYATOA_FIGURES": os.path.join(args.output_dir, "figures"),
        "MISFIT_FILE": os.path.join(args.output_dir, "misfits.json"),
        "ADJ_TRACES": os.path.join(args.current_dir, "traces", "adj"), 
        "SYN_TRACES": os.path.join(args.current_dir, "traces", "syn"),
        "EVENT_DATA": os.path.join(args.current_dir, "DATA"),
        "STATIONS": os.path.join(args.current_dir, "DATA", "STATIONS")
             } 

    # get paths from seisflows outputs, add them to existing paths
    sfpaths = os.path.join(args.working_dir, "output", "seisflows_paths.json")
    seisflows_paths = json.load(open(sfpaths, "r"))
    for item in seisflows_paths.keys():
        paths[item] = seisflows_paths[item] 
     
    return paths


def finalize(ds, model, args, paths):
    """
    Before finishing the workflow, create the requisite files for
    Specfem and Seisflows
    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset that now contains misfit windows/adjoint srcs
    :type model: str
    :param model: model number, e.g. "m00"
    :type paths: dict
    :param paths: return of get_paths() containing relevant paths 
    """
    print("finalizing")
    # add statistics to auxiliary_data
    write_stats_to_asdf(ds, model, args.step_count)
    
    # create the .sem ascii files required by specfem
    write_adj_src_to_ascii(ds, model, paths["ADJ_TRACES"])
    
    # create the STATIONS_ADJOINT file required by specfem
    create_stations_adjoint(ds, model, paths["EVENT_DATA"])
    
    # sum and write misfits and statistics information to master text file 
    write_misfit_json(ds, model, args.step_count, paths["MISFIT_FILE"]) 

    # tile and combine .png images into a composite .pdf for easy fetching
    tile_and_combine(ds, model, "s{:0<2}".format(args.step_count),
                     paths["PYATOA_FIGURES"], purge_originals=True,
                     purge_tiles=True
                     )    


def process_data(args):
    """
    Main workflow calling on the core functionality of Pyatoa to process 
    observed and synthetic waveforms and perform misfit quantification
   
    :type args: argparse.Parser 
    :param args: arguments passed in from Seisflows
    """
    print("initiating")
    paths = get_paths(args)   

    # set the Pyatoa Config object for misfit qu
    config = pyatoa.Config(
        event_id=args.event_id,
        model_number=args.model_number,
        min_period=10,
        max_period=30,
        filter_corners=4,
        rotate_to_rtz=False,
        unit_output="DISP",
        pyflex_config="UAF",
        adj_src_type="multitaper_misfit",
        paths_to_synthetics=[paths["SYN_TRACES"]],
        paths_to_waveforms=[],
        paths_to_responses=[]
        )

    # save output by event id
    print("running")
    ds_name = os.path.join(paths["PYATOA_DATA"], 
                           "{}.h5".format(config.event_id)
                           )
    with pyasdf.ASDFDataSet(ds_name) as ds:
        clean_ds(ds, config.model_number)
        config.write_to_asdf(ds)

        # calculate misfit by station, get stations from Specfem STATIONS file
        mgmt = pyatoa.Manager(config=config, ds=ds)
        stations = np.loadtxt(paths["STATIONS"], usecols=[0,1], dtype=str)
        for station in stations:
            sta, net = station
            print("{}.{}".format(net, sta))
            try:
                mgmt.reset()
                mgmt.gather_data(station_code="{net}.{sta}.{loc}.{cha}".format(
                                           net=net, sta=sta, loc="*", cha="HH?")
                                 )
                mgmt.preprocess()
                mgmt.run_pyflex()
                mgmt.run_pyadjoint()
                mgmt.plot_wav(save=os.path.join(
                    paths["EVENT_FIGURES"], "wav_{}".format(sta)), show=False)
                mgmt.plot_map(save=os.path.join(
                    paths["EVENT_FIGURES"], "map_{}".format(sta)), show=False)
            except Exception as e:
                traceback.print_exc()
                print("\n")
                continue

        # generate output files for specfem and seisflows
        finalize(ds, args.model_number, args, paths)

    
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    print("Ignoring warnings due to H5PY deprecation warning")
    
    args = initialize_parser()
    if args.logging:
        logger = logging.getLogger("pyatoa")
        logger.setLevel(logging.DEBUG)
    
    process_data(args)


