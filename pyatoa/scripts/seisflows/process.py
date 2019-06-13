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
from pyatoa.utils.asdf.additions import write_stats_to_asdf
from pyatoa.utils.operations.file_generation import create_stations_adjoint, \
                   write_adj_src_to_ascii, write_misfit_json, write_misfit_stats
from pyatoa.visuals.convert_images import tile_and_combine
from pyatoa.scripts.seisflows.sfconfig import sfconfig


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

    return parser.parse_args()


def get_paths(args, usrcfg):
    """
    Make the necessary output directories for Pyatoa
    :type args: argparse.Parser
    :param args: arguments passed in from seisflows
    :rtype *_paths: dict
    :return *_paths: dictionaries for given pyatoa and seisflows pathing
    """
    # pyatoa specific directories
    fig_dir = os.path.join(args.output_dir, usrcfg["figure_dir"],
                           args.model_number, args.event_id
                           )
    data_dir = os.path.join(args.output_dir, usrcfg["data_dir"])
    misfit_dir = os.path.join(args.output_dir, usrcfg["misfit_dir"])
    for d in [fig_dir, data_dir, misfit_dir]:
        if not os.path.exists(d):
            os.makedirs(d)
   
    paths = {
        "EVENT_FIGURES": fig_dir, 
        "PYATOA_DATA": data_dir,
        "PYATOA_MISFITS": misfit_dir,
        "PYATOA_FIGURES": os.path.join(args.output_dir, usrcfg["figure_dir"]),
        "MISFIT_FILE": os.path.join(args.output_dir, usrcfg["misfits_json"]),
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


def finalize(ds, model, args, paths, usrcfg):
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
    # REQUIRED: add statistics to auxiliary_data
    write_stats_to_asdf(ds, model, args.step_count)
    
    # REQUIRED: create the .sem ascii files required by specfem
    write_adj_src_to_ascii(ds, model, paths["ADJ_TRACES"])
    
    # REQUIRED: create the STATIONS_ADJOINT file required by specfem
    create_stations_adjoint(ds, model, paths["EVENT_DATA"])
    
    # REQUIRED: write misfits for seisflows into individual text files
    write_misfit_stats(ds, model, paths["PYATOA_MISFITS"])

    # OPTIONAL: sum and write misfits information to a json file
    if usrcfg["write_mifsits_json"]:
        write_misfit_json(ds, model, args.step_count, paths["MISFIT_FILE"])

    # OPTIONAL: combine .png images into a composite .pdf for easy fetching
    if usrcfg["tile_and_combine"]:
        tile_and_combine(ds, model, "s{:0<2}".format(args.step_count),
                         paths["PYATOA_FIGURES"],
                         purge_originals=usrcfg["purge_originals"],
                         purge_tiles=usrcfg["purge_tiles"]
                         )


def process_data(args, usrcfg):
    """
    Main workflow calling on the core functionality of Pyatoa to process 
    observed and synthetic waveforms and perform misfit quantification
   
    :type args: argparse.Parser 
    :param args: arguments passed in from Seisflows
    :type usrcfg: dict
    :param usrcfg: user configuration, from sfconfig.py
    """
    print("initiating")
    paths = get_paths(args, usrcfg)

    # Set the Pyatoa Config object for misfit quantification
    config = pyatoa.Config(
        event_id=args.event_id,
        model_number=args.model_number,
        min_period=usrcfg["min_period"],
        max_period=usrcfg["max_period"],
        filter_corners=usrcfg["filter_corners"],
        rotate_to_rtz=usrcfg["rotate_ro_rtz"],
        unit_output=usrcfg["unit_output"],
        pyflex_config=usrcfg["pyflex_config"],
        adj_src_type=usrcfg["adj_src_type"],
        paths_to_synthetics=[paths["SYN_TRACES"]],
        paths_to_waveforms=[usrcfg["paths_to_waveforms"]],
        paths_to_responses=[usrcfg["paths_to_responses"]]
        )

    # Save HDF5 output by event id
    print("running")
    ds_name = os.path.join(paths["PYATOA_DATA"], 
                           "{}.h5".format(config.event_id)
                           )
    with pyasdf.ASDFDataSet(ds_name) as ds:
        # Make sure the ASDFDataSet doesn't already contain auxiliary_data
        # that will be collected here.
        clean_ds(ds, config.model_number, usrcfg["fix_windows"])

        # Write the Config to auxiliary_data for provenance
        config.write_to_asdf(ds)

        # Calculate misfit by station, get stations from Specfem STATIONS file
        mgmt = pyatoa.Manager(config=config, ds=ds)
        stations = np.loadtxt(paths["STATIONS"], usecols=[0, 1], dtype=str)

        # Loop through stations and invoke Pyatoa workflow
        for station in stations:
            sta, net = station
            print("{}.{}".format(net, sta))
            try:
                mgmt.reset()
                mgmt.gather_data(station_code="{net}.{sta}.{loc}.{cha}".format(
                                           net=net, sta=sta, loc="*", cha="HH?")
                                 )
                mgmt.preprocess()

                # Misfit Window Fixing:
                # if 'fix_windows==False', OR 'fix_windows==True' BUT
                # misfit windows do not exist yet, run Pyflex to create windows
                if not usrcfg["fix_windows"] or not \
                        hasattr(ds.auxiliary_data.MisfitWindows,
                                config.model_number):
                    mgmt.run_pyflex()

                mgmt.run_pyadjoint()
                mgmt.plot_wav(save=os.path.join(
                    paths["EVENT_FIGURES"], "wav_{}".format(sta)), show=False)
                mgmt.plot_map(save=os.path.join(
                    paths["EVENT_FIGURES"], "map_{}".format(sta)), show=False)
            # Traceback ensures more detailed error tracking, but do not allow
            # error to stop the workflow
            except Exception:
                traceback.print_exc()
                print("\n")
                continue

        # Generate output files for specfem and seisflows
        finalize(ds, args.model_number, args, paths, usrcfg)

    
if __name__ == "__main__":
    try:
        # Ignoring warnings due to H5PY deprecation warning
        warnings.filterwarnings("ignore")

        # Arguments passed in by Seisflows, related to model number, step count
        # and relevant pathing
        args = initialize_parser()

        # User defined configuration for Pyatoa, controlling processing parameters
        # pathing and switches for outputs and logging
        usrcfg = sfconfig()

        # Set logging based on user defined parameter
        if usrcfg["logging"]:
            logger = logging.getLogger("pyatoa")
            logger.setLevel(logging.DEBUG)

        # Run Pyatoa, return successful exit code
        process_data(args, usrcfg)
        sys.exit(0)
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)


