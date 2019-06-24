"""
A plugin to seisflows solver.specfem3d_nz.eval_func() to use in evaluating the
misfit functional within the automated workflow, in the context of
requirements mandated by seisflows. Seisflows interacts with Pyatoa via an
argument parser and a user defined configuration dictionary.
"""
import os

import sys
import json
import glob
import pyasdf
import pyatoa
import logging
import argparse
import warnings
import traceback
import numpy as np

from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.asdf.additions import write_stats_to_asdf
from pyatoa.utils.asdf.extractions import windows_from_ds
from pyatoa.utils.operations.file_generation import create_stations_adjoint, \
                write_misfit_json, write_adj_src_to_ascii, write_misfit_stats


def initialize_parser():
    """
    Seisflows calls pyatoa via subprocess, so we use an argparser to provide
    Pyatoa with information that it requires
        
    :return:
    """
    parser = argparse.ArgumentParser(description="Inputs from Seisflows")

    # Required arguments for all modes
    parser.add_argument("-m", "--mode", type=str,
                        help="Allows for single Seisflows entry point,"
                             "available: initialize, process, finalize")
    parser.add_argument("--working_dir", type=str,
                        help="Working dir: main Seisflows directory")

    # Processing arguments
    parser.add_argument("--event_id", type=str, default=None,
                        help="Event Identifier")
    parser.add_argument("--model_number", type=str, default=None,
                        help="Model Number, e.g. 'm00'")
    parser.add_argument("--step_count", type=int, default=None,
                        help="Step Count, e.g. 0")
    parser.add_argument("--current_dir", type=str,
                        help="Event Specfem directory, with 'DATA', 'traces'")
    parser.add_argument("--suffix", type=str, default=None,
                        help="Seisflows suffix")

    return parser.parse_args()


def initialize(args):
    """
    Create necessary directories for Pyatoa to operate with Seisflows

    :type args: arparse.Parser
    :param args: arguments from Seisflows
    :type usrcfg: dict
    :param usrcfg: user configuration, from sfconfig.py
    """
    pt_paths, _ = pyatoa_paths(args)

    for key in pt_paths.keys():
        # Don't make any filenames
        if "FILE" in key:
            continue
        elif not os.path.exists(pt_paths[key]):
            os.makedirs(pt_paths[key])


def finalize(args):
    """
    Before finishing the iteration, create some final objects

    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset that now contains misfit windows/adjoint srcs
    :type model: str
    :param model: model number, e.g. "m00"
    :type paths: dict
    :param paths: return of get_paths() containing relevant paths
    """
    paths, usrcfg = assemble_paths(args)

    # Loop through all available datasets
    for fid in glob.glob(os.path.join(paths["PYATOA_DATA"], "*.h5")):
        with pyasdf.ASDFDataSet(fid) as ds:
            # Generate .vtk files for given source and receivers
            # TO DO: only create one .vtk file for all datasets
            if usrcfg["create_srcrcv_vtk"]:
                from pyatoa.utils.operations.file_generation import \
                    create_srcrcv_vtk
                create_srcrcv_vtk(ds=ds, model=model,
                                  path_out=paths["PYATOA_VTKS"],
                                  event_separate=usrcfg["create_src_vtk"],
                                  )


def pyatoa_paths(args):
    """
    Pyatoa requires pre-built directories to save data and figures. This only
    needs to be run once at the beginning of the inversion workflow.

    :type args: arparse.Parser
    :param args: arguments from Seisflows
    :type usrcfg: dict
    :param usrcfg: user configuration, from sfconfig.py
    """
    # Get paths from Seisflows outputs .json file
    sf_json = os.path.join(args.working_dir, "output", "seisflows_paths.json")
    sf_paths = json.load(open(sf_json, "r"))

    # User defined configuration for Pyatoa, controlling processing params
    # pathing and switches for outputs and logging. Import from user supplied
    # output directory
    os.chdir(sf_paths["PYATOA"])
    from config import sfconfig
    usrcfg = sfconfig()

    # Set Pyatoa paths based on user configurations
    figs = os.path.join(sf_paths["PYATOA"], usrcfg["figure_dir"])
    vtks = os.path.join(figs, "vtks")

    data = os.path.join(sf_paths["PYATOA"], usrcfg["data_dir"])
    misfits = os.path.join(data, "misfits")

    misfit_file = os.path.join(sf_paths["PYATOA"], usrcfg["misfits_json"])

    pt_paths = {"PYATOA_FIGURES": figs, "PYATOA_DATA": data,
                "PYATOA_VTKS": vtks, "PYATOA_MISFITS": misfits,
                "MISFIT_FILE": misfit_file
                }

    return pt_paths, usrcfg


def assemble_paths(args, usrcfg):
    """
    Make the necessary output directories for Pyatoa
    :type args: argparse.Parser
    :param args: arguments passed in from seisflows
    :rtype *_paths: dict
    :return *_paths: dictionaries for given pyatoa and seisflows pathing
    """
    pt_paths, usrcfg = pyatoa_paths(args, usrcfg)

    event_paths = {
        "ADJ_TRACES": os.path.join(args.current_dir, "traces", "adj"),
        "SYN_TRACES": os.path.join(args.current_dir, "traces", "syn"),
        "EVENT_DATA": os.path.join(args.current_dir, "DATA"),
        "STATIONS": os.path.join(args.current_dir, "DATA", "STATIONS"),
        "EVENT_FIGURES": os.path.join(pt_paths["PYATOA_FIGURES"],
                                      args.model_number, args.event_id
                                      )
    }

    paths = {**pt_paths, **event_paths}

    return paths, usrcfg


def process(args, usrcfg):
    """
    Main workflow calling on the core functionality of Pyatoa to process 
    observed and synthetic waveforms and perform misfit quantification
   
    :type args: argparse.Parser 
    :param args: arguments passed in from Seisflows
    :type usrcfg: dict
    :param usrcfg: user configuration, from sfconfig.py
    """
    paths, usrcfg = assemble_paths(args, usrcfg)

    # Set loggging output
    if usrcfg["set_logging"]:
        logger = logging.getLogger("pyatoa")
        logger.setLevel(logging.DEBUG)

    # Set the Pyatoa Config object for misfit quantification
    config = pyatoa.Config(
        event_id=args.event_id,
        model_number=args.model_number,
        min_period=usrcfg["min_period"],
        max_period=usrcfg["max_period"],
        filter_corners=usrcfg["filter_corners"],
        rotate_to_rtz=usrcfg["rotate_to_rtz"],
        unit_output=usrcfg["unit_output"],
        pyflex_config=usrcfg["pyflex_config"],
        adj_src_type=usrcfg["adj_src_type"],
        paths_to_synthetics=[paths["SYN_TRACES"]],
        paths_to_waveforms=usrcfg["paths_to_waveforms"],
        paths_to_responses=usrcfg["paths_to_responses"]
        )

    # The trial step number is only used sparsely, e.g. the Statistics subgroup
    step_number = "s{:0>2}".format(args.step_count)

    # Save HDF5 output by event id
    ds_name = os.path.join(paths["PYATOA_DATA"], 
                           "{}.h5".format(config.event_id)
                           )
    with pyasdf.ASDFDataSet(ds_name) as ds:
        # Make sure the ASDFDataSet doesn't already contain auxiliary_data
        # because it will be collected in this workflow
        clean_ds(ds=ds, model=config.model_number, step=step_number,
                 fix_windows=usrcfg["fix_windows"])

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
                
                # Gather data, searching internal pathways, else fetching from
                # external pathways if possible. Preprocess identically
                mgmt.gather_data(station_code="{net}.{sta}.{loc}.{cha}".format(
                                 net=net, sta=sta, loc="*", cha="HH[NZE]")
                                 )
                mgmt.preprocess()
                
                # Either no fixed misfit windows or no windows exist yet
                if not usrcfg["fix_windows"] or \
                        not hasattr(ds.auxiliary_data.MisfitWindows,
                                    config.model_number):
                    mgmt.run_pyflex()
                else:
                    # If windows exist and fixed windows, grab from ASDF dataset
                    misfit_windows = windows_from_ds(
                                              ds, config.model_number, net, sta)
                    mgmt.overwrite_windows(misfit_windows)

                mgmt.run_pyadjoint()

                # Plot waveforms with misfit windows and adjoint sources
                if usrcfg["plot_waveforms"]:
                    mgmt.plot_wav(
                        save=os.path.join(
                                  paths["EVENT_FIGURES"], "wav_{}".format(sta)), 
                        show=False
                        )
               
                # Plot source-receiver maps  
                if usrcfg["plot_maps"]:
                    mgmt.plot_map(
                        save=os.path.join(
                                  paths["EVENT_FIGURES"], "map_{}".format(sta)),
                        show=False
                        )
                print("\n")
            # Traceback ensures more detailed error tracking
            except Exception:
                traceback.print_exc()
                print("\n")
                continue

        # Add statistics to auxiliary_data
        write_stats_to_asdf(ds, config.model_number, args.step_count)

        # Create the .sem ascii files required by specfem
        write_adj_src_to_ascii(ds, config.model_number, paths["ADJ_TRACES"])

        # Create the STATIONS_ADJOINT file required by specfem
        create_stations_adjoint(ds, config.model_number, paths["EVENT_DATA"])

        # Write misfits for seisflows into individual text files
        write_misfit_stats(ds, config.model_number, paths["PYATOA_MISFITS"])

        # Sum and write misfits information to a JSON file
        write_misfit_json(ds, args.model_number, args.step_count,
                          paths["MISFIT_FILE"])

        # Combine .png images into a composite .pdf for easy fetching
        if usrcfg["tile_and_combine"]:
            from pyatoa.utils.visuals.convert_images import tile_and_combine
            tile_and_combine(ds=ds, model=args.model_number, step=step_number,
                             figure_path=paths["PYATOA_FIGURES"],
                             purge_originals=usrcfg["purge_originals"],
                             purge_tiles=usrcfg["purge_tiles"]
                             )

    
if __name__ == "__main__":
    # Ignoring warnings due to H5PY deprecation warning
    warnings.filterwarnings("ignore")

    # Arguments passed in by Seisflows
    args = initialize_parser()

    # Initialize Pyatoa directory structure
    if args.mode == "initialize":
        initialize(args)

    # Run some cleanup scripts at the end of an iteration
    elif args.mode == "finalize":
        finalize(args)

    # Process misfit values
    elif args.mode == "process":
        # Run Pyatoa, return successful exit code
        try:
            process(args)
            sys.exit(0)
        except Exception as e:
            traceback.print_exc()
            sys.exit(1)

    else:
        print("invalid 'mode' argument")
        sys.exit(1)
