"""
A plugin to seisflows solver.specfem3d_nz.eval_func() to use in evaluating the
misfit functional within the automated workflow, in the context of
requirements mandated by Seisflows.

Seisflows is written in Python 2, so Pyatoa interacts with it via an
argument parser and a user defined configuration file.
"""
import os

import sys
import time
import yaml
import json
import glob
import pyasdf
import pyatoa
import shutil
import logging
import argparse
import warnings
import traceback
import numpy as np

from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.asdf.additions import write_stats_to_asdf
from pyatoa.utils.asdf.extractions import windows_from_ds
from pyatoa.utils.tools.io import create_stations_adjoint, write_misfit_json, \
    write_adj_src_to_ascii, write_misfit_stats, tile_combine_imgs


def initialize_parser(args):
    """
    Seisflows calls pyatoa via subprocess, so we use an argparser to provide
    Pyatoa with information that it requires
    
    :rtpye argparser.args
    :return: arguments passed in from external call
    """
    argparser = argparse.ArgumentParser(description="Inputs from Seisflows")

    # Required arguments for all modes
    argparser.add_argument("-m", "--mode", type=str,
                           help="Allows for single Seisflows entry point,"
                                "available: initialize, process, finalize")
    argparser.add_argument("--working_dir", type=str,
                           help="Working dir: main Seisflows directory")
    # Processing arguments
    argparser.add_argument("--event_id", type=str, default=None,
                           help="Event Identifier")
    argparser.add_argument("--model_number", type=str, default=None,
                           help="Model Number, e.g. 'm00'")
    argparser.add_argument("--current_dir", type=str,
                           help="Event Specfem directory, w/ 'DATA', 'traces'")
    argparser.add_argument("--suffix", type=str, default=None,
                           help="Seisflows suffix")
    # Processing and Finalize arguments
    argparser.add_argument("--step_count", type=str, default=None,
                           help="Step Count, e.g. s00")

    return argparser.parse_args(args)


def assemble_paths(parser, mode=''):
    """
    Make the necessary output directories for Pyatoa.
    All paths should be defined and referenced from this function.

    :type parser: argparse.Parser
    :param parser: arguments passed in from seisflows
    :type mode: str
    :param mode: if process, extra paths are appended
    :rtype *_paths: dict
    :return *_paths: dictionaries for given pyatoa and seisflows pathing
    """
    # Hardcoded directory structure, paths and parameters located in this yaml
    config_fid = os.path.join(parser.working_dir, "parameters.yaml")

    # User defined configuration for Pyatoa, controlling processing params
    # pathing and switches for outputs and logging. This is set by the 
    # parameters.yaml file in Seisflows
    with open(config_fid, "r") as f:
        sfparams = yaml.load(f, Loader=yaml.Loader)
    
    # Pyatoa parameters supplied by User in Seisflows paths
    usrcfg = sfparams["PYATOA"]
    pyatoa_io = sfparams["PATHS"]["PYATOA_IO"]

    # Set Pyatoa paths, this is hardcoded directory structure naming
    paths = {
        "CONFIG_FILE": config_fid,
        "PYATOA_FIGURES": os.path.join(pyatoa_io, "figures"), 
        "PYATOA_DATA": os.path.join(pyatoa_io, "data"), 
        "PYATOA_MISFITS": os.path.join(pyatoa_io, "data", "misfits"),
        "PYATOA_MAPS": os.path.join(pyatoa_io, "figures", "maps"),
        "PYATOA_VTKS": os.path.join(pyatoa_io, "figures", "vtks"), 
        "PYATOA_COMPOSITES": os.path.join(pyatoa_io, "figures", "composites"), 
        "PYATOA_SNAPSHOTS": os.path.join(pyatoa_io, "data", "snapshot"),
        "MISFIT_FILE": os.path.join(pyatoa_io, "misfits.json"),
        }


    # Processing requires path names set by Seisflows
    if mode == "process":
        event_paths = {
            "ADJ_TRACES": os.path.join(parser.current_dir, "traces", "adj"),
            "SYN_TRACES": os.path.join(parser.current_dir, "traces", "syn"),
            "OBS_TRACES": os.path.join(parser.current_dir, "traces", "obs"),
            "EVENT_DATA": os.path.join(parser.current_dir, "DATA"),
            "STATIONS": os.path.join(parser.current_dir, "DATA", "STATIONS"),
            "EVENT_FIGURES": os.path.join(paths["PYATOA_FIGURES"],
                                          parser.model_number, parser.event_id),
            "EVENT_MAPS": os.path.join(paths["PYATOA_MAPS"], parser.event_id)
        }
        # Make process specific directories if necessary
        for key in ["EVENT_MAPS", "EVENT_FIGURES"]:
            if not os.path.exists(event_paths[key]):
                os.makedirs(event_paths[key])

        paths = {**paths, **event_paths}

    return paths, usrcfg


def initialize(parser):
    """
    Create necessary directories for Pyatoa to operate with Seisflows

    :type parser: argparse.Parser
    :param parser: arguments from Seisflows
    """
    paths, _ = assemble_paths(parser)

    for key, item in paths.items():
        # Don't make any filenames
        if "FILE" in key:
            continue
        elif not os.path.exists(item):
            os.makedirs(item)


def _snapshot(paths):    
    """
    Internal function to be called by finalize()

    Create a copy of each HDF5 file in the data/snapshot directory
    
    :type paths: dict
    :parma paths: paths from assemble_paths()
    """
    if not os.path.exists(paths["PYATOA_SNAPSHOTS"]):
        os.makedirs(paths["PYATOA_SNAPSHOTS"])
    srcs = glob.glob(os.path.join(paths["PYATOA_DATA"], "*.h5"))
    for src in srcs:
        shutil.copy(src, os.path.join(paths["PYATOA_SNAPSHOTS"], 
                                      os.path.basename(src))
                    )


def _misfit_maps(paths, parser, usrcfg):
    """
    Internal function to be called by finalize()

    Create misfit maps for each of the HDF5 datasets available

    :type paths: dict
    :parma paths: paths from assemble_paths()
    :type parser: argparse.Parser
    :param parser: arguments passed in from Seisflows
    :type usrcfg: dict
    :param usrcfg: user configuration from .yaml file
    """
    from pyatoa.utils.visuals.mapping import event_misfit_map
    from pyatoa.utils.visuals.combine_imgs import combine_images

    name_template = "{eid}_{m}_{s}_misfit_map.png"
    file_ids = []
    # Loop through each available dataset to create misfit map
    datasets = glob.glob(os.path.join(paths["PYATOA_DATA"], "*.h5"))
    for dataset in datasets:
        with pyasdf.ASDFDataSet(dataset) as ds:
            event_id = os.path.basename(ds.filename).split('.')[0]
             
            # Save figures into event directories
            event_figures = os.path.join(paths["PYATOA_FIGURES"], 
                                         parser.model_number, event_id
                                         )
            # Save the fid based on event id, model number, step count
            fidout = os.path.join(event_figures, 
                                  name_template.format(eid=event_id, 
                                                       m=parser.model_number, 
                                                       s=parser.step_count)
                                  )
            file_ids.append(fidout)

            # Use the average misfit to normalize the misfit map
            stats = ds.auxiliary_data.Statistics[
                parser.model_number][parser.step_count].parameters
            average_misfit = stats['average_misfit']
            
            event_misfit_map(map_corners=usrcfg.map_corners,
                             ds=ds, model=parser.model_number,
                             step=parser.step_count,
                             normalize=average_misfit,
                             annotate_station_info='simple',
                             contour_overlay=True, filled_contours=True,
                             show=False, save=fidout
                             )

    # Combine all the misfit maps into a single pdf
    save_to = os.path.join(
                    paths["PYATOA_COMPOSITES"], 
                    f"{parser.model_number}_{parser.step_count}_misfitmaps.pdf"
    )
    combine_images(file_ids=file_ids, save_to=save_to, purge=False)


def finalize(parser):
    """
    Before finishing the iteration, create some final objects

    :type parser: argparse.Parser
    :param parser: arguments passed in from Seisflows
    """
    paths, usrcfg = assemble_paths(parser)

    # Plot the output.optim file outputted by Seisflows
    from pyatoa.utils.visuals.statistics import plot_output_optim
    plot_output_optim(
        path_to_optim=os.path.join(parser.working_dir, "output.optim"),
        save=os.path.join(paths["PYATOA_FIGURES"], "output_optim.png")
    )

    # Generate .vtk files for given source and receivers
    if usrcfg["create_srcrcv_vtk"]:
        from pyatoa.utils.tools.io import create_srcrcv_vtk_multiple
        create_srcrcv_vtk_multiple(
            pathin=paths["PYATOA_DATA"], pathout=paths["PYATOA_VTKS"],
            model=parser.model_number
        )

    # Create copies of .h5 files at the end of each iteration, because .h5
    # files are easy to corrupt so it's good to have a backup
    if usrcfg["snapshot"]:
        _snapshot(paths)

    # Create misfit maps for each event with contour overlay showing misfit
    # Only create misfit maps for the first step count
    if usrcfg["plot_misfit_maps"] and (parser.step_count == "s00"):
        _misfit_maps(paths, parser, usrcfg)


def process(parser):
    """
    Main workflow calling on the core functionality of Pyatoa to process
    observed and synthetic waveforms and perform misfit quantification

    :type parser: argparse.Parser
    :param parser: arguments passed in from Seisflows
    """
    paths, usrcfg = assemble_paths(parser, mode="process")

    # Set logging output for Pyflex and Pyatoa, less output using 'info'
    if usrcfg["set_logging"]:
        logger_pyatoa = logging.getLogger("pyatoa")
        logger_pyflex = logging.getLogger("pyflex")
        if usrcfg["set_logging"].lower() == "info":
            logger_pyatoa.setLevel(logging.INFO)
            logger_pyflex.setLevel(logging.INFO)
        else:
            logger_pyatoa.setLevel(logging.DEBUG)
            logger_pyflex.setLevel(logging.DEBUG)

    # Read in the Pyatoa Config object and set some attributes based on workflow
    config = pyatoa.Config(yaml_fid=paths["CONFIG_FILE"])
    config.event_id = parser.event_id
    config.model_number = parser.model_number
    config.synthetic_tag = f"synthetic_{parser.model_number}"

    # Make sure Pyatoa knows to look in the Seisflows directories for data
    config.cfgpaths["synthetics"].append(paths["SYN_TRACES"])
    config.cfgpaths["waveforms"].append(paths["OBS_TRACES"])

    # Save HDF5 output by event id
    ds_name = os.path.join(paths["PYATOA_DATA"], f"{config.event_id}.h5")
    with pyasdf.ASDFDataSet(ds_name) as ds:
        # Make sure the ASDFDataSet doesn't already contain auxiliary_data
        # because it will be collected in this workflow
        clean_ds(ds=ds, model=config.model_number, step=parser.step_count,
                 fix_windows=usrcfg["fix_windows"])

        # Write the Config to auxiliary_data for provenance
        config.write(write_to=ds)

        # Instantiate the Manager
        mgmt = pyatoa.Manager(config=config, ds=ds)

        # Get stations from Specfem STATIONS file in form NET STA LAT LON ...
        stations = np.loadtxt(paths["STATIONS"], usecols=[0, 1, 2, 3],
                              dtype=str)
        coords = stations[:, 2:]

        # Loop through stations and invoke Pyatoa workflow
        for station in stations:
            sta, net = station[:2]
            print(f"{net}.{sta}")
            try:
                mgmt.reset()

                # Gather data, searching internal pathways, else fetching from
                # external pathways if possible. Preprocess identically
                mgmt.gather(station_code=f"{net}.{sta}.*.HH*")
                mgmt.standardize()
                mgmt.preprocess()

                # Either no fixed misfit windows or no windows exist yet
                if not usrcfg["fix_windows"] or \
                        not hasattr(ds.auxiliary_data.MisfitWindows,
                                    config.model_number):
                    mgmt.window()
                else:
                    # If windows exist and fixed windows, grab from ASDF dataset
                    misfit_windows = windows_from_ds(
                                              ds, config.model_number, net, sta)
                    mgmt.windows = misfit_windows

                mgmt.measure()

                # Plot waveforms with misfit windows and adjoint sources
                if usrcfg["plot_waveforms"]:
                    # Format some strings to append to the waveform plot title
                    append_title = (
                        f"\n{config.model_number}{parser.step_count} "
                        f"pyflex={config.pyflex_map}, "
                        f"pyadjoint={config.adj_src_type},")
                    if mgmt.misfit is not None:
                        append_title = " ".join([append_title, 
                                                 f"misfit={mgmt.misfit:.2E}"])
                    f = mgmt.plot(append_title=append_title, save=os.path.join(
                                  paths["EVENT_FIGURES"], f"wav_{sta}"),
                                  show=False, return_figure=True
                                  )

                # Plot source-receiver maps, don't make a map if no wav data
                # Don't make the map if the map has already been made
                if usrcfg["plot_srcrcv_maps"] and f:
                    map_fid = os.path.join(paths["EVENT_MAPS"], f"map_{sta}")
                    if not os.path.exists(map_fid):
                        mgmt.srcrcvmap(stations=coords, save=map_fid, 
                                       show=False)

                print("\n")
            # Traceback ensures more detailed error tracking
            except Exception:
                traceback.print_exc()
                print("\n")
                continue

        print("writing stats to ASDF file...")
        write_stats_to_asdf(ds, config.model_number, parser.step_count)

        print("writing adjoint sources to .sem? files...")
        write_adj_src_to_ascii(ds, config.model_number, paths["ADJ_TRACES"])

        print("creating STATIONS_ADJOINT file...")
        create_stations_adjoint(ds, config.model_number,
                                specfem_station_file=paths["STATIONS"],
                                pathout=paths["EVENT_DATA"])

        print("writing individual misfit to file...")
        write_misfit_stats(ds, config.model_number, paths["PYATOA_MISFITS"])

        print("writing misfits.json file...")
        write_misfit_json(ds, parser.model_number, parser.step_count,
                          paths["MISFIT_FILE"])

        # Only run this for the first 'step', otherwise we get too many pdfs
        if usrcfg["combine_imgs"] and (parser.step_count == "s00"):
            print("creating composite pdf...")

            # Create the name of the pdf to save to
            save_to = os.path.join(
                    paths["PYATOA_COMPOSITES"], 
                    f"{config.event_id}_{config.model_number}_"
                    f"{parser.step_count}_wavmap.pdf"
                    )
            tile_combine_imgs(ds=ds, save_pdf_to=save_to,
                              wavs_path=paths["EVENT_FIGURES"],
                              maps_path=paths["EVENT_MAPS"],
                              purge_wavs=usrcfg["purge_waveforms"],
                              purge_tiles=usrcfg["purge_tiles"]
                              )


if __name__ == "__main__":
    # TO DO: change this, ignoring all warnings is pretty dangerous
    # Ignoring warnings due to H5PY deprecation warning
    warnings.filterwarnings("ignore")

    # Arguments passed in by Seisflows
    parser = initialize_parser(sys.argv[1:])

    # Initialize Pyatoa directory structure
    if parser.mode == "initialize":
        try:
            initialize(parser)
            sys.exit(0)
        except Exception as e:
            traceback.print_exc()
            sys.exit(1)
    # Run some cleanup scripts at the end of an iteration
    elif parser.mode == "finalize":
        try:
            finalize(parser)
            sys.exit(0)
        except Exception as e:
            traceback.print_exc()
            sys.exit(1)
    # Process misfit values
    elif parser.mode == "process":
        # Run Pyatoa, return successful exit code. Time for posterity
        try:
            tstart = time.time()
            process(parser)
            print(f"{(time.time() - tstart) / 60.:.2f}m elapsed")
            sys.exit(0)
        except Exception as e:
            traceback.print_exc()
            sys.exit(1)
    else:
        print("invalid 'mode' argument")
        sys.exit(1)
