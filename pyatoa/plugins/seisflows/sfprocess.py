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
import shutil
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


def initialize_parser(args):
    """
    Seisflows calls pyatoa via subprocess, so we use an argparser to provide
    Pyatoa with information that it requires

    :return:
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
    argparser.add_argument("--step_count", type=int, default=None,
                           help="Step Count, e.g. 0")
    argparser.add_argument("--current_dir", type=str,
                           help="Event Specfem directory, w/ 'DATA', 'traces'")
    argparser.add_argument("--suffix", type=str, default=None,
                           help="Seisflows suffix")

    return argparser.parse_args(args)


def assemble_paths(parser, mode=''):
    """
    Make the necessary output directories for Pyatoa
    :type parser: argparse.Parser
    :param parser: arguments passed in from seisflows
    :type mode: str
    :param mode: if process, extra paths are appended
    :rtype *_paths: dict
    :return *_paths: dictionaries for given pyatoa and seisflows pathing
    """
    # Get paths from Seisflows outputs .json file
    sf_json = os.path.join(parser.working_dir, "output", "seisflows_paths.json")
    with open(sf_json, "r") as f:
        sf_paths = json.load(f)
    
    # Pyatoa input/output directory, supplied by User in seisflows paths 
    pyatoa_io = sf_paths["PYATOA_IO"]

    # User defined configuration for Pyatoa, controlling processing params
    # pathing and switches for outputs and logging
    usrcfg_json = os.path.join(pyatoa_io, "sfconfig.json")
    with open(usrcfg_json, "r") as f:
        usrcfg = json.load(f)

    # Set Pyatoa paths, these is a hardcoded directory structure
    figs = os.path.join(pyatoa_io, "figures")
    vtks = os.path.join(figs, "vtks")
    data = os.path.join(pyatoa_io, "data")
    misfits = os.path.join(data, "misfits")
    misfit_file = os.path.join(pyatoa_io, "misfits.json")

    paths = {"PYATOA_FIGURES": figs, "PYATOA_DATA": data,
             "PYATOA_VTKS": vtks, "PYATOA_MISFITS": misfits,
             "MISFIT_FILE": misfit_file
             }

    # Processing requires extra process dependent paths
    if mode == "process":
        event_paths = {
            "ADJ_TRACES": os.path.join(parser.current_dir, "traces", "adj"),
            "SYN_TRACES": os.path.join(parser.current_dir, "traces", "syn"),
            "OBS_TRACES": os.path.join(parser.current_dir, "traces", "obs"),
            "EVENT_DATA": os.path.join(parser.current_dir, "DATA"),
            "STATIONS": os.path.join(parser.current_dir, "DATA", "STATIONS"),
            "EVENT_FIGURES": os.path.join(paths["PYATOA_FIGURES"],
                                          parser.model_number, parser.event_id
                                          )
        }
        # Make the event figure if necessary
        if not os.path.exists(event_paths["EVENT_FIGURES"]):
            os.makedirs(event_paths["EVENT_FIGURES"])

        paths = {**paths, **event_paths}

    return paths, usrcfg


def initialize(parser):
    """
    Create necessary directories for Pyatoa to operate with Seisflows

    :type parser: argparse.Parser
    :param parser: arguments from Seisflows
    """
    paths, _ = assemble_paths(parser)

    for key in paths.keys():
        # Don't make any filenames
        if "FILE" in key:
            continue
        elif not os.path.exists(paths[key]):
            os.makedirs(paths[key])


def finalize(parser):
    """
    Before finishing the iteration, create some final objects

    :type parser: argparse.Parser
    :param parser: arguments passed in from Seisflows
    """
    paths, usrcfg = assemble_paths(parser)

    # Generate .vtk files for given source and receivers
    if usrcfg["create_srcrcv_vtk"]:
        from pyatoa.utils.operations.file_generation import \
            create_srcrcv_vtk_multiple
        create_srcrcv_vtk_multiple(
            pathin=paths["PYATOA_DATA"], pathout=paths["PYATOA_VTKS"],
            model=parser.model_number
        )
    # Create copies of .h5 files at the end of each iteration, because .h5
    # files are easy to corrupt so it's good to have a backup
    if usrcfg["snapshot"]:
        snapshot_path = os.path.join(paths["PYATOA_DATA"], "snapshot")
        if not os.path.exists(snapshot_path):
            os.makedirs(snapshot_path)
        srcs = glob.glob(os.path.join(paths["PYATOA_DATA"], "*.h5"))
        for src in srcs:
            shutil.copy(src, os.path.join(snapshot_path, os.path.basename(src)))

    # Create misfit maps for each event with contour overlay showing misfit
    if usrcfg["plot_misfit_maps"]:
        from pyatoa.utils.visuals.mapping import event_misfit_map

        # Set the proper pathing so the misfit map can be found easily
        map_dir = os.path.join(paths["PYATOA_FIGURES"], "maps")
        if not os.path.exists(map_dir):
            os.makedirs(map_dir)

        step_number = "s{:0>2}".format(parser.step_count)
        name_template = "{eid}_{m}_{s}_misfit_map.png"

        # Loop through each available dataset to create misfit map
        datasets = glob.glob(os.path.join(paths["PYATOA_DATA"], "*.h5"))
        for dataset in datasets:
            with pyasdf.ASDFDataSet(dataset) as ds:
                fidout = os.path.join(
                    map_dir, name_template.format(
                        eid=ds.events[0].resource_id.id.split('/')[-1],
                        m=parser.model_number, s=step_number
                    )
                )
                # TO DO: set map corners using usrcfg, and not default Config?
                event_misfit_map(map_corners=pyatoa.Config().map_corners,
                                 ds=ds, model=parser.model_number,
                                 step=step_number,
                                 annotate_station_info='simple',
                                 contour_overlay=True, filled_contours=True,
                                 show=False, save=fidout
                                 )


def process(parser):
    """
    Main workflow calling on the core functionality of Pyatoa to process
    observed and synthetic waveforms and perform misfit quantification

    :type parser: argparse.Parser
    :param parser: arguments passed in from Seisflows
    """
    paths, usrcfg = assemble_paths(parser, mode="process")

    # Set loggging output, allow user to specify less output using 'info'
    if usrcfg["set_logging"]:
        logger_pyatoa = logging.getLogger("pyatoa")
        if usrcfg["set_logging"].lower() == "info":
            logger_pyatoa.setLevel(logging.INFO)
        else:
            logger_pyatoa.setLevel(logging.DEBUG)
            logger_pyflex = logging.getLogger("pyflex")
            logger_pyflex.setLevel(logging.DEBUG)

    # Set the Pyatoa Config object for misfit quantification
    config = pyatoa.Config(
        event_id=parser.event_id,
        model_number=parser.model_number,
        min_period=usrcfg["min_period"],
        max_period=usrcfg["max_period"],
        filter_corners=usrcfg["filter_corners"],
        rotate_to_rtz=usrcfg["rotate_to_rtz"],
        unit_output=usrcfg["unit_output"],
        pyflex_config=usrcfg["pyflex_config"],
        adj_src_type=usrcfg["adj_src_type"],
        synthetics_only=usrcfg["synthetics_only"],
        window_amplitude_ratio=usrcfg["window_amplitude_ratio"],
        cfgpaths={"synthetics": paths["SYN_TRACES"],
                  "waveforms": usrcfg["paths_to_waveforms"] +
                               [paths["OBS_TRACES"]],
                  "responses": usrcfg["paths_to_responses"]
                  }
        )

    # The trial step number is only used sparsely, e.g. the Statistics subgroup
    step_number = "s{:0>2}".format(parser.step_count)

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

        # Station file has the form NET STA LAT LON ...
        stations = np.loadtxt(paths["STATIONS"], usecols=[0, 1, 2, 3],
                              dtype=str)
        coords = stations[:, 2:]

        # Loop through stations and invoke Pyatoa workflow
        for station in stations:
            sta, net = station[:2]
            print("{}.{}".format(net, sta))
            try:
                mgmt.reset()

                # Gather data, searching internal pathways, else fetching from
                # external pathways if possible. Preprocess identically
                mgmt.gather_data(station_code="{net}.{sta}.{loc}.{cha}".format(
                                 net=net, sta=sta, loc="*", cha="HH*")
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
                    mgmt.windows = misfit_windows

                mgmt.run_pyadjoint()

                # Plot waveforms with misfit windows and adjoint sources
                if usrcfg["plot_waveforms"]:
                    # Format some strings to append to the waveform plot title
                    append_title = (
                        "\n{md}{sn} pyflex={pf}, pyadjoint={pa},".format(
                            md=config.model_number, sn=step_number, 
                            pf=config.pyflex_config[0],
                            pa=config.pyadjoint_config[0])
                    )
                    if mgmt.total_misfit is not None:
                        append_title = " ".join([
                            append_title,
                            "misfit={:.2E}".format(mgmt.total_misfit)]
                        )
                    f = mgmt.plot_wav(
                        append_title=append_title,
                        save=os.path.join(
                                  paths["EVENT_FIGURES"], "wav_{}".format(sta)),
                        show=False, return_figure=True
                        )

                # Plot source-receiver maps, don't make a map if no wav data
                if usrcfg["plot_maps"] and f:
                    mgmt.plot_map(
                        stations=coords, save=os.path.join(
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
        write_stats_to_asdf(ds, config.model_number, parser.step_count)

        # Create the .sem ascii files required by specfem
        write_adj_src_to_ascii(ds, config.model_number, paths["ADJ_TRACES"])

        # Create the STATIONS_ADJOINT file required by specfem
        create_stations_adjoint(ds, config.model_number,
                                specfem_station_file=paths["STATIONS"],
                                pathout=paths["EVENT_DATA"])

        # Write misfits for seisflows into individual text files
        write_misfit_stats(ds, config.model_number, paths["PYATOA_MISFITS"])

        # Sum and write misfits information to a JSON file
        write_misfit_json(ds, parser.model_number, parser.step_count,
                          paths["MISFIT_FILE"])

        # Combine .png images into a composite .pdf for easy fetching
        if usrcfg["tile_and_combine"]:
            from pyatoa.utils.visuals.convert_images import tile_and_combine
            tile_and_combine(ds=ds, model=parser.model_number, step=step_number,
                             figure_path=paths["PYATOA_FIGURES"],
                             purge_originals=usrcfg["purge_originals"],
                             purge_tiles=usrcfg["purge_tiles"]
                             )


if __name__ == "__main__":
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
        # Run Pyatoa, return successful exit code
        try:
            process(parser)
            sys.exit(0)
        except Exception as e:
            traceback.print_exc()
            sys.exit(1)

    else:
        print("invalid 'mode' argument")
        sys.exit(1)
