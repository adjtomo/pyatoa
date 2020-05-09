#!/usr/bin/env python3
"""
Pyaflowa - Pyatoa's Seisflows plugin class

This Seisflows plugin class that allows easy scripting of Pyatoa functionality
into a Seisflows workflow. It takes care of input from the Seisflows master
directory, and output into a Pyatoa-specific auxiliary directory with a
pre-defined directory structure that allows for easy navigation of outputs.
"""
import os
# cheeky shorthand for cleaner calls
from os.path import join as oj
import sys
import glob
import pyasdf
import pyatoa
import shutil
import logging
import traceback

from pyatoa import logger
from pyatoa.utils.read import read_stations
from pyatoa.utils.form import model_number, step_count, event_name
from pyatoa.utils.asdf.clean import clean_ds
from pyatoa.visuals.statistics import plot_output_optim
from pyatoa.utils.write import (write_stations_adjoint, write_adj_src_to_ascii,
                                write_misfit, tile_combine_imgs,
                                src_vtk_from_specfem, rcv_vtk_from_specfem
                                )

# A list of key word arguments that are accepted by Pyaflowa but are only
# listed in Seisflows' parameters.yaml file. Listed here so that Pyatoa
# knows that these arguments are acceptable.
pyaflowa_kwargs = ["set_logging", "win_amp_ratio", "fix_windows", "snapshot",
                   "srcrcv_vtk", "plot_wav", "plot_map", "make_pdf"]


class Pyaflowa:
    """
    The plugin object that is created to exist within Seisflows, keep track of
    the Seisflows workflow, and create the necessary outputs when requested
    """
    def __init__(self, par, paths):
        """
        Pyaflowa only needs to know what Seisflows knows.
        With this information it can create the internal directory
        structure that it uses to navigate around Seisflows.

        :type par: dict
        :param par: a dictionary of the Seisflows parameters contained in the
            `PAR` variable. should be passed here as vars(PAR)
        :type paths: dict
        :param paths: a dictionary of the Seisflows paths contained in the
            `PATH` variable. should be passed here as vars(PATH)
        """
        # Distribute the relative Seisflows paramaters to Pyaflowa
        self.__dict__ = par["PYATOA"]
        self.title = par["TITLE"]

        # Grab relevant external paths
        assert("PYATOA_IO" in paths.keys())
        pyatoa_io = paths["PYATOA_IO"]
        self.work_dir = paths["WORKDIR"]
        self.specfem_data = paths["SPECFEM_DATA"]
        self.config_file = oj(self.work_dir, "parameters.yaml")

        # Distribute internal paths
        self.data_dir = oj(pyatoa_io, "data")
        self.misfits_dir = oj(pyatoa_io, "data", "misfits")
        self.snapshots_dir = oj(pyatoa_io, "data", "snapshot")
        self.figures_dir = oj(pyatoa_io, "figures")
        self.maps_dir = oj(pyatoa_io, "figures", "maps")
        self.vtks_dir = oj(pyatoa_io, "figures", "vtks")

        # Create Pyatoa directory structure
        for fid in [self.figures_dir, self.data_dir, self.misfits_dir,
                    self.maps_dir, self.vtks_dir, self.snapshots_dir]:
            if not os.path.exists(fid):
                os.makedirs(fid)

        # Set some attributes that will be set/used during the workflow
        self.iteration = 0
        self.step_count = 0
        self.synthetics_only = bool(par["CASE"].lower() == "synthetic")

        # Set logging output level for all packages within Pyatoa
        for log, level in self.set_logging.items():
            if level:
                logger_ = logging.getLogger(log)
                if level == "info":
                    logger_.setLevel(logging.INFO)
                elif level == "debug":
                    logger_.setLevel(logging.DEBUG)

        self._check_parameters(par)

    def __str__(self):
        """
        String representation, only really used in Seisflows debug mode so
        not very pretty, but should have some useful information
        """
        # An ugly way to format the first two lines the same as the rest
        str_out = "\n".join([
            "PYAFLOWA", "\t{:<25}{}".format("model:", self.model),
            "\t{:<25}{}".format("step:", self.step), ""]
        )
        for key, item in vars(self).items():
            if isinstance(item, dict):
                continue
            str_out += f"\t{key+':':<25}{item}\n"
        return str_out

    @property
    def model(self):
        """
        The model number is based on the current iteration, which starts at 1
        so model number, starting from 0, lags behind by 1. Formatted like 'm00'
        """
        return model_number(max(self.iteration - 1, 0))

    @property
    def step(self):
        """
        Formatted step count, e.g. 's00'
        """
        return step_count(self.step_count)

    def _check_parameters(self, ext_par):
        """
        Perform some sanity checks upon initialization. If they fail, hard exit
        so that Seisflows crashes, that way things don't crash after jobs have
        been submitted etc.

        :type ext_par: dict
        :param ext_par: parameter dictionary from Seisflows
        """
        # Ensure that the gathered seismogram length is greater than the
        # length of synthetics
        if ext_par["DT"] * ext_par["NT"] >= self.start_pad + self.end_pad:
            logger.warning("length of gathered observed waveforms will be less "
                           "than the length of synthetics... exiting")
            sys.exit(-1)

        # Check that fix windows parameter set correctly
        assert(self.fix_windows in [True, False, "iter"]), \
            "fix_windows must be bool or 'iter'"

        # Try to feed the parameter file into Config to see if it throws
        # any ValueErrors from incorrect arguments
        try:
            pyatoa.Config(yaml_fid=self.config_file)
        except ValueError as e:
            logger.warning(e)
            logger.warning("Config encountered unexpected arguments... exiting")
            sys.exit(-1)

    def _check_for_fixed_windows(self, ds):
        """
        Determine if window fixing is required. This can be done by step count,
        by iteration, or not at all.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to check for misfit windows
        :rtype: bool
        :return: if fix window True, don't reevaluate window unless m00s00
                 if fix by iteration, return False if s00, True otherwise
                 if fix window False, reevaluate windows each step
        """
        if isinstance(self.fix_windows, bool):
            return self.fix_windows
        elif self.fix_windows == "iter":
            # False if this is the first step in the iteration, True else
            try:
                return hasattr(ds.auxiliary_data.MisfitWindows, self.model)
            except AttributeError:
                # If no MisfitWindows aux data, ds has just been instantiated
                return False
        else:
            raise TypeError("fix_windows should be bool or 'iter'")

    def set(self, **kwargs):
        """
        High level function to interact with Seisflows.

        Set internal parameters before calling other functions. Allows Seisflows
        to communicate where in the inversion it is before calling Pyaflowa.

        Overwrite internally used attributes using kwargs. Ensure that
        attributes other than the ones set in __init__ are allowed.
        """
        for key in list(kwargs.keys()):
            if not hasattr(self, key):
                logger.warning(f"Pyaflowa has no attribute '{key}', ignoring")
                del kwargs[key]
        self.__dict__.update(kwargs)

    def eval_func(self, cwd, event_id, overwrite=None):
        """
        High level function to interact with Seisflows.

        Reads data, applies preprocessing, writes residuals and adjoint sources.

        :type cwd: str
        :param cwd: path to the seisflows Solver current working directory
        :type event_id: str
        :param event_id: source name used for labelling data from `cwd`
        :type overwrite: function
        :param overwrite: preprocessing overwrite function that can be passed
            from Seisflows, if None, default Pyatoa preproc function used.
        """
        # Set up the machinery for a single workflow instnace
        config, paths = self.prepare_event(cwd, event_id)
        with pyasdf.ASDFDataSet(paths["dataset"]) as ds:
            status = self.process_event(ds, config, paths, overwrite)
            if status:
                self.prepare_eval_grad(ds)

                if self.make_pdf:
                    self.make_event_pdf(ds, paths)

                logger.info("writing event misfit to disk")
                write_misfit(ds, self.model, self.step, path=self.misfits_dir)

    def finalize(self):
        """
        High level function to interact with Seisflows.

        Finalization function to be run after each iteration to allow
        Pyaflowa to clean up intermediate files and create any optional
        output files

        At the end of an iteration, clean up working directory and create final
        objects if requested by the User. This includes statistical plots
        VTK files for model visualizations, and backups of the data.
        """
        # Delete the temporary stored misfit data to avoid over-printing
        for fid in glob.glob(oj(self.misfits_dir, "*")):
            os.remove(fid)

        # Create copies of .h5 files at the end of each iteration, because .h5
        # can be corrupted if open during crashes so it's good to have a backup
        if self.snapshot:
            srcs = glob.glob(oj(self.data_dir, "*.h5"))
            for src in srcs:
                shutil.copy(src, oj(self.snapshots_dir, os.path.basename(src)))

        # Plot the output.optim file outputted by Seisflows
        plot_output_optim(path_to_optim=oj(self.work_dir, "output.optim"),
                          save=oj(self.figures_dir, "output_optim.png")
                          )

        # Generate .vtk files for given source and receivers for model 0
        if self.srcrcv_vtk and self.iteration == 1:
            for func in [src_vtk_from_specfem, rcv_vtk_from_specfem]:
                func(path_to_data=self.specfem_data, path_out=self.vtks_dir)

    def prepare_event(self, cwd, event_id):
        """
        Mid level function to set up an embarassingly parallelizable workflow
        instance by establishing event-dependent directory structure.
        Creating matching Pyatoa Config object.

        :type cwd: str
        :param cwd: current working directory for this instance of Pyatoa
        :type event_id: str
        :param event_id: event identifier tag for file naming etc.
        :rtype config: pyatoa.core.Config
        :return config: Configuration object to control the workflow
        :rtype paths: dict
        :return paths: event-specific paths used for I/O
        """
        # Process specific internal directories for the processing
        paths = {"cwd": cwd,
                 "dataset": oj(self.data_dir, f"{event_id}.h5"),
                 "maps": oj(self.maps_dir, event_id),
                 "figures": oj(self.figures_dir, self.model, event_id),
                 }
        for key in ["maps", "figures"]:
            if not os.path.exists(paths[key]):
                os.makedirs(paths[key])

        # Config object from .yaml file, event specific trace directories
        config = pyatoa.Config(
            yaml_fid=self.config_file, event_id=event_id, model=self.model,
            step=self.step, synthetics_only=self.synthetics_only,
            cfgpaths={"synthetics": oj(cwd, "traces", "syn"),
                      "waveforms": oj(cwd, "traces", "obs")}
        )

        return config, paths

    def process_event(self, ds, config, paths, overwrite=None):
        """
        Mid-level functionality to gather, preprocess data for a given dataset

        Main workflow calling on the core functionality of Pyatoa to process
        observed and synthetic waveforms and perform misfit quantification.

        :type ds: pyasdf.ASDFDataSet
        :param ds: the dataset that will be used for collecting and storing data
        :type config: pyatoa.core.Config
        :param config: Configuration object to control the workflow
        :type paths: dict
        :param paths: event-specific paths used for I/O
        :type overwrite: function
        :param overwrite: a preprocessing function to overwite the default
            Pyatoa preprocessing function. Must be specified in Seisflows.
        """
        # Count number of successful processes
        processed = 0
        fix_windows = self._check_for_fixed_windows(ds)
        logger.info(f"Fix windows: {fix_windows}")

        # Make sure the ASDFDataSet doesn't already contain auxiliary_data
        clean_ds(ds=ds, model=self.model, step=self.step)

        # Set up the Manager and get station information
        config.write(write_to=ds)
        mgmt = pyatoa.Manager(config=config, ds=ds)
        inv = read_stations(oj(paths["cwd"], "DATA", "STATIONS"))

        # Loop through stations and invoke Pyatoa workflow
        for net in inv:
            for sta in net:
                logger.info(f"{net.code}.{sta.code}")
                try:
                    processed += mgmt.flow(
                            station_code=f"{net.code}.{sta.code}.*.HH*",
                            preprocess_overwrite=overwrite,
                            fix_windows=fix_windows
                            )
                    if self.plot_wav:
                        mgmt.plot(save=oj(paths["figures"], f"wav_{sta.code}"),
                                  show=False, return_figure=False
                                  )
                    if self.plot_map:
                        # Only plot maps once since they won't change
                        map_fid = oj(paths["maps"], f"map_{sta.code}.png")
                        if not os.path.exists(map_fid):
                            mgmt.srcrcvmap(stations=inv, show=False,
                                           save=map_fid)
                except pyatoa.ManagerError:
                    continue
                except Exception:
                    traceback.print_exc()
                    continue
        logger.info(f"Pyaflowa processed {processed} stations")
        return bool(processed)

    def prepare_eval_grad(self, ds, paths):
        """
        Prepare for gradient evaluation by exporting the adjoint traces from the
        ASDFDataSet and generating a stations file required for Specfem's
        adjoint simulations

        :type ds: pyasdf.ASDFDataSet
        :param ds: the dataset that will be used for collecting and storing data
        :type paths: dict
        :param paths: event-specific paths used for I/O
        """
        logger.info("writing adjoint sources")
        write_adj_src_to_ascii(ds=ds, model=self.model, step=self.step,
                               pathout=oj(paths["cwd"], "traces", "adj")
                               )

        logger.info("creating STATIONS_ADJOINT")
        write_stations_adjoint(
            ds=ds, model=self.model, step=self.step,
            specfem_station_file=oj(paths["cwd"], "DATA", "STATIONS"),
            pathout=oj(paths["cwd"], "DATA")
        )

    def make_event_pdf(self, ds, paths, purge=True):
        """
        Make a PDF of waveforms and maps tiled together for easy analysis of
        waveform fits, default purge the original files for space-saving

        :type ds: pyasdf.ASDFDataSet
        :param ds: the dataset that will be used for collecting and storing data
        :type paths: dict
        :param paths: event-specific paths used for I/O
        :type purge: bool
        :param purge: delete the orginal waveform and intermediate tile .png
            files and only retain the resultant PDF. Optional but hidden as
            purging is preferable to avoid too many files.
        """
        event_id = event_name(ds)

        # Establish correct directory and file name
        # path/to/figures/m??/s??/eid_m??_s??.pdf
        outfile = oj(self.figures_dir, self.model, self.step,
                     f"{event_id}_{self.model}{self.step}.pdf")

        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))

        tile_combine_imgs(
            ds=ds, save_pdf_to=outfile, wavs_path=paths["figures"],
            maps_path=paths["maps"], purge_wavs=purge, purge_tiles=purge
        )

        # if purged, remove the empty event directory
        if not glob.glob(oj(paths["figures"], "*")):
            os.rmdir(paths["figures"])

