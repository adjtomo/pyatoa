#!/usr/bin/env python3
"""
Pyaflowa - Pyatoa's Seisflows plugin class

This Seisflows plugin class that allows easy scripting of Pyatoa
functionality into a Seisflows workflow. Pre-written functionalities simplify
calls made in Seisflows to Pyatoa, to reduce clutter inside the workflow.
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
import numpy as np

from pyatoa import logger
from pyatoa.utils.form import model_number, step_count
from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.visuals.statistics import plot_output_optim
from pyatoa.utils.write import (create_stations_adjoint, write_adj_src_to_ascii, 
                                write_misfit_stats, tile_combine_imgs, 
                                src_vtk_from_specfem, rcv_vtk_from_specfem
                                )

# Overwrite the preprocessing function
from pyatoa.plugins.new_zealand.process import preproc

# A list of key word arguments that are accepted by Pyaflowa but are only
# listed in Seisflows' parameters.yaml file. Listed here so that Pyatoa
# knows that these arguments are acceptable.
pyaflowa_kwargs = ["set_logging", "window_amp_ratio", "fix_windows", "snapshot",
                   "write_misfit_json", "create_srcrcv_vtk", "plot_waveforms",
                   "plot_srcrcv_maps", "combine_imgs", "inspect"]


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
        self.stats_dir = oj(pyatoa_io, "figures", "stats")

        # Create Pyatoa directory structure
        for fid in [self.figures_dir, self.data_dir, self.misfits_dir,
                    self.maps_dir, self.vtks_dir, self.snapshots_dir]:
            if not os.path.exists(fid):
                os.makedirs(fid)

        # Set some attributes that will be set/used during the workflow
        self.iteration = 0
        self.step_count = 0
        self.synthetics_only = bool(par["CASE"].lower() == "synthetic")

        self._check(par)

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

    def _check(self, ext_par):
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

    def fixwin(self, ds):
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
        Convenience function to easily set multiple parameters before calling
        other functions.

        Overwrite internally used attributes using kwargs. Ensure that
        attributes other than the ones set in __init__ are allowed.
        """
        for key in list(kwargs.keys()):
            if not hasattr(self, key):
                logger.warning(f"Pyaflowa has no attribute '{key}', ignoring")
                del kwargs[key]
        self.__dict__.update(kwargs)

    def setup(self, cwd, event_id=None):
        """
        Set up one embarassingly parallelizable workflow by creating individual
        process dependent pathways, and creating indvidual Pyatoa Config objects
        that will control the worklow.

        :type cwd: str
        :param cwd: current working directory for this instance of Pyatoa
        :type event_id: str
        :param event_id: event identifier tag for file naming etc.
        """
        # Default event id is the name of the current working directory
        if event_id is None:
            event_id = os.path.basename(cwd)

        # Process specific internal directories for the processing
        ev_paths = {"stations": oj(cwd, "DATA", "STATIONS"),
                    "maps": oj(self.maps_dir, event_id),
                    "figures": oj(self.figures_dir, self.model, event_id),
                    }

        # Create the process specific event directories
        for key, item in ev_paths.items():
            if not os.path.exists(item):
                os.makedirs(item)

        # Set logging output level for all packages
        for log, level in self.set_logging.items():
            if level:
                logger_ = logging.getLogger(log)
                if level == "info":
                    logger_.setLevel(logging.INFO)
                elif level == "debug":
                    logger_.setLevel(logging.DEBUG)

        # Read in the Pyatoa Config object from the .yaml file, with
        # additional parameter set by the individual process
        config = pyatoa.Config(
            yaml_fid=self.config_file, event_id=event_id, model=self.model,
            step=self.step, synthetics_only=self.synthetics_only,
            cfgpaths={"synthetics": oj(cwd, "traces", "syn"),
                      "waveforms": oj(cwd, "traces", "obs")}
        )

        return config, ev_paths

    def process(self, cwd, event_id=None, overwrite=None):
        """
        Main workflow calling on the core functionality of Pyatoa to process
        observed and synthetic waveforms and perform misfit quantification.

        :type cwd: str
        :param cwd: current working directory for this instance of Pyatoa
        :type event_id: str
        :param event_id: event identifier tag for file naming etc.
        """
        # Run the setup and standardize some names
        config, ev_paths = self.setup(cwd, event_id)
        ds_name = oj(self.data_dir, f"{config.event_id}.h5")

        # Count number of successful processes
        processed = 0
        with pyasdf.ASDFDataSet(ds_name) as ds:
            fix_windows = self.fixwin(ds)
            logger.info(f"Fix windows: {fix_windows}")

            # Make sure the ASDFDataSet doesn't already contain auxiliary_data
            # for the model_number/step_count
            clean_ds(ds=ds, model=self.model, step=self.step)

            # Set up the manager and get station information
            config.write(write_to=ds)
            mgmt = pyatoa.Manager(config=config, ds=ds)
            stations = np.loadtxt(ev_paths["stations"], usecols=[0, 1, 2, 3],
                                  dtype=str)
            coords = stations[:, 2:]

            # Loop through stations and invoke Pyatoa workflow
            for station in stations:
                sta, net = station[:2]
                logger.info(f"{net}.{sta}")
                try:
                    mgmt.reset()
                    mgmt.gather(station_code=f"{net}.{sta}.*.HH*")
                    mgmt.standardize()
                    mgmt.preprocess(overwrite=overwrite)
                    mgmt.window(fix_windows=fix_windows)
                    mgmt.measure()

                    # Plot waveforms with misfit windows and adjoint sources
                    if self.plot_waveforms:
                        # Format some strings to append to the waveform title
                        tit = " ".join([
                            f"\n{config.model}{self.step}",
                            f"pyflex={config.pyflex_preset},",
                            f"pyadjoint={config.adj_src_type},",
                            f"misfit={mgmt.misfit:.2E}"
                        ])
                        mgmt.plot(append_title=tit,
                                  save=oj(ev_paths["figures"], f"wav_{sta}"),
                                  show=False, return_figure=False
                                  )

                    # Only plot maps once since they won't change
                    if self.plot_srcrcv_maps:
                        map_fid = oj(ev_paths["maps"], f"map_{sta}.png")
                        if not os.path.exists(map_fid):
                            mgmt.srcrcvmap(stations=coords, show=False,
                                           save=map_fid)

                    processed += 1
                # Use traceback ensures more detailed error tracking
                except Exception:
                    logger.debug(traceback.print_exc())
                    continue

            # Run finalization procedures for processing iff gathered waveforms
            if processed:
                logger.info(f"Pyaflowa processed {processed} stations")
                self.finalize_process(ds=ds, cwd=cwd, ev_paths=ev_paths,
                                      config=config)
            else:
                logger.info("Pyaflowa processed 0 stations, skipping finalize")

    def finalize_process(self, cwd, ds, ev_paths, config):
        """
        After all waveforms have been windowed and measured, run some functions
        that create output files useful for Specfem, or for the User.

        :type cwd: str
        :param cwd: current working directory of solver
        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset contianing the waveforms and misfit for this solver
        :type ev_paths: dict
        :param ev_paths: dictionary of event/solver specific paths
        :type config: pyatoa.core.config.Config
        :param config: Pyatoa config object containing parameters needed for
            finalization of workflow
        """
        logger.info("writing adjoint sources")
        write_adj_src_to_ascii(ds, self.model, self.step,
                               oj(cwd, "traces", "adj"))

        logger.info("creating STATIONS_ADJOINT")
        create_stations_adjoint(ds, model=self.model, 
                                step=self.step,
                                specfem_station_file=ev_paths["stations"],
                                pathout=oj(cwd, "DATA")
                                )

        logger.info("writing event misfit to disk")
        write_misfit_stats(ds, self.model, self.step, pathout=self.misfits_dir)

        # Combine images into a pdf for easier visualization, will delete .png's
        if self.combine_imgs:
            logger.info("creating composite pdf")
            # path/to/figures/m??/s??/eid_m??_s??.pdf
            fig_path = oj(self.figures_dir, self.model, self.step)
            fid = f"{config.event_id}_{self.model}{self.step}.pdf"
            if not os.path.exists(fig_path):
                os.makedirs(fig_path)
            tile_combine_imgs(ds=ds, save_pdf_to=oj(fig_path, fid),
                              wavs_path=ev_paths["figures"],
                              maps_path=ev_paths["maps"],
                              purge_wavs=True, purge_tiles=True
                              )
            # remove the empty event directory which has been purged
            if not glob.glob(oj(ev_paths["figures"], "*")):
                os.rmdir(ev_paths["figures"])

    def finalize(self):
        """
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
        if self.create_srcrcv_vtk and self.iteration == 1:
            for func in [src_vtk_from_specfem, rcv_vtk_from_specfem]:
                func(path_to_data=self.specfem_data, path_out=self.vtks_dir)

        # Run the Inspector class to analyze the misfit behavior of inversion
        if self.inspect:
            insp = pyatoa.Inspector(path=self.data_dir)
            insp.save(tag=self.title, path=self.data_dir)

            # Create a misfit histogram for the initial and final model
            for choice, binsize in zip(["cc_shift_sec", "dlna"], [0.5, 0.25]):
                insp.misfit_histogram(model=insp.models[0], choice=choice,
                                      model_comp=insp.models[-1], show=False,
                                      binsize=binsize,
                                      save=oj(
                                          self.stats_dir,
                                          f"misfithisto_{insp.models[-1]}.png")
                                      )


