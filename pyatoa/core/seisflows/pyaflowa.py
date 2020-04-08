#!/usr/bin/env python3
"""
Pyaflowa

The Seisflows plugin class that allows easy scripting of Pyatoa
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
from pyatoa.utils.form import model, step
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
        # Ensure that the inputs are accessible by the class
        self.par = par["PYATOA"]
        self.ext_paths = paths
        self.ext_par = par

        # Distribute internal hardcoded path structure
        assert("PYATOA_IO" in self.ext_paths.keys())
        pyatoa_io = self.ext_paths["PYATOA_IO"]

        # Tag the external files that will need to be used throughout
        self.config_file = oj(self.ext_paths["WORKDIR"], "parameters.yaml")
        self.misfit_file = oj(pyatoa_io, "misfits.json")

        # Distribute internal paths
        self.data_dir = oj(pyatoa_io, "data")
        self.misfits_dir = oj(pyatoa_io, "data", "misfits")
        self.snapshots_dir = oj(pyatoa_io, "data", "snapshot")
        self.figures_dir = oj(pyatoa_io, "figures")
        self.maps_dir = oj(pyatoa_io, "figures", "maps")
        self.vtks_dir = oj(pyatoa_io, "figures", "vtks")
        self.stats_dir = oj(pyatoa_io, "figures", "stats")
        self.composites_dir = oj(pyatoa_io, "figures", "composites")

        # Create Pyatoa directory structure
        for fid in [self.figures_dir, self.data_dir, self.misfits_dir,
                    self.maps_dir, self.vtks_dir, self.composites_dir,
                    self.snapshots_dir]:
            if not os.path.exists(fid):
                os.makedirs(fid)

        # Set some attributes that will be set/used during the workflow
        self.iteration = 0
        self.step = 0
        self.fix_windows = self.par["fix_windows"]
        self.synthetics_only = bool(par["CASE"].lower() == "synthetic")

        self._check()

    def __str__(self):
        """
        String representation, only really used in Seisflows debug mode so
        not very pretty, but should have some useful information
        """
        # An ugly way to format the first two lines the same as the rest
        str_out = "\n".join([
            "PYAFLOWA", "\t{:<25}{}".format("model_number:", self.model_number),
            "\t{:<25}{}".format("step_count:", self.step_count), ""]
        )
        for key, item in vars(self).items():
            if isinstance(item, dict):
                continue
            str_out += f"\t{key+':':<25}{item}\n"
        return str_out

    def __repr__(self):
        """Simple point to __str__()"""
        return self.__str__()

    @property
    def model_number(self):
        """
        The model number is based on the current iteration, which starts at 1
        so model number, starting from 0, lags behind by 1.
        """
        return model(max(self.iteration - 1, 0))

    @property
    def step_count(self):
        """
        Formatted step count, e.g. 's00'
        """
        return step(self.step)

    def _check(self):
        """
        Perform some sanity checks upon initialization. If they fail, hard exit
        so that Seisflows crashes, that way things don't crash after jobs have
        been submitted etc.
        """
        # Ensure that the gathered seismogram length is greater than the
        # length of synthetics
        if (self.ext_par["DT"] * self.ext_par["NT"] >=
                self.par['start_pad'] + self.par['end_pad']):
            logger.warning("length of gathered observed waveforms will be less "
                           "than the length of synthetics... exiting")
            sys.exit(-1)

        # Try to feed the parameter file into Config to see if it throws
        # any ValueErrors from incorrect arguments
        try:
            pyatoa.Config(yaml_fid=self.config_file)
        except ValueError as e:
            logger.warning(e)
            logger.warning("Config encountered unexpected arguments... exiting")
            sys.exit(-1)

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

    def setup_process(self, cwd, event_id=None):
        """
        Set up the workflow by creating process dependent pathways, and creating
        the Pyatoa Config object that will control the worklow

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
                    "figures": oj(self.figures_dir, self.model_number,
                                  event_id),
                    }

        # Create the process specific event directories
        for key, item in ev_paths.items():
            if not os.path.exists(item):
                os.makedirs(item)

        # Set logging output level for all packages
        for log, level in self.par["set_logging"].items():
            if level:
                logger_ = logging.getLogger(log)
                if level == "info":
                    logger_.setLevel(logging.INFO)
                elif level == "debug":
                    logger_.setLevel(logging.DEBUG)

        # Read in the Pyatoa Config object from the .yaml file, with
        # additional parameter set by the individual process
        config = pyatoa.Config(
            yaml_fid=self.config_file, event_id=event_id,
            model_number=self.model_number, step_count=self.step_count,
            synthetics_only=self.synthetics_only,
            cfgpaths={"synthetics": oj(cwd, "traces", "syn"),
                      "waveforms": oj(cwd, "traces", "obs")}
        )

        return config, ev_paths

    def process(self, cwd, event_id=None):
        """
        Main workflow calling on the core functionality of Pyatoa to process
        observed and synthetic waveforms and perform misfit quantification

        :type cwd: str
        :param cwd: current working directory for this instance of Pyatoa
        :type event_id: str
        :param event_id: event identifier tag for file naming etc.
        """
        # Run the setup and standardize some names
        config, ev_paths = self.setup_process(cwd, event_id)
        ds_name = oj(self.data_dir, f"{config.event_id}.h5")

        # Count number of successful processes
        processed, config_written = 0, False
        with pyasdf.ASDFDataSet(ds_name) as ds:
            # Make sure the ASDFDataSet doesn't already contain auxiliary_data
            # for the model_number/step_count
            clean_ds(ds=ds, model=self.model_number, step=self.step_count,
                     fix_windows=self.fix_windows)

            # Set up the manager and get station information
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
                    mgmt.preprocess(overwrite=preproc)
                    mgmt.window(fix_windows=self.fix_windows)
                    mgmt.measure()

                    # Plot waveforms with misfit windows and adjoint sources
                    if self.par["plot_waveforms"]:
                        # Format some strings to append to the waveform title
                        tit = " ".join([
                            f"\n{config.model_number}{self.step_count}",
                            f"pyflex={config.pyflex_preset},",
                            f"pyadjoint={config.adj_src_type},",
                            f"misfit={mgmt.misfit:.2E}"
                        ])
                        mgmt.plot(append_title=tit,
                                  save=oj(ev_paths["figures"], f"wav_{sta}"),
                                  show=False, return_figure=False
                                  )

                    # Only plot maps once since they won't change
                    if self.par["plot_srcrcv_maps"] and \
                            self.model_number == "m00" and \
                            self.step_count == "s00":
                        mgmt.srcrcvmap(stations=coords, show=False,
                                       save=oj(ev_paths["maps"], f"map_{sta}"))

                    # Just once, grab the processing stats from the Streams and
                    # append them to the Config object and save. A sort of
                    # hacky way to retain processing information from old runs.
                    if not config_written and mgmt.st_obs is not None:
                        setattr(config, "obs_processing",
                                mgmt.st_syn[0].stats.processing)
                        setattr(config, "syn_processing",
                                mgmt.st_obs[0].stats.processing)
                        config.write(write_to=ds)
                        config_written = True

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
        write_adj_src_to_ascii(ds, self.model_number, self.step_count, 
                               oj(cwd, "traces", "adj"))

        logger.info("creating STATIONS_ADJOINT")
        create_stations_adjoint(ds, model=self.model_number, 
                                step=self.step_count,
                                specfem_station_file=ev_paths["stations"],
                                pathout=oj(cwd, "DATA")
                                )

        logger.info("writing event misfit to disk")
        write_misfit_stats(ds, self.model_number, self.step_count, 
                           pathout=self.misfits_dir)

        # Combine images into a pdf for easier visualization
        if self.par["combine_imgs"]:
            # path/to/eventid_modelstep_wavmap.png
            logger.info("creating composite pdf")
            save_to = oj(self.composites_dir,
                         f"{config.event_id}_{self.model_step}.pdf")
            tile_combine_imgs(ds=ds, save_pdf_to=save_to,
                              wavs_path=ev_paths["figures"],
                              maps_path=ev_paths["maps"],
                              purge_wavs=True, purge_tiles=True
                              )

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
        if self.par["snapshot"]:
            srcs = glob.glob(oj(self.data_dir, "*.h5"))
            for src in srcs:
                shutil.copy(src, oj(self.snapshots_dir, os.path.basename(src)))

        # Plot the output.optim file outputted by Seisflows
        plot_output_optim(path_to_optim=oj(self.ext_paths["WORKDIR"],
                                           "output.optim"),
                          save=oj(self.figures_dir, "output_optim.png")
                          )

        # Generate .vtk files for given source and receivers for model 0
        if self.par["create_srcrcv_vtk"] and self.iteration == 1:
            for func in [src_vtk_from_specfem, rcv_vtk_from_specfem]:
                func(path_to_data=self.ext_paths["SPECFEM_DATA"],
                     path_out=self.vtks_dir)

        # Run the Inspector class to analyze the misfit behavior of inversion
        if self.par["inspect"]:
            insp = pyatoa.Inspector(path=self.data_dir)
            insp.save(tag=f"{self.ext_par['TITLE']}", path=self.data_dir)

            # Create a misfit histogram for the initial and final model
            for choice, binsize in zip(["cc_shift_sec", "dlna"], [0.5, 0.25]):
                insp.misfit_histogram(model=insp.models[0], choice=choice,
                                      model_comp=insp.models[-1], show=False,
                                      binsize=binsize,
                                      save=oj(
                                          self.stats_dir,
                                          f"misfithisto_{insp.models[-1]}.png")
                                      )


