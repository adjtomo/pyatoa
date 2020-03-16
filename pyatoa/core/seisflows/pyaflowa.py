#!/usr/bin/env python3
"""
Pyaflowa

The Seisflows plugin class that allows easy scripting of Pyatoa
functionality into a Seisflows workflow. Pre-written functionalities simplify
calls made in Seisflows to Pyatoa, to reduce clutter inside the workflow
"""
import os
import glob
import pyasdf
import pyatoa
import shutil
import logging
import warnings
import traceback
import numpy as np

from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.asdf.additions import write_stats_to_asdf
from pyatoa.visuals.statistics import plot_output_optim
from pyatoa.utils.io import (create_stations_adjoint, write_misfit_json,
                             write_adj_src_to_ascii, write_misfit_stats,
                             tile_combine_imgs,
                             create_srcrcv_vtk_multiple)

# Overwrite the preprocessing function
from pyatoa.plugins.new_zealand.process import preproc


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
        # Ensure that necessary inputs are accessible by the class
        self.ext_paths = paths
        self.par = par["PYATOA"]

        # Distribute internal hardcoded path structure
        assert("PYATOA_IO" in self.ext_paths.keys())
        pyatoa_io = self.ext_paths["PYATOA_IO"]

        self.int_paths = {
            "config_file": os.path.join(self.ext_paths["WORKDIR"],
                                        "parameters.yaml"),
            "misfit_file": os.path.join(pyatoa_io, "misfits.json"),
            "figures": os.path.join(pyatoa_io, "figures"),
            "data": os.path.join(pyatoa_io, "data"),
            "misfits": os.path.join(pyatoa_io, "data", "misfits"),
            "maps": os.path.join(pyatoa_io, "figures", "maps"),
            "vtks": os.path.join(pyatoa_io, "figures", "vtks"),
            "composites": os.path.join(pyatoa_io, "figures", "composites"),
            "snapshots": os.path.join(pyatoa_io, "data", "snapshot"),
            }

        # Create Pyatoa directory structure
        for key, item in self.int_paths.items():
            if "file" not in key:
                if not os.path.exists(item):
                    os.makedirs(item)

        # Set some attributes that will be set/used during the workflow
        self.iteration = 0
        self.step = 0
        self.synthetics_only = bool(par["CASE"].lower() == "synthetic")
        
    @property
    def model_number(self):
        """
        The model number is based on the current iteration
        """
        return f"m{max(self.iteration - 1, 0):0>2}"

    @property
    def step_count(self):
        """
        Step count based on
        """
        return f"s{self.step:0>2}"

    def set(self, **kwargs):
        """
        Convenience function to easily set multiple parameters before calling
        other functions.

        Overwrite internally used attributes using kwargs. Ensure that
        attributes other than the ones set in __init__ are allowed.
        """
        for key in list(kwargs.keys()):
            if not hasattr(self, key):
                warnings.warn(f"Pyaflowa has no attribute '{key}', ignoring")
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
        ev_paths = {"stations": os.path.join(cwd, "DATA", "STATIONS"),
                    "maps": os.path.join(self.int_paths["maps"], event_id),
                    "figures": os.path.join(self.int_paths["figures"],
                                            self.model_number, event_id),
                    }
        
        # Make the process specific event directories
        for key, item in ev_paths.items():
            if not os.path.exists(item):
                os.makedirs(item)

        # Set logging output for Pyflex and Pyatoa, less output using 'info'
        for log, level in self.par["set_logging"].items():
            if level:
                logger = logging.getLogger(log)
                if level == "info":
                    logger.setLevel(logging.INFO)
                elif level == "debug":
                    logger.setLevel(logging.DEBUG)

        # Read in the Pyatoa Config object and set attributes based on workflow
        config = pyatoa.Config(yaml_fid=self.int_paths["config_file"])
        setattr(config, "event_id", event_id)
        setattr(config, "model_number", self.model_number)
        setattr(config, "synthetic_tag", f"synthetic_{self.model_number}")
        setattr(config, "synthetics_only", self.synthetics_only)

        # Make sure Pyatoa knows to look in the Seisflows directories for data
        config.cfgpaths["synthetics"].append(os.path.join(cwd, "traces", "syn"))
        config.cfgpaths["waveforms"].append(os.path.join(cwd, "traces", "obs"))

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
        ds_name = os.path.join(self.int_paths["data"],
                               f"{config.event_id}.h5")

        # Count number of successful processes
        successes = 0
        with pyasdf.ASDFDataSet(ds_name) as ds:
            # Make sure the ASDFDataSet doesn't already contain auxiliary_data
            clean_ds(ds=ds, model=self.model_number, step=self.step_count,
                     fix_windows=self.par["fix_windows"])

            # Write the Config to auxiliary_data for provenance
            config.write(write_to=ds)
            mgmt = pyatoa.Manager(config=config, ds=ds)

            # Get stations from Specfem STATIONS file in form NET STA LAT LON ..
            stations = np.loadtxt(ev_paths["stations"], usecols=[0, 1, 2, 3],
                                  dtype=str)
            coords = stations[:, 2:]

            # Loop through stations and invoke Pyatoa workflow
            for station in stations:
                sta, net = station[:2]
                print(f"{net}.{sta}")
                try:
                    mgmt.reset()
                    mgmt.gather(station_code=f"{net}.{sta}.*.HH*")
                    mgmt.standardize()
                    mgmt.preprocess(overwrite=preproc)
                    mgmt.window(fix_windows=self.par["fix_windows"])
                    mgmt.measure()

                    # Plot waveforms with misfit windows and adjoint sources
                    if self.par["plot_waveforms"]:
                        # Format some strings to append to the waveform title
                        tit = " ".join([
                            f"\n{config.model_number}{self.step_count}",
                            f"pyflex={config.pyflex_map},",
                            f"pyadjoint={config.adj_src_type},",
                            f"misfit={mgmt.misfit:.2E}"
                        ])
                        mgmt.plot(append_title=tit,
                                  save=os.path.join(
                                      ev_paths["figures"], f"wav_{sta}"),
                                  show=False, return_figure=False
                                  )

                        # Plot source-receiver maps of waveforms were plotted
                        if self.par["plot_srcrcv_maps"]:
                            map_fid = os.path.join(ev_paths["maps"],
                                                   f"map_{sta}")
                            if not os.path.exists(map_fid):
                                mgmt.srcrcvmap(stations=coords, save=map_fid,
                                               show=False)
                    print("\n")
                    successes += 1
                # Traceback ensures more detailed error tracking
                except Exception:
                    traceback.print_exc()
                    print("\n")
                    continue

            # Run finalization procedures for processing if gathered waveforms
            if successes:
                self.finalize_process(ds=ds, cwd=cwd, ev_paths=ev_paths,
                                      config=config)
            else:
                print("Pyaflowa processed no data, skipping finalize")

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
        # Write adjoint sources directly to the Seisflows traces/adj dir
        print("exporting files to Specfem3D")
        print("\twriting adjoint sources to .sem? files...")
        write_adj_src_to_ascii(ds, config.model_number,
                               os.path.join(cwd, "traces", "adj"))

        # Write the STATIONS_ADJOINT file to the DATA directory of cwd
        print("\tcreating STATIONS_ADJOINT file...")
        create_stations_adjoint(ds, config.model_number,
                                specfem_station_file=ev_paths["stations"],
                                pathout=os.path.join(cwd, "DATA")
                                )

        print("exporting files to Seisflows")
        print("\twriting event misfit to disk...")
        write_misfit_stats(ds, config.model_number, self.int_paths["misfits"])

        print("writing files for internal use")
        print("\twriting stats to ASDF file...")
        write_stats_to_asdf(ds, config.model_number, self.step_count)

        print("\twriting misfits.json to disk...")
        write_misfit_json(ds, self.model_number, self.step_count,
                          self.int_paths["misfit_file"])

        # Only run this for the first 'step', otherwise we get too many pdfs
        if self.par["combine_imgs"] and (self.step_count == "s00"):
            print("\tcreating composite pdf...")

            # Create the name of the pdf to save to
            save_to = os.path.join(self.int_paths["composites"],
                                   f"{config.event_id}_{config.model_number}_"
                                   f"{self.step_count}_wavmap.pdf"
                                   )
            tile_combine_imgs(ds=ds, save_pdf_to=save_to,
                              wavs_path=ev_paths["figures"],
                              maps_path=ev_paths["maps"],
                              purge_wavs=self.par["purge_waveforms"],
                              purge_tiles=self.par["purge_tiles"]
                              )

    def finalize(self):
        """
        At the end of an iteration, clean up working directory and create final
        objects if requested by the User
        """
        # Plot the output.optim file outputted by Seisflows
        plot_output_optim(path_to_optim=os.path.join(self.ext_paths["WORKDIR"],
                                                     "output.optim"),
                          save=os.path.join(self.int_paths["figures"],
                                            "output_optim.png")
                          )

        # Generate .vtk files for given source and receivers
        if self.par["create_srcrcv_vtk"]:
            create_srcrcv_vtk_multiple(pathin=self.int_paths["data"],
                                       pathout=self.int_paths["vtks"],
                                       model=self.model_number
                                       )

        # Create copies of .h5 files at the end of each iteration, because .h5
        # files are easy to corrupt so it's good to have a backup
        if self.par["snapshot"]:
            srcs = glob.glob(os.path.join(self.int_paths["data"], "*.h5"))
            for src in srcs:
                shutil.copy(src, os.path.join(self.int_paths["snapshots"],
                                              os.path.basename(src))
                            )

