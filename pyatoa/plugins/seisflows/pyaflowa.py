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
from pyatoa.utils.asdf.extractions import windows_from_ds
from pyatoa.utils.visuals.statistics import plot_output_optim
from pyatoa.utils.visuals.mapping import event_misfit_map
from pyatoa.utils.visuals.plot_tools import imgs_to_pdf
from pyatoa.utils.tools.io import (create_stations_adjoint, write_misfit_json,
                                   write_adj_src_to_ascii, write_misfit_stats,
                                   tile_combine_imgs,
                                   create_srcrcv_vtk_multiple)


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
        # Ensure that these are accessible by the class
        self.par = par["PYATOA"]
        self.ext_paths = paths

        # Distribute internal hardcoded path structure
        assert("PYATOA_IO" in self.ext_paths.keys())
        pyatoa_io = self.ext_paths["PYATOA_IO"]
        self.int_paths = {
            "config_file": os.path.join(self.ext_paths["WORKDIR"],
                                        "parameters.yaml"),
            "misfit_file": os.path.join(pyatoa_io, "misfits.json"),
            "figures": os.path.join(pyatoa_io, "figures"),
            "pyatoa_data": os.path.join(pyatoa_io, "data"),
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

    def _setup_process(self, cwd, event_id=None):
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
        event = {"stations": os.path.join(cwd, "DATA", "STATIONS"),
                 "maps": os.path.join(self.int_paths["PYATOA_MAPS"], event_id)
                 "figures": os.path.join(self.int_paths["figures"],
                                         self.model_number, event_id),
                 }

        # Set logging output for Pyflex and Pyatoa, less output using 'info'
        if self.par["set_logging"]:
            logger_pyatoa = logging.getLogger("pyatoa")
            logger_pyflex = logging.getLogger("pyflex")
            if self.par["set_logging"] == "info":
                logger_pyatoa.setLevel(logging.INFO)
                logger_pyflex.setLevel(logging.INFO)
            else:
                logger_pyatoa.setLevel(logging.DEBUG)
                logger_pyflex.setLevel(logging.DEBUG)

        # Read in the Pyatoa Config object and set attributes based on workflow
        config = pyatoa.Config(yaml_fid=self.int_paths["config_file"])
        setattr(config, "event_id", event_id)
        setattr(config, "model_number", self.model_number)
        setattr(config, "synthetic_tag", f"synthetic_{self.model_number}")

        # Make sure Pyatoa knows to look in the Seisflows directories for data
        config.cfgpaths["synthetics"].append(os.path.join(cwd, "traces", "syn"))
        config.cfgpaths["waveforms"].append(os.path.join(cwd, "traces", "obs"))

        return config, event

    def process(self, cwd, event_id=None):
        """
        Main workflow calling on the core functionality of Pyatoa to process
        observed and synthetic waveforms and perform misfit quantification

        :type cwd: str
        :param cwd: current working directory for this instance of Pyatoa
        :type event_id: str
        :param event_id: event identifier tag for file naming etc.
        """
        config, event = self._setup_process(cwd, event_id)

        ds_name = os.path.join(self.int_paths["pyatoa_data"],
                               f"{config.event_id}.h5")
        with pyasdf.ASDFDataSet(ds_name) as ds:
            # Make sure the ASDFDataSet doesn't already contain auxiliary_data
            clean_ds(ds=ds, model=self.model_number, step=self.step_count,
                     fix_windows=self.par["fix_windows"])

            # Write the Config to auxiliary_data for provenance
            config.write(write_to=ds)

            # Instantiate the Manager
            mgmt = pyatoa.Manager(config=config, ds=ds)

            # Get stations from Specfem STATIONS file in form NET STA LAT LON ..
            stations = np.loadtxt(event["stations"], usecols=[0, 1, 2, 3],
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
                    if not self.pyatoa_par["fix_windows"] or \
                            not hasattr(ds.auxiliary_data.MisfitWindows,
                                        config.model_number):
                        mgmt.window()
                    else:
                        # If windows exist and fixed windows, grab from ASDF
                        misfit_windows = windows_from_ds(
                            ds, config.model_number, net, sta)
                        mgmt.windows = misfit_windows

                    mgmt.measure()

                    # Plot waveforms with misfit windows and adjoint sources
                    if self.pars["plot_waveforms"]:
                        # Format some strings to append to the waveform plot title
                        append_title = (
                            f"\n{config.model_number}{parser.step_count} "
                            f"pyflex={config.pyflex_map}, "
                            f"pyadjoint={config.adj_src_type},")
                        if mgmt.misfit is not None:
                            append_title = " ".join(
                                [append_title, f"misfit={mgmt.misfit:.2E}"])
                        f = mgmt.plot(
                            append_title=append_title,
                            save=os.path.join(event["figures"], f"wav_{sta}"),
                            show=False, return_figure=True
                        )

                    # Plot source-receiver maps, don't make a map if no wav data
                    # Don't make the map if the map has already been made
                    if self.par["plot_srcrcv_maps"] and f:
                        map_fid = os.path.join(event["maps"], f"map_{sta}")
                        if not os.path.exists(map_fid):
                            mgmt.srcrcvmap(stations=coords, save=map_fid,
                                           show=False)
                    print("\n")
                # Traceback ensures more detailed error tracking
                except Exception:
                    traceback.print_exc()
                    print("\n")
                    continue

            # Run finalization procedures for processing
            self._finalize_process(ds=ds, cwd=cwd, event=event, config=config)



    def _export_specfem3d(self, **kwargs):
        """
        Create necessary input files for the adjoint solver of Specfem3D
        Place them into the proper directories within Seisflows, so that
        Seisflows can find these input files
        """
        ds = kwargs.get("ds", None)
        cwd = kwargs.get("cwd", None)
        config = kwargs.get("config", None)

        # Write adjoint sources directly to the Seisflows traces/adj dir
        print("writing adjoint sources to .sem? files...")
        write_adj_src_to_ascii(ds, config.model_number,
                               os.path.join(cwd, "traces", "adj"))

        # Write the STATIONS_ADJOINT file to the DATA directory of cwd
        print("creating STATIONS_ADJOINT file...")
        create_stations_adjoint(ds, config.model_number,
                                specfem_station_file=self.int_paths["stations"],
                                pathout=os.path.join(cwd, "DATA")
                                )

    def _finalize_process(self, **kwargs):
        """
        After all waveforms have been windowed and measured, run some functions
        that create output files useful for Specfem, or for the User.

        Pass arguemnts as kwargs to give some flexibility to input parameters
        """
        ds = kwargs.get("ds", None)
        event = kwargs.get("event", None)
        config = kwargs.get("config", None)

        self._export_specfem3d(**kwargs)

        print("writing stats to ASDF file...")
        write_stats_to_asdf(ds, config.model_number, self.step_count)

        print("writing individual misfit to file...")
        write_misfit_stats(ds, config.model_number,
                           self.int_paths["misfits"])

        print("writing misfits.json file...")
        write_misfit_json(ds, self.model_number, self.step_count,
                          self.int_paths["misfit_file"])

        # Only run this for the first 'step', otherwise we get too many pdfs
        if self.pars["combine_imgs"] and (self.step_count == "s00"):
            print("creating composite pdf...")

            # Create the name of the pdf to save to
            save_to = os.path.join(self.int_paths["composites"],
                                   f"{config.event_id}_{config.model_number}_"
                                   f"{parser.step_count}_wavmap.pdf"
                                   )
            tile_combine_imgs(ds=ds, save_pdf_to=save_to,
                              wavs_path=event["figures"], maps_path=event["maps"],
                              purge_wavs=self.pars["purge_waveforms"],
                              purge_tiles=self.pars["purge_tiles"]
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

        # Create misfit maps for each event with contour overlay showing misfit
        if self.par["plot_misfit_maps"] and (self.step_count == "s00"):
            name_template = "{eid}_{m}_{s}_misfit_map.png"
            file_ids = []
            # Loop through each available dataset to create misfit map
            datasets = glob.glob(os.path.join(self.int_paths["data"], "*.h5"))
            for dataset in datasets:
                with pyasdf.ASDFDataSet(dataset) as ds:
                    event_id = os.path.basename(ds.filename).split('.')[0]

                    # Save figures into event directories
                    event_figures = os.path.join(self.int_paths["figures"],
                                                 self.model_number, event_id
                                                 )
                    # Save the fid based on event id, model number, step count
                    fidout = os.path.join(
                        event_figures, name_template.format(eid=event_id,
                                                            m=self.model_number,
                                                            s=self.step_count)
                    )
                    file_ids.append(fidout)

                    # Use the average misfit to normalize the misfit map
                    stats = ds.auxiliary_data.Statistics[self.model_number][
                        self.step_count].parameters
                    average_misfit = stats['average_misfit']

                    event_misfit_map(map_corners=self.pars["map_corners"],
                                     ds=ds, model=self.model_number,
                                     step=self.step_count,
                                     normalize=average_misfit,
                                     annotate_station_info='simple',
                                     contour_overlay=True, filled_contours=True,
                                     show=False, save=fidout
                                     )

            # Combine all the misfit maps into a single pdf
            save_to = os.path.join(
                self.int_paths["composites"],
                f"{self.model_number}_{self.step_count}_misfitmaps.pdf"
            )
            imgs_to_pdf(fids=file_ids, fid_out=save_to)

