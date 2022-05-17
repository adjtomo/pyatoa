#!/usr/bin/env python
"""
A SeisFlows3 plugin which simplifies how SeisFlows3 calls Pyatoa during an
active workflow. Pyaflowa is meant to be called from within the
SeisFlows3.preprocess.pyatoa module.

Pyaflowa is in charge of navigating through the SeisFlows3 directory structure,
collecting required data including waveforms and metadata, processing waveforms
to generate misfit windows and adjoint sources, store data, metadata and results
in ASDFDataSets, and generate figures.

To call Pyaflowa from inside an active SeisFlows3 working environment:
.. rubric::
    from pyatoa import Pyaflowa
    pyaflowa = Pyaflowa(sfpar=PAR, sfpath=PATH)  # PAR and PATH defined by SF3
    pyaflowa.setup(source_name=solver.source_name[0],
                   iteration=optimize.iter,
                   step_count=optimize.line_search.step_count,
                   loc="*", cha="*")
    misfit = pyaflowa.process()
"""
import os
import logging
import time
import random
import shutil
import numpy as np
import pyatoa
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from pyasdf import ASDFDataSet

from pyatoa.utils.images import imgs_to_pdf, merge_pdfs
from pyatoa.utils.read import read_station_codes
from pyatoa.utils.asdf.clean import clean_dataset


class Paths(dict):
    """Dictionary that prints keys like attributes for cleaner calls and allows
    on the fly formatting of all values within the dict. Assumes that all
    values are strings.
       
    .. note::
        Had some trouble pickling this class, fixed see: 
        https://stackoverflow.com/questions/42272335/
                         how-to-make-a-class-which-has-getattr-properly-pickable
    """
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)

    @staticmethod
    def _make_nonexistent_dir(key, val):
        """Checks if a path is a directory by whether it has a file ext, and
        whether or not the dir actually exists"""
        if "_file" in key:
            make_dir = False
        else:
            if os.path.splitext(val)[1]:
                make_dir = False
            else:
                if os.path.exists(val):
                    make_dir = False
                else:
                    make_dir = True
        if make_dir:
            os.makedirs(val)

    def format(self, **kwargs):
        """
        Format any open string formatters with any available kwargs. Allow
        for list-like values where all values in the list will be formatted. 
        Also allow for making whatever looks like a directory (no file ext) 
        to establish the path structure on disk
        """
        path_out = self
        for key, val in self.items():
            try:
                new_val = val.format(**kwargs)
                path_out[key] = new_val
                self._make_nonexistent_dir(key, new_val)
            # Lists will not format by themselves, need to go in and format
            # each entry separately.
            except AttributeError:
                new_vals = []
                for item in val:
                    new_val = item.format(**kwargs)
                    self._make_nonexistent_dir(key, new_val)
                    new_vals.append(new_val)
                path_out[key] = new_vals
        return path_out


class Pyaflowa:
    """
    A class that simplifies calling the Pyatoa waveform misfit quantification
    workflow inside an active SeisFlows3 workflow
    """
    def __init__(self, sfpar, sfpath):
        """
        Initialize Pyaflowa based on the SeisFlows3 PAR and PATH dictionaries,
        which define the parameter suite and required paths used in the
        workflow.

        :type sfpar: seisflows.config.Dict
        :param sfpar: Parameter definitions defined by var `PAR` in SeisFlows3
        :type sfpath: seisflows.config.Dict
        :param sfpath: Path definitions defined by variable `PATH` in SeisFlows3
        """
        self.sfpar = sfpar
        self.sfpath = sfpath

        # Attributes to be created and formatted by setup()
        self.paths = self.define_path_structure(self.sfpath)  # Unformatted
        self.config = None
        self.logger = None
        self.codes = None
        self.nproc = os.cpu_count()  # default behavior: fill out the node

        # Internal stats attributes used by process_single_station()
        self.fix_windows = False

    @staticmethod
    def define_path_structure(sfpath):
        """
        Define the internal workflow directory structure based on the SeisFlows3
        PATH dictionary. Allows Pyaflowa to know how to navigate around and
        where to find data and save results.

        .. note::
            Path definitions are as follows:
            * cwd: The individual event working directory. This is modelled
              after a SPECFEM working directory, and for SeisFlows points
              directly to each of the solver working directories
            * datasets: The path where Pyaflowa is allowed to read and write
              ASDFDataSets which contain all the working data.
            * figures: Path where Pyaflowa is allowed to save the PDFS/PNGS
              that result from the workflow
            * tmplogs: Path where temporary log files are stored, which are
              later collected to the main log file and then deleted
            * logs: Path where Pyaflowa can store the log outputs from each
              of the individual event workflows
            * responses: The path to response files stored on disk in SEED fmt
              that the Gatherer object will search in order to find
              corresponding StationXML files. Can be left blank if StationXML
              files will be queried from FDSN
            * waveforms: Path to waveform files stored on disk in SEED format,
              same caveat as responses
            * synthetics: Paths to synthetic waveforms generated by SPECFEM.
              Need to be stored in directories pertaining to source names
              e.g. path/to/synthetics/SOURCE_NAME/*semd
            * stations_file: Path to the STATIONS file that dictates the
              receivers being used in the processing step. This file needs
              to match the SPECFEM3D format, and be saved into the source dir.
              e.g. path/to/SOURCE_NAME/STATIONS
            * adjsrcs: Path to save the resulting adjoint sources that will be
              generated during the processing workflow.

        :type sfpath: seisflows.config.Dict
        :param sfpath: Path definitions defined by variable `PATH` in SeisFlows3
        """
        # DATA path may be None, e.g., during synthetic inversion
        if sfpath.DATA is None:
            path_data = os.path.join(sfpath.PREPROCESS, "data")
        else:
            path_data = sfpath.DATA

        paths = Paths(
            workdir=sfpath.WORKDIR,
            cwd=os.path.join(sfpath.SOLVER, "{source_name}"),
            data=os.path.join(sfpath.SOLVER, "{source_name}", "DATA"),
            datasets=os.path.join(sfpath.PREPROCESS, "datasets"),
            figures=os.path.join(sfpath.PREPROCESS, "figures"),
            logs=os.path.join(sfpath.PREPROCESS, "logs"),
            tmplogs=os.path.join(sfpath.PREPROCESS, "logs", "{source_name}"),
            ds_file=os.path.join(sfpath.PREPROCESS, "datasets",
                                 "{source_name}.h5"),
            event_figures=os.path.join(sfpath.PREPROCESS, "figures",
                                       "{source_name}"),
            responses=os.path.join(path_data, "seed"),
            stations_file=os.path.join(sfpath.SOLVER, "{source_name}", "DATA",
                                       "STATIONS"),
            synthetics=os.path.join(sfpath.SOLVER, "{source_name}",
                                    "traces", "syn"),
            adjsrcs=os.path.join(sfpath.SOLVER, "{source_name}", "traces",
                                 "adj"),
            waveforms=[os.path.join(path_data, "mseed"),
                       os.path.join(sfpath.SOLVER, "{source_name}", "traces",
                                    "obs")]
        )
        return paths

    def setup(self, source_name, iteration, step_count, loc="*", cha="*"):
        """
        Format the paths and config object with the current place in the
        workflow. This must be called each time we use Pyaflowa, as there are
        source and evaluation-specific filenames and paths that need to be
        generated each time.

        :type source_name: str
        :param source_name: the name of the event which is used to tag the
            solver working directory. In SeisFlows3 this is usually defined
            by `solver.source_names`
        :type iteration: int
        :param iteration: The current iteration of the SeisFlows3 workflow,
            within SeisFlows3 this is defined by `optimize.iter`
        :type step_count: int
        :param step_count: Current line search step count within the SeisFlows3
            workflow. Within SeisFlows3 this is defined by
            `optimize.line_search.step_count`
        :type loc: str
        :param loc: How to format location codes as in NET.STA.LOC.CHA for data
            fetching. Since SPECFEM STATIONS files do not specify locations,
            we must do it ourselves. Defaults to a wildcard '*'.
        :param cha: How to format channel codes as in NET.STA.LOC.CHA for data
            fetching. Since SPECFEM STATIONS files do not specify channel,
            we must do it ourselves. Defaults to a wildcard '*'.
        """
        # Set path structure, formatted by the specific event tag
        unformatted_paths = self.define_path_structure(self.sfpath)
        self.paths = unformatted_paths.format(source_name=source_name)

        # Generate the configuration parameter which will control processing
        self.config = self.generate_config(source_name, iteration, step_count)

        # Tag and create a unique event log file which will track processing.
        # Write a nice log header to make things easy to identify
        self.logger = self.generate_logger(self.config)
        _header = f"EVENT: {source_name} // EVAL: {self.config.eval_tag}"
        self.logger.info(f"\n{'=' * 80}\n{_header:^80}\n{'=' * 80}\n")

        # Codes define the stations used in this workflow
        self.codes = read_station_codes(self.paths.stations_file, loc=loc,
                                        cha=cha)

        # Figure out how to address whether we re-use misfit windows
        self.fix_windows = self.check_fix_windows(iteration, step_count)

        with ASDFDataSet(self.paths.ds_file) as ds:
            # Make sure Dataset does not have previous processing information
            clean_dataset(ds, iteration=iteration, step_count=step_count)
            self.config.write(write_to=ds)

            # Gather event from solver/*/DATA/ directory and save to ASDFDataSet
            mgmt = pyatoa.Manager(ds=ds, config=self.config)
            mgmt.gather(choice=["event"], event_id="",
                        prefix=self.sfpar.SOURCE_PREFIX)

    def generate_config(self, source_name, iteration, step_count):
        """
        Generate a Pyatoa.Config object based on the current location in the
        workflow and on the pre-defined path structure.

        :type source_name: str
        :param source_name: the name of the event which is used to tag the
            solver working directory. In SeisFlows3 this is usually defined
            by `solver.source_names`
        :type iteration: int
        :param iteration: The current iteration of the SeisFlows3 workflow,
            within SeisFlows3 this is defined by `optimize.iter`
        :type step_count: int
        :param step_count: Current line search step count within the SeisFlows3
            workflow. Within SeisFlows3 this is defined by
            `optimize.line_search.step_count`
        """
        config = pyatoa.Config(seisflows_par=self.sfpar, iteration=iteration,
                               step_count=step_count, event_id=source_name,
                               paths={"responses": self.paths.responses,
                                      "waveforms": self.paths.waveforms,
                                      "synthetics": self.paths.synthetics,
                                      "events": self.paths.data},
                               )
        # Only query FDSN at the very first function evaluation
        if iteration != 1 and step_count != 0:
            config.client = None

        return config

    def generate_logger(self, config, prefix="", tmpdir=False):
        """
        Create a log file to spit out any log statements that occur during
        processing. Removes handlers that are by default defined in the Pyatoa
        __init__ script so that log statements don't end up printing to stdout
        during a SeisFlows3 workflow.

        :type config: pyatoa.core.config.Config
        :param config: Pyatoa Config object from which we take the iteration,
            step count and event name
        :type prefix: str
        :param prefix: optional prefix for logs, to create process-specific
            log files (i.e., for every source-receiver pair)
        :type tmpdir: bool
        :param tmpdir: if True, redirect log files to a temporary directory
            stored in the `path` attribute, to be collected later. This is for
            station specific logging
        :rtype: logging.Logger
        :return: an individualized logging handler
        """
        # Generate the name of the config based on our place in the workflow
        if config.eval_tag is not None:
            log_fid = f"{prefix}{config.eval_tag}_{config.event_id}.txt"
        else:
            log_fid = f"{prefix}{config.event_id}.txt"

        if tmpdir:
            log_path = os.path.join(self.paths.tmplogs, log_fid)
        else:
            log_path = os.path.join(self.paths.logs, log_fid)

        # Convenience to not carry around long attribute call
        log_level = self.sfpar.PYATOA_LOG_LEVEL.upper()

        # Propagate logging to individual log files, always overwrite
        handler = logging.FileHandler(log_path, mode="w")

        # Maintain the same look as the standard console log messages
        logfmt = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
        formatter = logging.Formatter(logfmt, datefmt="%Y-%m-%d %H:%M:%S")
        handler.setFormatter(formatter)
        handler.setLevel(log_level)

        for log in ["pyflex", "pyadjoint", "pyatoa"]:
            # Set the overall log level
            logger = logging.getLogger(log)
            # Turn off any existing handlers (stream and file)
            while logger.hasHandlers():
                logger.removeHandler(logger.handlers[0])

            logger.setLevel(log_level)
            logger.addHandler(handler)

        return logger

    def check_fix_windows(self, iteration, step_count):
        """
        Determine how to address re-using misfit windows during an inversion
        workflow. Throw some log messages out to let the User know whether or
        not misfit windows will be re used throughout an inversion.

        Options for SeisFlows3 PAR.FIX_WINDOWS:
            True: Always fix windows except for i01s00 because we don't have any
                  windows for the first function evaluation
            False: Don't fix windows, always choose a new set of windows
            Iter: Pick windows only on the initial step count (0th) for each
                  iteration. WARNING - does not work well with Thrifty Inversion
                  because the 0th step count is usually skipped
            Once: Pick new windows on the first function evaluation and then fix
                  windows. Useful for when parameters have changed, e.g. filter
                  bounds

        :type iteration: int
        :param iteration: The current iteration of the SeisFlows3 workflow,
            within SeisFlows3 this is defined by `optimize.iter`
        :type step_count: int
        :param step_count: Current line search step count within the SeisFlows3
            workflow. Within SeisFlows3 this is defined by
            `optimize.line_search.step_count`
        :rtype: bool
        :return: bool on whether to use windows from the previous step
        """
        fix_windows = False
        # First function evaluation never fixes windows
        if iteration == 1 and step_count == 0:
            fix_windows = False
            self.logger.info("new windows; first evaluation")
        elif isinstance(self.sfpar.FIX_WINDOWS, str):
            # By 'iter'ation only pick new windows on the first step count
            if self.sfpar.FIX_WINDOWS.upper() == "ITER":
                if step_count == 0:
                    fix_windows = False
                    self.logger.info("new windows; first step count")
                else:
                    fix_windows = True
                    self.logger.info("fix windows; mid line search")
            # 'Once' picks windows only for the first function evaluation of
            # the current set of iterations.
            elif self.sfpar.FIX_WINDOWS.upper() == "ONCE":
                if iteration == self.sfpar.BEGIN and step_count == 0:
                    fix_windows = False
                    self.logger.info("new windows; first workflow evaluation")
                else:
                    fix_windows = True
                    self.logger.info("fix windows; mid workflow")
        # Bool fix windows simply sets the parameter
        elif isinstance(self.sfpar.FIX_WINDOWS, bool):
            fix_windows = self.sfpar.FIX_WINDOWS
            self.logger.info(f"fixed windows flag set: "
                             f"{self.sfpar.FIX_WINDOWS}")

        return fix_windows

    def process(self, nproc=None):
        """
        Calculate misfit given a set of data-synthetic waveforms for a single
        event and a list of stations. Save all of the output into ASDFDataSets.
        Generate the required STATIONS_ADJOINT file and adjoint source files.

        .. note::
            For serial processing set `nproc`==1

        :type nproc: int
        :param nproc: number of processors to use concurrently for parallel
            processing. If not given, defaults to internal attribute nproc
            which will try to fill out the entire CPU
        """
        _start = time.time()
        mgmt_stats = {}
        # Use concurrent futures to do station processing in parallel
        with ProcessPoolExecutor(max_workers=nproc or self.nproc) as executor:
            for code, stats in zip(self.codes,
                                    executor.map(self._process_single_station,
                                                 self.codes)):
                mgmt_stats[code] = stats

        # Calculate total misfit based on number of windows collected
        misfit, nwin = 0, 0
        for sta, stats in mgmt_stats.items():
            if stats.status_ok:
                misfit += stats.misfit
                nwin += stats.nwin or 0  # deal with NoneType nwin values
        scaled_misfit = self.calculate_misfit(misfit=misfit, nwin=nwin)

        # Remainder of these tasks are run in serial; collect processing stats
        self.make_event_figure_pdf()
        self.write_stations_adjoint()

        # Finalization log statement for the user
        _end = time.time()
        status_str = self._parse_stats_for_log(mgmt_stats)
        msg = (f"\n{'=' * 80}\n{'SUMMARY':^80}\n{'=' * 80}\n"
               f"NPROC: {nproc} // TIME: {_end - _start:.2f}s\n"
               f"SOURCE NAME: {self.config.event_id}\n"
               f"WINDOWS: {nwin}\n"
               f"RAW MISFIT: {misfit:.4f}\n"
               f"{status_str}"
               f"\n{'=' * 80}\n{'RAW LOGS':^80}\n{'=' * 80}"
               )
        self.logger.info(msg)
        self.collect_tmp_log_files()

        return scaled_misfit

    def collect_tmp_log_files(self):
        """
        Each source-receiver pair has made its own log file. This function
        collects these files and writes their content back into the main log.
        This is a lot of IO but should be okay since the files are small.

        .. note::
            This was the most foolproof method for having multiple parallel
            processes write to the same file. I played around with StringIO
            buffers and file locks, but they became overly complicated and
            ultimately did not work how I wanted them to. This function trades
            filecount and IO overhead for simplicity.

        .. warning::
            The assumption here is that the number of source-receiver pairs
            is manageable (in the thousands). If we start reaching file count
            limits on the cluster then this method for logging may have to be
            re-thought. See link for example:

            https://stackless.readthedocs.io/en/3.7-slp/howto/
              logging-cookbook.html#using-concurrent-futures-processpoolexecutor
        """
        tmp_logs = sorted(glob(os.path.join(self.paths.tmplogs, "*")))
        with open(self.logger.handlers[0].baseFilename, "a") as fw:
            for tmp_log in tmp_logs:
                with open(tmp_log, "r") as fr:
                    fw.writelines(fr.readlines())
        shutil.rmtree(self.paths.tmplogs)

    def _process_single_station(self, code):
        """
        Process a single seismic station for a given event. Multiple error
        catching chunks to ensure that a failed processing for a single station
        won't kill the entire job.

        .. note::
            Each station gets tagged with a status that let's the user know how
            processing went. Manager status returns are:
            * 'STATION_OK': Processing completed with expected behavior
            * 'STATION_NO_DATA': Processing error; failed to return some data
            * 'STATION_CONTROLLED_ERROR': Processing error but within expected
                failure criteria for processing (e.g., large misfit)
            * 'STATION_UNCONTROLLED_ERROR': Processing error; unexpected failure
                criteria should probably be investigated by user

        :type code: str
        :param code: Pyatoa station code (formatted NN.SSS.LL.CCC) used to
            select which station we are processing
        :rtype: pyatoa.core.manager.ManagerStats
        :return: manager stats dictionary which contains information about
            the misfit, number of windows, and processing status
        """
        net, sta, loc, cha = code.split(".")
        # Station-specific temporary log file which will be re-added to main
        logger = self.generate_logger(config=self.config,
                                      prefix=f"{net}_{sta}_", tmpdir=True)

        logger.info(f"\n{'/' * 80}\n{code:^80}\n{'/' * 80}")
        mgmt = pyatoa.Manager(config=self.config)
        # Data gathering chunk; if fail, do not continue
        try:
            mgmt.gather(code=code,  event_id="",
                        prefix=self.sfpar.SOURCE_PREFIX)
        except pyatoa.ManagerError as e:
            logger.warning(e)
            mgmt.stats.status_ok = "STATION_NO_DATA"
            return mgmt.stats
        # Data processing chunk; if fail, continue to plotting
        try:
            mgmt.flow(fix_windows=self.fix_windows)
            mgmt.stats.status_ok = "STATION_OK"
        except pyatoa.ManagerError as e:
            logger.warning(e)
            mgmt.stats.status_ok = "STATION_CONTROLLED_ERROR"
        except Exception as e:
            # Uncontrolled exceptions should be noted in more detail
            logger.warning(e, exc_info=True)
            mgmt.stats.status_ok = "STATION_UNCONTROLLED_ERROR"

        # Data plotting chunk; if fail continue to writing
        if self.sfpar.PLOT:
            # fid is e.g. path/i01s00_NZ_BFZ.png
            plot_fid = "_".join([mgmt.config.iter_tag,  # e.g., i01
                                 mgmt.config.step_tag,  # e.g., s00
                                 net,  # e.g., NZ
                                 sta + ".png"  # e.g., BFZ.pdf
                                 ])
            save = os.path.join(self.paths.event_figures, plot_fid)
            # Mapping may fail if no Inventory or Event object is attached
            # Waveform figures are expected to work if we have gotten this far
            try:
                mgmt.plot(choice="both", show=False, save=save)
            except pyatoa.ManagerError as e:
                logger.warning("cannot plot map, plotting waveform only")
                mgmt.plot(choice="wav", show=False, save=save)
            logger.info(f"plotting figure and saving to: {save}")

        # Final cleanup steps only if processing is successful
        if mgmt.stats.status_ok:
            # SPECFEM wants adjsrcs for each comp, regardless if it has data
            mgmt.write_adjsrcs(path=self.paths.adjsrcs, write_blanks=True)

        # Wait until the very end to write to the HDF5 file, then do it
        # pseudo-serially to get around trying to parallel write to HDF5 file
        while True:
            try:
                with ASDFDataSet(self.paths.ds_file) as ds:
                    mgmt.write(ds=ds)
                break
            except BlockingIOError:
                # Random sleep time [0,1]s to decrease chances of two processes
                # attempting to access at exactly the same time
                time.sleep(random.random())

        return mgmt.stats

    @staticmethod
    def calculate_misfit(misfit, nwin):
        """
        Calculate the misfit based on the number of windows. Equation from
        Tape et al. (2010). If no windows, misfit is simply raw misfit

        :type misfit: float
        :param misfit: total RAW event misfit for all data-synthetic pairs
        :type nwin: int
        :param nwin: number of misfit windows collected for all data-synthetic
            pairs for a given event
        """
        # Scale the raw event misfit by
        try:
            scaled_misfit = 0.5 * misfit / nwin
        # Dealing with the cases where 1) nwin==0 (signifying either no windows
        # found, or calc'ing misfit on whole trace) and 2) when misfit and nwin
        # are None (no misfit found for event)
        except (TypeError, ZeroDivisionError):
            scaled_misfit = misfit
        return scaled_misfit

    def _parse_stats_for_log(self, mgmt_stats):
        """
        Count and list all the station status codes to let the user know
        how each stations was processed and what status category it falls under

        :type mgmt_stats: pyatoa.core.manager.ManagerStats
        :param mgmt_stats: manager stats dictionary which contains information
            processing status
        :rtype: str
        :return: string with all statuses, their respective counts,
            and a list of stations which falls under each status. To be
            appended to the final log message
        """
        # Count the different output states and related codes to relay to user
        status_str = ""
        status_list = np.array([s.status_ok for s in mgmt_stats.values()])
        code_list = np.array([*mgmt_stats])  # list of station codes

        statuses, counts = np.unique(status_list, return_counts=True)

        for i, status in enumerate(statuses):
            status_str += (f"{status}: {counts[i]} / {len(status_list)}\n"
                           f"{code_list[np.where(status_list == status)[0]]}"
                           )
        return status_str

    def write_stations_adjoint(self):
        """
        Create the STATIONS_ADJOINT file required by SPECFEM to run an adjoint
        simulation. Should be run after all processing has occurred. Works by
        checking what stations have adjoint sources available and re-writing the
        existing STATIONS file that is on hand.
        """
        self.logger.info("writing SPECFEM file: STATIONS_ADJOINT")
        # These paths follow the structure of SeisFlows and SPECFEM
        adjoint_traces = glob(os.path.join(self.paths.adjsrcs, "*.adj"))
        # Simply append to end of file name e.g. "path/to/STATIONS" + "_ADJOINT"
        stations_adjoint = self.paths.stations_file + "_ADJOINT"
        # Determine the network and station names for each adjoint source
        # Note this will contain redundant values but thats okay because we're
        # only going to check membership using it, e.g. [['NZ', 'BFZ']]
        adjoint_stations = [os.path.basename(_).split(".")[:2] for _ in
                            adjoint_traces]
        # The STATION file is already formatted so leverage for STATIONS_ADJOINT
        lines_in = open(self.paths.stations_file, "r").readlines()
        with open(stations_adjoint, "w") as f_out:
            for line in lines_in:
                # Station file line format goes: STA NET LAT LON DEPTH BURIAL
                # and we only need: NET STA
                check = line.split()[:2][::-1]
                if check in adjoint_stations:
                    f_out.write(line)

    def make_event_figure_pdf(self, delete_originals=True):
        """
        Combine a list of single source-receiver PNGs into a single PDF file
        for the given event. Mostly a convenience function to make it easier
        to ingest waveform figures during a workflow.

        :type delete_originals: bool
        :param delete_originals: delete original pdf files after mergin
        """
        # e.g., i01s00_NZ_BFZ.png
        input_fids = glob(os.path.join(self.paths.event_figures, "*.png"))
        if not input_fids:
            self.logger.warning(f"no event figures found searching glob for: "
                                f"{self.paths.event_figures}/*.png")
            return
        # e.g. i01s00_2018p130600.pdf
        output_fid = (f"{self.config.iter_tag}{self.config.step_tag}_" 
                      f"{self.config.event_id}.pdf")
        self.logger.info(f"merge {len(input_fids)} figures into "
                         f"composite: '{output_fid}'")

        # Merge all output pdfs into a single pdf, delete originals
        save = os.path.join(self.paths.event_figures, output_fid)
        imgs_to_pdf(fids=sorted(input_fids), fid_out=save)
        if delete_originals:
            for fid in input_fids:
                os.remove(fid)

    def make_evaluation_composite_pdf(self, delete_originals=True):
        """
        Utility function to combine all PDFs generated by 
        make_event_figure_pdf() into a single PDF tagged by the current 
        evaluation (iteration, step count). This is meant to make it easier
        for the User to scroll through figures. Option to delete the original
        event-specific PDFs which are now redundant

        .. note::
            This can be run without running Pyaflowa.format()

        :type delete_originals: bool
        :param delete_originals: delete original pdf files after mergin
        """
        event_figures = glob(os.path.join(self.paths.figures, "*", "*.pdf"))
        if not event_figures:
            self.logger.warning(f"no composite PDFs found searching glob for: "
                                f"{self.paths.figures}/*/*.pdf")
            return
        # Collecting evaluation tags, e.g., ['i01s00', 'i01s01']
        tags = set([os.path.basename(_).split("_")[0] for _ in event_figures])
        for tag in tags:
            fids = [fid for fid in event_figures if tag in fid]
            fid_out = os.path.join(self.paths.figures, f"{tag}.pdf")

            merge_pdfs(fids=sorted(fids), fid_out=fid_out)
            if delete_originals:
                for fid in fids:
                    os.remove(fid)
