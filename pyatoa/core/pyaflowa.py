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
import pyatoa
import logging
from glob import glob
from pyasdf import ASDFDataSet

from pyatoa.utils.images import merge_pdfs
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

    def format(self, mkdir=True, **kwargs):
        """
        Format any open string formatters with any available kwargs. Allow
        for list-like values where all values in the list will be formatted. 
        Also allow for making whatever looks like a directory (no file ext) 
        to establish the path structure on disk

        :type mkdir: bool
        :param mkdir: make any directory-like value (no file ext) after 
            formatting 
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

        # Internal stats attributes used by process_single_station()
        self.mgmt = None
        self.misfit = None
        self.nwin = None
        self.stations = 0
        self.processed = 0
        self.exceptions = 0
        self.plot_fids = []
        self.fix_windows = False
        self.scaled_misfit = 0

    def reset(self):
        """
        Although it is better practice to initiate a NEW Pyaflowa instance each
        time processing is required, the reset function also allows clearing out
        any saved values from a previous run.
        """
        self.config = None
        self.logger = None
        self.codes = None
        self.mgmt = None
        self.misfit = None
        self.nwin = None
        self.stations = 0
        self.processed = 0
        self.exceptions = 0
        self.plot_fids = []
        self.fix_windows = False
        self.scaled_misfit = 0

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

        # Tag and create a unique event log file which will track processing
        log_name = self._create_logger_name(self.config)
        self.logger = self.generate_logger(log_name)

        # Codes define the stations used in this workflow
        self.codes = read_station_codes(self.paths.stations_file, loc=loc,
                                        cha=cha)

        # Figure out how to address whether or not we re-use misfit windows
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

    def _create_logger_name(self, config):
        """
        Generate the name of the log file which is dependent on the current
        evaluation as well as the event name

        :type config: pyatoa.core.config.Config
        :param config: Pyatoa Config object from which we take the iteration,
            step count and event name
        :rtype: str
        :return: full path and name of logging file
        """
        # If we have iteration or step_count set, use to tag the log file
        if config.eval_tag is not None:
            log_fid = f"{config.eval_tag}_{config.event_id}.txt"
        else:
            log_fid = f"{config.event_id}.txt"
        log_path = os.path.join(self.paths.logs, log_fid)
        return log_path

    def generate_logger(self, log_path):
        """
        Create a log file to spit out any log statements that occur during
        processing. Removes handlers that are by default defined in the Pyatoa
        __init__ script so that log statements don't end up printing to stdout
        during a SeisFlows3 workflow.

        :type log_path: str
        :param log_path: full path and name of the outputted log file
        :rtype: logging.Logger
        :return: an individualized logging handler
        """
        # Propagate logging to individual log files, always overwrite
        handler = logging.FileHandler(log_path, mode="w")

        # Maintain the same look as the standard console log messages
        logfmt = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
        formatter = logging.Formatter(logfmt, datefmt="%Y-%m-%d %H:%M:%S")
        handler.setFormatter(formatter)
        handler.setLevel(self.sfpar.PYATOA_LOG_LEVEL.upper())

        for log in ["pyflex", "pyadjoint", "pyatoa"]:
            # Set the overall log level
            logger = logging.getLogger(log)
            # Turn off any existing handlers (stream and file)
            while logger.hasHandlers():
                logger.removeHandler(logger.handlers[0])

            logger.setLevel(self.sfpar.PYATOA_LOG_LEVEL.upper())
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
            self.logger.info("first evaluation, gathering new misfit windows")
        elif isinstance(self.sfpar.FIX_WINDOWS, str):
            # By 'iter'ation only pick new windows on the first step count
            if self.sfpar.FIX_WINDOWS.upper() == "ITER":
                if step_count == 0:
                    fix_windows = False
                    self.logger.info("first step count of iteration, gathering "
                                     "new misfit windows")
                else:
                    fix_windows = True
                    self.logger.info("mid line search, re-using misfit windows")
            # 'Once' picks windows only for the first function evaluation of
            # the current set of iterations.
            elif self.sfpar.FIX_WINDOWS.upper() == "ONCE":
                if iteration == self.sfpar.BEGIN and step_count == 0:
                    fix_windows = False
                    self.logger.info("first evaluation, gathering new misfit "
                                     "windows")
                else:
                    fix_windows = True
                    self.logger.info("mid workflow, re-using misfit windows")
        # Bool fix windows simply sets the parameter
        elif isinstance(self.sfpar.FIX_WINDOWS, bool):
            fix_windows = self.sfpar.FIX_WINDOWS
            self.logger.info(f"option for re-using misfit windows is set to: "
                             f"{self.sfpar.FIX_WINDOWS}")

        return fix_windows

    def process(self):
        """
        The main processing function for Pyaflowa misfit quantification.

        Processes waveform data for all stations related to a given event,
        produces waveform and map plots during the processing step, saves data
        to an ASDFDataSet and writes adjoint sources and STATIONS_ADJOINT file,
        required by SPECFEM3D's adjoint simulations, to disk.

        :rtype: float or None
        :return: the total scaled misfit collected during the processing chain,
            scaled_misfit will return None if no windows have been found or
            no misfit was calculated
        """
        # Open the dataset as a context manager and process all events in serial
        with ASDFDataSet(self.paths.ds_file) as ds:
            self.mgmt = pyatoa.Manager(ds=ds, config=self.config)
            for code in self.codes:
                self.process_single_station(code=code)
        if self.plot_fids:
            self.make_event_figure_pdf()
        self.write_stations_adjoint()
        self.scaled_misfit = self.calculate_misfit()
        self.log_stats()
        return self.scaled_misfit

    def process_single_station(self, code):
        """
        Process a single seismic station for a given event. Multiple error
        catching chunks to ensure that a failed processing for a single station
        won't kill the entire job.

        :type code: str
        :param code: Pyatoa station code (formatted NN.SSS.LL.CCC) used to
            select which station we are processing
        """
        self.logger.info(f"\n{'=' * 80}\n\n{code}\n\n{'=' * 80}")
        self.stations += 1
        self.mgmt.reset()
        # Data gathering chunk; if fail, do not continue
        try:
            self.mgmt.gather(code=code)
        except pyatoa.ManagerError as e:
            self.logger.warning(e)
            return None
        # Data processing
        status_ok = self._manager_flow()
        # Data plotting
        if self.sfpar.PLOT:
            self._plot_waveform_and_map(code)
        # Final cleanup steps only if processing is successful
        if status_ok:
            self._record_stats_after_processing()
            # SPECFEM wants adjsrcs for each comp, regardless if it has data
            self.mgmt.write_adjsrcs(path=self.paths.adjsrcs,
                                    write_blanks=True)

    def _manager_flow(self):
        """
        Attempt to process data for the given Pyaflowa state. Allow unctronlled
        exceptions through so as to not break the workflow. Return a status
        to let the main function know if things should proceed.

        :rtype: bool
        :return: The status of the processing step. True means processing
            completed nominally (without error), False means processing failed.
        """
        # Data processing chunk; if fail, continue to plotting
        try:
            self.mgmt.flow(fix_windows=self.fix_windows)
            status_ok = True
        except pyatoa.ManagerError as e:
            self.logger.warning(e)
            status_ok = False
        except Exception as e:
            # Uncontrolled exceptions should be noted in more detail
            self.logger.warning(e, exc_info=True)
            self.exceptions += 1
            status_ok = False
        return status_ok

    def _plot_waveform_and_map(self, code):
        """
        Attempt to plot the internal Manager as a waveform + map figure with
        a unique tag that defines the evaluation number, event id etc.
        If Mapping fails (e.g., when we have no inventory or event object)
        try to fall back to waveform plotting only

        :type code: str
        :param code: Pyatoa station code, NN.SSS.LL.CCC
        """
        net, sta, loc, cha = code.split(".")
        # fid is e.g. path/i01s00_NZ_BFZ.pdf
        plot_fid = "_".join([self.mgmt.config.iter_tag,  # e.g., i01
                             self.mgmt.config.step_tag,  # e.g., s00
                             net,  # e.g., NZ
                             sta + ".pdf"  # e.g., BFZ.pdf
                             ])
        save = os.path.join(self.paths.event_figures, plot_fid)
        # Mapping may fail if no Inventory or Event object is attached
        # Waveform figures are expected to work if we have gotten this far
        try:
            self.mgmt.plot(choice="both", show=False, save=save)
        except pyatoa.ManagerError as e:
            self.logger.warning("cannot plot map, plotting waveform only")
            self.mgmt.plot(choice="wav", show=False, save=save)
        # If a plot is made, keep track so it can be merged later on
        self.logger.info(f"plotting figure and saving to: {save}")
        self.plot_fids.append(save)

    def _record_stats_after_processing(self):
        """
        Keep track of processing stats for each source-receiver pair. Thsi will
        be used to calculate the scaled misfit later

        .. note::
            We allow misfit and nwin to be None to signify that no misfits
            have been calculated. But if we calculate something then we
            need to overwrite the None values. Should only happen once.
        """
        if self.misfit is None:
            self.misfit = 0
        if self.nwin is None:
            self.nwin = 0
        # Raw misfit which will be scaled by the number of misfit windows
        self.misfit += self.mgmt.stats.misfit
        # The deal with the case where window selection is skipped and
        # mgmt.stats.nwin == None
        num_new_win = self.mgmt.stats.nwin or 0
        self.nwin += num_new_win
        self.processed += 1

    def log_stats(self):
        """
        A final log message after all the processing is done to wrap up log file
        """
        # Finalization log statement for the user
        msg = (f"\n{'=' * 80}\n\nSUMMARY\n\n{'=' * 80}\n"
               f"SOURCE NAME: {self.config.event_id}\n"
               f"STATIONS: {self.processed} / {self.stations}\n"
               f"WINDOWS: {self.nwin}\n"
               f"RAW MISFIT: {self.misfit}\n"
               f"UNEXPECTED ERRORS: {self.exceptions}"
               )
        self.logger.info(msg)

    def calculate_misfit(self):
        """
        Calculate the misfit based on the number of windows. Equation from
        Tape et al. (2010).
        """
        # Scale the raw event misfit by
        try:
            scaled_misfit = 0.5 * self.misfit / self.nwin
        # Dealing with the cases where 1) nwin==0 (signifying either no windows
        # found, or calc'ing misfit on whole trace) and 2) when misfit and nwin
        # are None (no misfit found for event)
        except (TypeError, ZeroDivisionError):
            scaled_misfit = self.misfit
        return scaled_misfit

    def write_stations_adjoint(self):
        """
        Create the STATIONS_ADJOINT file required by SPECFEM to run an adjoint
        simulation. Should be run after all processing has occurred. Works by
        checking what stations have adjoint sources available and re-writing the
        existing STATIONS file that is on hand.
        """
        self.logger.info("generating STATIONS_ADJOINT file for SPECFEM")
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

    def make_event_figure_pdf(self):
        """
        Combine a list of single source-receiver PDFS into a single PDF file
        for the given event. Mostly a convenience function to make it easier
        to ingest waveform figures during a workflow.
        """
        # e.g. i01s00_2018p130600.pdf
        output_fid = (f"{self.config.iter_tag}{self.config.step_tag}_" 
                      f"{self.config.event_id}.pdf")
        self.logger.info(f"merging {len(self.plot_fids)} output figures into a "
                         f"single pdf: '{output_fid}'")

        # Merge all output pdfs into a single pdf, delete originals
        save = os.path.join(self.paths.event_figures, output_fid)
        merge_pdfs(fids=sorted(self.plot_fids), fid_out=save)
        for fid in self.plot_fids:
            os.remove(fid)

    def make_evaluation_composite_pdf(self, delete_originals=True):
        """
        Utility function to combine all PDFs generated by 
        make_event_figure_pdf() into a single PDF tagged by the current 
        evaluation (iteration, step count). This is meant to make it easier
        for the User to scroll through figures. Option to delete the original
        event-specific PDFs which are now redundant

        .. note::
            This can be run without running format()

        :type delete_originals: bool
        :param delete_originals: delete original pdf files after mergin
        """
        event_figures = glob(os.path.join(self.paths.figures, "*", "*.pdf"))
        # Collecting evaluation tags, e.g., ['i01s00', 'i01s01']
        tags = set([os.path.basename(_).split("_")[0] for _ in event_figures])
        for tag in tags:
            fids = [fid for fid in event_figures if tag in fid]
            fid_out = os.path.join(self.paths.figures, f"{tag}.pdf")

            merge_pdfs(fids=sorted(fids), fid_out=fid_out)
            if delete_originals:
                for fid in fids:
                    os.remove(fid)






