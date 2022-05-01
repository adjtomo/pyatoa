"""
A class and associated functions that simplify calling Pyatoa functionality
with a SeisFlows workflow. Includes multiprocessing functionality to run Pyatoa
processing in parallel.
"""
import os
import pyatoa
import logging
import warnings
from glob import glob
from time import sleep
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor
from pyasdf import ASDFDataSet

from pyatoa.utils.images import merge_pdfs
from pyatoa.utils.read import read_station_codes
from pyatoa.utils.asdf.clean import clean_dataset


class IO(dict):
    """
    Dictionary with accessible attributes, used to simplify access to dicts.
    """
    def __init__(self, event_id, iter_tag, step_tag, paths, logger, codes,
                 misfit=None, nwin=None, stations=0, processed=0, exceptions=0, 
                 plot_fids=None, fix_windows=False):
        """
        Hard set required parameters here, that way the user knows what is
        expected of the IO class during the workflow.

        :type event_id: str
        :param event_id: event identifier to be passed in from Config
        :type iter_tag: str
        :param iter_tag: iteration identifier to be passed in from Config, 
            e.g., i00
        :type paths: pyatoa.core.pyaflowa.PathStructure
        :param paths: The specific path structure that Pyaflowa will use to
            navigate the filesystem, gather inputs and produce outputs.
        :type logger: logging.Logger
        :param logger: An individual event-specific log handler so that log
            statements can be made in parallel if required
        :type codes: list
        :param codes: a list of station codes that will be passed to the 
            Manager when gathering waveform data
        :type cfg: pyatoa.core.config.Config
        :param cfg: The event specific Config object that will be used to
            control the processing during each pyaflowa workflow.
        :type misfit: int
        :param misfit: output storage to keep track of the total misfit accrued
            for all the stations processed for a given event
        :type nwin: int
        :param nwin: output storage to keep track of the total number of windows
            accrued during event processing. Will be used to scale raw misfit
        :type stations: int
        :param stations: output storage to keep track of the total number of
            stations where processing was attempted. For output log statement
        :type processed: int
        :param processed: output storage to keep track of the total number of
            stations where SUCCESSFULLY processed. For output log statement
        :type exceptions: int
        :param exceptions: output storage to keep track of the total number of
            stations where processing hit an unexpected exception.
            For output log statement
        :type plot_fids: list
        :param plot_fids: output storage to keep track of the the output .pdf
            files created for each source-receiver pair. Used to merge all pdfs
            into a single output pdf at the end of the processing workflow.
        :type fix_windows: bool
        :param fix_windows: tells the processing function within the Manager
            whether or not to re-use misfit windows from a previous evaluation.
        """
        self.event_id = event_id
        self.iter_tag = iter_tag
        self.step_tag = step_tag
        self.paths = paths
        self.codes = codes
        self.logger = logger
        self.misfit = misfit
        self.nwin = nwin
        self.stations = stations
        self.processed = processed
        self.exceptions = exceptions
        self.plot_fids = plot_fids or []
        self.fix_windows = fix_windows

    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


class PathStructure:
    """
    Generalizable path structure that Pyaflowa requires to work.
    The idea is that we hardcode paths into a separate class so that the
    functionality of Pyaflowa remains independent of the path structure,
    allowing Pyaflowa to operate e.g. with SeisFlows, or independently.
    """
    def __init__(self, structure="standalone", **kwargs):
        """
        Define the necessary path structure for Pyaflowa

        .. note::
            Pyaflowa mandates the following required directory structure:

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

        :type structure: str
        :param structure: the choice of PathStructure

            * 'standalone': The default path structure that is primarily used
              for running Pyaflowa standalone, without other workflow tools.
            * 'seisflows': The path structure required when Pyaflowa is called
              by SeisFlows. Paths are hardcoded here based on the SeisFlows 
              directory structure.
        """
        # This 'constants' list mandates that the following paths exist.
        # The Pyaflowa workflow assumes that it can read/write from all of the
        # paths associated with these keys
        self._REQUIRED_PATHS = ["cwd", "data", "datasets", "figures", "logs", 
                                "ds_file", "stations_file", "responses", 
                                "waveforms", "synthetics", "adjsrcs", 
                                "event_figures"]

        # Call the available function using its string representation,
        # which sets the internal path structure.
        try:
            getattr(self, structure)(**kwargs)
        except AttributeError as e:
            raise AttributeError(
                f"{structure} is not a valid path structure") from e

    def __str__(self):
        """String representation for PathStructure, print out dict"""
        maxkeylen=max([len(_) for _ in self._REQUIRED_PATHS])

        str_out = ""
        for path in self._REQUIRED_PATHS:
            str_out += f"{path:<{maxkeylen}}: '{getattr(self, path)}'\n"
        return str_out
            
    def __repr__(self):
        """Simple call string representation"""
        return self.__str__()

    def standalone(self, workdir=None, data=None, datasets=None, figures=None,
                   logs=None, responses=None, waveforms=None, 
                   synthetics=None, adjsrcs=None, stations_file=None, 
                   event_file=None, **kwargs):
        """
        If Pyaflowa should be used in a standalone manner without external
        workflow management tools. Simply creates the necessary directory
        structure within the current working directory.
        Attributes can be used to overwrite the default path names

        .. to do:
            !!! How to work in event file recovery in standalone mode?
        """
        # General directories that all processes will write to
        self.workdir = workdir or os.getcwd()
        self.cwd = self.workdir

        self.datasets = datasets or os.path.join(self.workdir, "datasets")
        self.figures = figures or os.path.join(self.workdir, "figures")
        self.logs = logs or os.path.join(self.workdir, "logs")

        # Event-specific directories that only certain processes will write to
        self.ds_file = os.path.join(self.datasets, "{source_name}.h5")
        self.event_figures = os.path.join(self.figures, "{source_name}")
        self.synthetics = synthetics or os.path.join(self.workdir, "input",
                                                     "synthetics",
                                                     "{source_name}"
                                                     )

        # General read-only directories that may or may not contain input data
        self.data = data or os.path.join(self.workdir, "input", "DATA")
        self.responses = responses or os.path.join(self.workdir, "input", 
                                                   "responses")
        self.waveforms = waveforms or os.path.join(self.workdir, "input", 
                                                   "waveforms")

        # General write-only directories that processes will output data to
        self.adjsrcs = adjsrcs or os.path.join(self.workdir, "adjsrcs", 
                                               "{source_name}")

        # Event-specific read-only STATIONS file that defines stations to be
        # used during the workflow
        self.stations_file = stations_file or os.path.join(self.workdir, 
                                                           "{source_name}",
                                                           "STATIONS")

    def seisflows(self, **kwargs):
        """
        Hard coded SeisFlows directory structure which is created based on the
        PATH seisflows.config.Dict class, central to the SeisFlows workflow.
        """
        # Be explicit about the required argument 'sfpaths'
        PATH = kwargs.get("sfpaths", None)
        if PATH is None:
            raise TypeError("Pyaflowa SeisFlows path structure requires the "
                            "positional argument 'sfpaths' which should point "
                            "to the global PATH attribute in SeisFlows")

        self.workdir = PATH.WORKDIR
        self.cwd = os.path.join(PATH.SOLVER, "{source_name}")
        self.data = os.path.join(self.cwd, "DATA")

        self.datasets = os.path.join(PATH.PREPROCESS, "datasets")
        self.figures = os.path.join(PATH.PREPROCESS, "figures")
        self.logs = os.path.join(PATH.PREPROCESS, "logs")
        
        self.ds_file = os.path.join(self.datasets, "{source_name}.h5")
        self.event_figures = os.path.join(self.figures, "{source_name}")

        self.responses = []
        self.waveforms = [os.path.join(self.cwd, "traces", "obs")]
        if PATH.DATA is not None: 
            self.responses.append(os.path.join(PATH.DATA, "seed"))
            self.waveforms.append(os.path.join(PATH.DATA, "mseed"))

        self.synthetics = [os.path.join(self.cwd, "traces", "syn")]
        self.adjsrcs = [os.path.join(self.cwd, "traces", "adj")]
        self.stations_file = os.path.join(self.cwd, "DATA", "STATIONS")

    def format(self, mkdirs=True, **kwargs):
        """
        Paths must be source dependent in order for the directory structure
        to remain dynamic and navigable. This function provides a blanket
        formatting to all required paths making it easier to pass this
        path structure around. Returns a copy of itself so that the internal
        structure is unchanged and can be re-used.

        :type mkdirs: bool
        :param mkdirs: make directories that don't exist. this should always be
            True otherwise Pyaflowa won't work as intended if directories are
            missing, but oh well it's here.
        :rtype: pyatoa.core.pyaflowa.PathStructure
        :return: a formatted PathStructure object
        """
        assert("source_name" in kwargs.keys()), \
                f"format() missing required argument 'source_name' (pos 1)"

        # Ensure we are not overwriting the template path structure
        path_structure_copy = deepcopy(self)

        for key in self._REQUIRED_PATHS:
            # Some paths may be lists so treat everything like a list
            req_paths = getattr(self, key)
            if isinstance(req_paths, str):
                req_paths = [req_paths]

            # Format using the kwargs, str with no format braces not affected
            fmt_paths = [os.path.abspath(_.format(**kwargs)) for _ in req_paths]
      
            # Files don't need to go through makedirs but need to exist. 
            if "_file" in key:
                # Dict comp to see which files in the list dont exist, if any
                file_bool = {f: os.path.exists(f) for f in fmt_paths 
                                                       if not os.path.exists(f)}
                # Kinda hacky way to skip over requiring dataset to exist,
                # maybe try to find a more elegant way to exclude it 
                if file_bool and key != "ds_file":
                    raise FileNotFoundError(f"Paths: {file_bool.keys()} "
                                            f"must exist and doesn't, please "
                                            f"check these files")

            # Required path structure must exist. Called repeatedly but cheap
            else:
                for path_ in fmt_paths:
                    if not os.path.exists(path_):
                        if mkdirs:
                            try:
                                os.makedirs(path_)  
                            except FileExistsError:
                                # Parallel processes trying to make the same dir
                                # will throw this error. Just let it pass as we
                                # only need to make them once.
                                continue
                        else:
                            raise IOError(f"{path_} is required but does not "
                                          f"exist and cannot be made")

            # Convert single lists back to str, ugh
            if len(fmt_paths) == 1:
                fmt_paths = fmt_paths[0]
            
            # Overwrite the template structure with the formatted one
            setattr(path_structure_copy, key, fmt_paths)

        return path_structure_copy


class Pyaflowa:
    """
    A class that simplifies calling the Pyatoa waveform misfit quantification
    workflow en-masse, i.e. processing, multiple stations and multiple events
    at once.
    """
    def __init__(self, structure="standalone", config=None, plot=True, 
                 map_corners=None, log_level="DEBUG",  **kwargs):
        """
        Initialize the flow. Feel the flow.
        
        :type paths: seisflows.config.Dict
        :param paths: PREPROCESS module specific tasks that should be defined
            by the SeisFlows preprocess class. Three required keys, data,
            figures, and logs
        :type par: seisflows.config.Dict
        :param par: Parameter list tracked internally by SeisFlows
        """
        # Establish the internal workflow directories based on chosen structure
        self.structure = structure.lower()

        self.plot = plot
        self.map_corners = map_corners
        self.log_level = log_level

        # To be created by setup()
        self.path_structure = None
        self.io = None
        self.config = None
        self.mgmt = None

    def copy(self):
        """
        Convenience copy function to provide a deep copy of Pyaflowa for use 
        when multiple instances of the same object are required
        """
        return deepcopy(self)

    def setup_seisflows(self, source_name, iteration=None, step_count=None, sfpar=None,
              source_prefix="CMTSOLUTION", loc="*", cha="*", fix_windows=False, 
              multiprocess=False, config=None):
        """
        Generate a config object for a specific event id / source name
        Set the correct paths and adjust a few parameters based on the 
        location in the inversion

        Sets up the IO attribute dictionary to be carried around through the 
        processing procedure.

        .. note::
            IO object is not made an internal attribute because multiprocessing
            may require multiple, different IO objects to exist simultaneously,
            so they need to be passed into each of the functions.

        :type cwd: str
        :param cwd: current event-specific working directory within SeisFlows
        :type multiprocess: bool
        :param multiprocess: flag to turn on parallel processing for a given 
            event --- uses concurrent futures to perform multiprocessing
        :type multiprocess: bool
        :param multiprocess: If intending to use concurrent futures to 
            multiprocess an event, setup needs to know as this affects how the
            logger is setup
        :type loc: str
        :param loc: if codes is None, Pyatoa will generate station codes based 
            on the SPECFEM STATIONS file, which does not contain location info.
            This allows user to set the location values manually when building
            the list of station codes. Defaults to wildcard '??', which is 
            usually acceptable
        :type cha: str
        :param cha: if codes is None, Pyatoa will generate station codes based
            on the SPECFEM STATIONS file, which does not contain channel info. 
            This variable allows the user to set channel searching manually,
            wildcards okay. Defaults to 'HH?' for high-gain, high-sampling rate
            broadband seismometers, but this is dependent on the available data.
        :type source_prefix: str
        :param source_prefix: How source files will be prefixed, e.g.,
            CMTSOLUTION_???????? or FORCESOLUTION_??????
        :rtype: pyatoa.core.pyaflowa.IO
        :return: dictionary like object that contains all the necessary
            information to perform processing for a single event
        """
        self.path_structure = PathStructure(self.structure, **kwargs)
        paths = self.path_structure.format(source_name=source_name)

        # Create a new instance of the internal config which keeps track of
        # our current location in the workflow. Or if provided by the user, 
        # create a new config based on SeisFlows3 parameter object
        if sfpar is not None:
            self.config = pyatoa.Config(seisflows_par=sfpar)
        elif config is None:
            warnings.warn("No Config object passed, initiating empty "
                          "Config", UserWarning)
            self.config = pyatoa.Config(iteration=1, step_count=0)
        else:
            self.config = config 

        # Setup some specific values for this Pyaflowa run
        self.config.iteration = iteration
        self.config.step_count = step_count
        self.config.event_id = source_name
        self.config.paths = {"responses": paths.responses,
                        "waveforms": paths.waveforms,
                        "synthetics": paths.synthetics,
                        "events": paths.data,
                        }

        # Event-specific log files to track processing workflow. If no iteration
        # given, dont tag with iter/step, likely not an inversion scenario
        log_fid = f"{config.event_id}.log"
        if config.iter_tag is not None:
            log_fid = f"{config.eval_tag}_{log_fid}"
        log_fid = os.path.join(paths.logs, log_fid)

        if not multiprocess:
            event_logger = self._create_event_log_handler(fid=log_fid)
        else:
            # Multiprocess logging is less verbose 
            event_logger = self._create_multiprocess_log_handler(fid=log_fid)

        # Only query FDSN at the very first function evaluation
        if config.iteration != 1 and config.step_count != 0:
            config.client = None

        # Clean out any existing dataset for the current evaluation
        with ASDFDataSet(paths.ds_file) as ds:
            # Enure the ASDFDataSet has no previous data 
            clean_dataset(ds, iteration=config.iteration, 
                          step_count=config.step_count) 
            config.write(write_to=ds)
            
            # Initiate the manager and gather event, searching for source prefix
            # only, e.g., CMTSOLUTION or FORCESOLUTION
            mgmt = pyatoa.Manager(ds=ds, config=config)
            mgmt.gather(choice="event", event_id="", prefix=source_prefix)

        codes = read_station_codes(paths.stations_file, loc=loc, cha=cha)

        # Dict-like object used to keep track of information for a single event
        # processing run, simplifies information passing between functions.
        io = IO(event_id=config.event_id, iter_tag=config.iter_tag, codes=codes,
                step_tag=config.step_tag, paths=paths, logger=event_logger,
                misfit=None, nwin=None, stations=0, processed=0, exceptions=0,
                plot_fids=[], fix_windows=fix_windows)
    
        # Set internal attributes 
        self.io = io
        self.config = config

    def process_event(self, station_code=None, **kwargs):
        """
        The main processing function for Pyaflowa misfit quantification. IO
        and config should be passed in from setup()

        Processes waveform data for all stations related to a given event,
        produces waveform and map plots during the processing step, saves data
        to an ASDFDataSet and writes adjoint sources and STATIONS_ADJOINT file,
        required by SPECFEM3D's adjoint simulations, to disk.

        Kwargs passed to pyatoa.Manager.flow() function.

        :type source_name: str
        :param source_name: event id to be used for data gathering, processing
        :type station_code: str
        :param station_code: used to limit processing to a single station,
            used mostly for debug purposes. Must match part of one of the codes
            defined by 'io.codes'. e.g., 'BFZ' to match 'NZ.BFZ.*.*'
        :rtype: float or None
        :return: the total scaled misfit collected during the processing chain,
            scaled_misfit will return None if no windows have been found or
            no misfit was calculated
        """
        # Open the dataset as a context manager and process all events in serial
        with ASDFDataSet(self.io.paths.ds_file) as ds:
            self.mgmt = pyatoa.Manager(ds=ds, config=self.config)
            for code in self.io.codes:
                # Allow user to process a single station, used for debugging
                if station_code and station_code not in code:
                    continue
                self.process_station(code=code, **kwargs)

        self.finalize()

        # Scale the raw event misfit by the number of windows according to
        # Tape et al. (2010).
        try:
            scaled_misfit = 0.5 * self.io.misfit / self.io.nwin
        # Dealing with the cases where 1) nwin==0 (signifying either no windows
        # found, or calc'ing misfit on whole trace) and 2) when misfit and nwin
        # are None (no misfit found for event)
        except (TypeError, ZeroDivisionError):
            scaled_misfit = self.io.misfit

        return scaled_misfit


    def finalize(self):
        """
        Wrapper for any finalization procedures after a single event workflow
        Returns total misfit calculated during process()

        :type io: pyatoa.core.pyaflowa.IO
        :param io: dict-like container that contains processing information
        :rtype: float or None
        :param: the scaled event-misfit, i.e. total raw misfit divided by
            number of windows. If no stations were processed, returns None
            because that means theres a problem. If 0 were returned that would
            give the false impression of 0 misfit which is wrong.
        """
        # Finalization log statement for the user
        io.logger.info(f"\n{'=' * 80}\n\nSUMMARY\n\n{'=' * 80}\n"
                       f"SOURCE NAME: {io.event_id}\n"
                       f"STATIONS: {io.processed} / {io.stations}\n"
                       f"WINDOWS: {io.nwin}\n"
                       f"RAW MISFIT: {io.misfit}\n"
                       f"UNEXPECTED ERRORS: {io.exceptions}"
                       )

        self._make_event_pdf_from_station_pdfs()
        self._write_specfem_stations_adjoint_to_disk()

    def process_station(self, code, **kwargs):
        """
        Process a single seismic station for a given event. Return processed
        manager and status describing outcome of processing. Multiple error 
        catching chunks to ensure that a failed processing for a single station
        won't kill the entire job. Needs to be called by process_event()

        Kwargs passed to pyatoa.core.manager.Manager.flow()

        :type mgmt: pyatoa.core.manager.Manager
        :param mgmt: Manager object to be used for data gathering
        :type code: str
        :param code: Pyatoa station code, NN.SSS.LL.CCC
        :type io: pyatoa.core.pyaflowa.IO
        :param io: dict-like object that contains the necessary information
            to process the station
        :type fix_windows: bool
        :param fix_windows: pased to the Manager flow function, determines 
            whether previously gathered windows will be used to evaluate the 
            current set of synthetics. First passed through an internal check
            function that evaluates a few criteria before continuing.
        :rtype tuple: (pyatoa.core.manager.Manager, pyatoa.core.pyaflowa.IO)
        :return: a processed manager class, and the IO attribute class
        """
        net, sta, loc, cha = code.split(".")

        self.io.logger.info(f"\n{'=' * 80}\n\n{code}\n\n{'=' * 80}")
        self.io.stations += 1
        self.mgmt.reset()

        # Data gathering chunk; if fail, do not continue
        try:
            self.mgmt.gather(code=code)
        except pyatoa.ManagerError as e:
            self.io.logger.warning(e)
            return None

        # Data processing chunk; if fail, continue to plotting
        try:
            self.mgmt.flow(fix_windows=self.io.fix_windows)
            status_ok = True
        except pyatoa.ManagerError as e:
            self.io.logger.warning(e)
            status_ok = False
            pass
        except Exception as e:
            # Uncontrolled exceptions should be noted in more detail
            self.io.logger.warning(e, exc_info=True)
            self.io.exceptions += 1
            status_ok = False
            pass

        # Plotting chunk; fid is e.g. path/i01s00_NZ_BFZ.pdf
        if self.plot:
            plot_fid = "_".join([self.mgmt.config.iter_tag,  # e.g., i01
                                 self.mgmt.config.step_tag,  # e.g., s00
                                 net,                        # e.g., NZ
                                 sta + ".pdf"                # e.g., BFZ.pdf
                                 ])
            save = os.path.join(self.io.paths.event_figures, plot_fid)
            self.io.logger.info(f"saving figure to: {save}")
            # Mapping may fail if no Inventory or Event object is attached
            # Waveform figures are expected to work if we have gotten this far
            try:
                self.mgmt.plot(corners=self.map_corners, show=False, save=save)
            except pyatoa.ManagerError as e:
                self.io.logger.warning("cannot plot map, plotting waveform only")
                self.mgmt.plot(choice="wav", show=False, save=save)

            # If a plot is made, keep track so it can be merged later on
            self.io.plot_fids.append(save)

        # Finalization chunk; only if processing is successful
        if status_ok:
            # We allow misfit and nwin to be None to signify that no misfits
            # have been calculated. But if we calculate something then we
            # need to overwrite the None values
            if self.io.misfit is None:
                self.io.misfit = 0
            if self.io.nwin is None:
                self.io.nwin = 0

            # Keep track of outputs for the final log summary and misfit value
            self.io.misfit += self.mgmt.stats.misfit
            # The deal with the case where window selection is skipped and 
            # mgmt.stats.nwin == None
            num_new_win = self.mgmt.stats.nwin or 0
            self.io.nwin += num_new_win
            self.io.processed += 1
            
            # SPECFEM wants adjsrcs for each comp, regardless if it has data
            self.mgmt.write_adjsrcs(path=self.io.paths.adjsrcs, 
                                    write_blanks=True)

    def _process_event_multiprocess_true(self, *args, **kwargs):
        """
        A hacky way to get around problem in passing additional kwargs through
        ProcessPoolExecutor. Simply define a new process_event function that 
        hardcodes the 'multiprocess' parameter controlling logging. To be called
        by multi_event_process() only.
        """
        return self.process_event(*args, **kwargs, multiprocess=True)
        
    def multi_event_process(self, source_names, max_workers=None, **kwargs):
        """
        Use concurrent futures to run the process() function in parallel.
        This is a multiprocessing function, meaning multiple instances of Python
        will be instantiated in parallel.

        :type source_names: list of str
        :param solver_dir: a list of all the source names to process. each will
            be passed to process()
        :type max_workers: int
        :param max_workers: maximum number of parallel processes to use. If
            None, automatically determined by system number of processors.
        """
        misfits = {}

        print(f"Beginning parallel processing of {len(source_names)} events...")
    
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            for source_name, misfit in zip(
                    source_names, 
                    executor.map(self._process_event_multiprocess_true,
                                 source_names)
                    ):
                misfits[os.path.basename(source_name)] = misfit
        return misfits

    def _write_specfem_stations_adjoint_to_disk(self):
        """
        Create the STATIONS_ADJOINT file required by SPECFEM to run an adjoint
        simulation. Should be run after all processing has occurred. Works by
        checking what stations have adjoint sources available and re-writing the
        existing STATIONS file that is on hand.

        :type cwd: str
        :param cwd: current SPECFEM run directory within the larger SeisFlows
            directory structure
        """
        self.io.logger.info("generating STATIONS_ADJOINT file for SPECFEM")

        # These paths follow the structure of SeisFlows and SPECFEM
        adjoint_traces = glob(os.path.join(self.io.paths.adjsrcs, "*.adj"))

        # Simply append to end of file name e.g. "path/to/STATIONS" + "_ADJOINT"
        stations_adjoint = self.io.paths.stations_file  + "_ADJOINT"

        # Determine the network and station names for each adjoint source
        # Note this will contain redundant values but thats okay because we're
        # only going to check membership using it, e.g. [['NZ', 'BFZ']]
        adjoint_stations = [os.path.basename(_).split(".")[:2] for _ in
                            adjoint_traces]

        # The STATION file is already formatted so leverage for STATIONS_ADJOINT
        lines_in = open(self.io.paths.stations_file, "r").readlines()

        with open(stations_adjoint, "w") as f_out:
            for line in lines_in:
                # Station file line format goes: STA NET LAT LON DEPTH BURIAL
                # and we only need: NET STA
                check = line.split()[:2][::-1]
                if check in adjoint_stations:
                    f_out.write(line)

    def _make_event_pdf_from_station_pdfs(self):
        """
        Combine a list of single source-receiver PDFS into a single PDF file
        for the given event.

        :type fids: list of str
        :param fids: paths to the pdf file identifiers
        :type output_fid: str
        :param output_fid: name of the output pdf, will be joined to the figures
            path in this function
        """
        if self.io.plot_fids:
            self.io.logger.info("creating single .pdf for all output figures")
            # e.g. i01s00_2018p130600.pdf
            output_fid = f"{io.iter_tag}{io.step_tag}_{io.event_id}.pdf"

            # Merge all output pdfs into a single pdf, delete originals
            save = os.path.join(self.io.paths.event_figures, output_fid)
            merge_pdfs(fids=sorted(self.io.plot_fids), fid_out=save)

            for fid in self.io.plot_fids:
                os.remove(fid)

    def _create_event_log_handler(self, fid):
        """
        Create a separate log file for each embarassingly parallel event.
        SeisFlows just runs individual Pyaflowa instances so we can use the 
        package-wide logger.

        :type fid: str
        :param fid: the name of the outputted log file
        :rtype: logging.Logger
        :return: an individualized logging handler
        """
        logger = pyatoa.logger

        # Propogate logging to individual log files, always overwrite
        handler = logging.FileHandler(fid, mode="w")

        # Maintain the same look as the standard console log messages
        logfmt = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
        formatter = logging.Formatter(logfmt, datefmt="%Y-%m-%d %H:%M:%S")
        handler.setFormatter(formatter)
        handler.setLevel(self.log_level.upper())

        for log in ["pyflex", "pyadjoint", "pyatoa"]:
            # Set the overall log level
            logger = logging.getLogger(log)
            # Turn off any existing handlers (stream and file) so that we can 
            # direct all logging to a single file
            while logger.hasHandlers():
                logger.removeHandler(logger.handlers[0])

            logger.setLevel(self.log_level.upper())
            logger.addHandler(handler)

        return logger

    def _create_multiprocess_log_handler(self, fid):
        """
        Create a separate logger and file for each multiprocessed event.

        .. note::
            Because of the way Pyatoa logging is written, we can't propagate the
            log statements within the package into Pyaflowa, it just leads to 
            all the processes writing everything to the same file.

            Instead we turn off the package logger completely and instead create
            a new 'Pyaflowa' logger which will be more restricted in the info
            it can provide, but will be able to log separately for each multi
            process

        :type fid: str
        :param fid: the name of the outputted log file
        :rtype: logging.Logger
        :return: an individualized logging handler
        """
        # Turn off the package-wide logger by setting to strictest mode and
        # removing all existing handlers
        logger = pyatoa.logger
        while logger.hasHandlers():
            logger.removeHandler(logger.handlers[0])
        logger.setLevel("CRITICAL")

        # Create a new file-specific logger
        logger = logging.getLogger(f"pyaflowa")
        logfmt = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
        formatter = logging.Formatter(logfmt, datefmt="%Y-%m-%d %H:%M:%S")

        # Set handler with specific file name
        handler = logging.FileHandler(fid, mode="w")
        handler.setFormatter(formatter)

        logger.setLevel(self.log_level)
        logger.addHandler(handler)

        return logger
