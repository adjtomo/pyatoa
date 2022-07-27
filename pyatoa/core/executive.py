#!/usr/bin/env python3
"""
A class and functionality to direct multiple Managers to parallelize data
processing in Pyatoa using the in-built concurrent.futures package

.. rubric::
    >>> from pyatoa import Executive
    >>> exc = Executive(event_ids=..., station_codes=..., config=...)
    >>> misfits = exc.process()

.. note::
    I was getting the following error randomly on my runs, I think it had to
    do with the fact that I was requesting all of my cores to do the job (16).
    Dropped the max number of jobs to 8 and things ran fine. Will investigate

    concurrent.futures.process.BrokenProcessPool:
    A process in the process pool was terminated abruptly while the future
    was running or pending.
"""
import os
import time
import logging
import random
from concurrent.futures import ProcessPoolExecutor
from logging.handlers import MemoryHandler
from pyasdf import ASDFDataSet
from pyatoa import Manager


class Executive:
    """
    The Executive is hierarchically above Pyatoa's core class, the Manager. 
    It sets up a simple framework to organize and parallelize misfit 
    quantification.
    """
    def __init__(self, event_ids, station_codes, config, max_stations=4,
                 max_events=1, cat="+", log_level="DEBUG", cwd=None,
                 datasets=None, figures=None, logs=None, adjsrcs=None,
                 ds_fid_template=None):
        """
        The Executor needs some key information before it can run processing

        :type event_ids: list
        :param event_ids: List of event identifiers used to tag outputs
        :type station_codes: list
        :param station_codes: list of station codes used to gather waveform and
            metadata. Needs to be in the format NN.SSS.LL.CCC where N=network,
            S=station, L=location, C=channel. Wildcards okay
        :type config: pyatoa.core.config.Config
        :param config: Config object which will control ALL processing
        :type max_stations: int
        :param max_stations: maximum number of stations to process
            simultaneously for one event
        :type max_events: int
        :param max_events: number of max events to process in parallel.
            Works together with `max_stations` such that the total number
            of concurrent processes is `max_events` * `max_stations`
        :type cat: str
        :param cat: a string conCATenator that will be used to join together
            event ids and station codes. Can be anything but try to make it
            uncommon as it is given to str.split()
        :type log_level: str
        :param log_level: log level to be given to all underlying loggers
        :type cwd: str
        :param cwd: active working directory to look for inputs and save output
        :type datasets: str
        :param datasets: path to save ASDFDataSets. defaults to a subdirectory
            'datasets', inside the current working directory.
        :type figures: str
        :param figures: path to save output figures. defaults to a subdirectory
            'figures', inside the current working directory.
        :type logs: str
        :param logs: path to save text log files. defaults to a subdirectory
            'logs', inside the current working directory.
        :type adjsrcs: str
        :param adjsrcs: path to save text adjoint source text files.
            defaults to a subdirectory 'adjsrcs', inside the current working
            directory.
        """
        self.config = config

        self.station_codes = station_codes
        self.event_ids = event_ids

        self.max_stations = max_stations
        self.max_events = max_events

        # Define a rudimentary path structure to keep main dir. light.
        # These can usually be defaults but if working with Pyaflowa, allow
        # them to be overwritten
        self.cwd = cwd or os.getcwd()
        self.datasets = datasets or os.path.join(self.cwd, "datasets")
        self.figures = figures or os.path.join(self.cwd, "figures")
        self.logs = logs or os.path.join(self.cwd, "logs")
        self.adjsrcs = adjsrcs or os.path.join(self.cwd, "adjsrcs")
        self.ds_fid_template = ds_fid_template or os.path.join(self.datasets,
                                                               "{event_id}.h5")

        for path in [self.datasets, self.figures, self.logs]:
            if not os.path.exists(path):
                os.mkdir(path)

        self.cat = cat
        self.log_level = log_level

        self.check()

    @property
    def codes(self):
        """
        Define a set of event-station codes that are used to traverse through
        all possible source receiver combinations.

        .. note::
            Workaround for having it be pretty difficult to pass multiple
            arguments into an executor. Just pass a list of strings that is
            split by the parallel processes
        """
        codes = []
        for sta in self.station_codes:
            for eid in self.event_ids:
                codes.append(f"{eid}{self.cat}{sta}")
        return sorted(codes)

    def check(self):
        """
        Parameter checking
        """
        # Ensure entries are lists
        if not isinstance(self.station_codes, list):
            self.station_codes = [self.station_codes]
        self.station_codes = sorted(self.station_codes)

        if not isinstance(self.event_ids, list):
            self.event_ids = [self.event_ids]
        self.event_ids = sorted(self.event_ids)

        for sta in self.station_codes:
            assert(len(sta.split(".")) == 4), (f"station codes must be in " 
                                               f"format: NN.SSS.LL.CCC")

    def process(self):
        """
        Process all events concurrently
        """
        event_misfits = {}
        with ProcessPoolExecutor(max_workers=self.max_events) as executor:
            for event_id, station_misfits in zip(
                    self.event_ids, executor.map(self.process_event,
                                                 self.event_ids)):
                event_misfits[event_id] = station_misfits

        return event_misfits

    def process_event(self, event_id):
        """
        Process all given stations concurrently for a single event

        :type event_id: str
        :param event_id: one value from the Executor.events list specifying
            a given event to process
        """
        station_misfits = {}

        # Subset internal code list by event id
        codes = [code for code in self.codes if event_id in code]
        with ProcessPoolExecutor(max_workers=self.max_stations) as executor:
            for code, misfit in zip(codes, executor.map(self.process_station,
                                                        codes)):
                # Tag dictionary entries by station code
                station_misfits[code.split(self.cat)[1]] = misfit

        return station_misfits

    def process_station(self, event_id_and_station_code):
        """
        Parallel process multiple Managers simultaneously, which is the biggest
        time sync. IO is done in serial to get around BlockingIO

        .. note::
            Very broad exceptions to keep process running smoothly, you will
            need to check log messages individually to figure out if and where
            things did not work

        .. note::
            Employs a workaround to inability to parallel write to HDF5 files
            BlockingIOError by doing the processing first, and then waiting
            for each process to finish writing before accessing.

        :type event_id_and_station_code: str
        :param event_id_and_station_code: a string concatenation of a given
            event id and station code, which will be used to process a single
            source receiver pair
        """
        # Using the event id and station code for indexing and tag information
        event_id, station_code = event_id_and_station_code.split(self.cat)
        net, sta, _, _ = station_code.split(".")
        filename = f"{event_id}_{net}_{sta}"
        idx = self.codes.index(event_id_and_station_code)
        rank = self._check_rank(event_id_and_station_code)

        print(f"processing {idx}/{len(self.codes)}: {event_id} {station_code}")

        logger, memhandler = self._generate_logger(
            os.path.join(self.logs, f"{filename}.txt")
        )

        config = self.config.copy()
        config.event_id = event_id

        mgmt = Manager(config=config)
        # Data gathering break will not allow further processing
        try:
            mgmt.gather(code=station_code, syn_dir_template=event_id)
        except Exception as e:
            logger.warning(e)
            return None
        # Processing break will allow writing waveforms and plotting
        try:
            mgmt.flow()
            mgmt.write_adjsrcs(path=self.paths.adjsrcs, write_blanks=True)
        except Exception as e:
            logger.warning(e)
            pass
        # Plotting break will allow writing waveforms
        try:
            mgmt.plot(choice="both", show=False,
                      save=os.path.join(self.figures, f"{filename}.png")
                      )
        except Exception as e:
            logger.warning(e)
            pass

        # Wait till the very end to write to the HDF5 file, then do it serially
        # as it should be quick
        while True:
            try:
                # Default dataset name needs to be formatted, but user-defined
                # filenames may not, and will not be affected by format()
                ds_fid = self.ds_fid_template.format(event_id=event_id)
                with ASDFDataSet(ds_fid) as ds:
                    mgmt.ds = ds
                    mgmt.write()
                    if rank == 0:
                        config.write(ds)
                    break
            except BlockingIOError:
                # Random sleep time to decrease chances of two processes
                # attempting to access at exactly the same time
                time.sleep(random.random())

        memhandler.flush()
        return mgmt.stats.misfit

    def _check_rank(self, event_id_and_station_code):
        """
        Poor man's method for determining the processor rank for a given event.
        Used so that processes that happen only once (e.g., writing config) are
        done consistently by one process

        :type event_id_and_station_code: str
        :param event_id_and_station_code: a string concatenation of a given
            event id and station code, which will be used to process a single
            source receiver pair
        :rtype: int
        :return: rank index in Executive.codes based on event and station
        """
        event_id, station_code = event_id_and_station_code.split(self.cat)
        event_codes = sorted([code for code in self.codes if event_id in code])
        return event_codes.index(event_id_and_station_code)

    def _generate_logger(self, log_path):
        """
        Create a log file for each source. No stream handler, only file output
        Also create a memory handler to dump all log messages at once, rather
        than as they happen, allowing multiple stations to write to the same
        file sequentially

        :type log_path: str
        :param log_path: path and filename to save log file
        """
        filehandler = logging.FileHandler(log_path, mode="w")
        memhandler = MemoryHandler(capacity=1024 * 100,
                                   flushLevel=self.log_level,
                                   target=filehandler
                                   )

        # Maintain the same look as the standard console log messages
        logfmt = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
        formatter = logging.Formatter(logfmt, datefmt="%Y-%m-%d %H:%M:%S")
        filehandler.setFormatter(formatter)
        filehandler.setLevel(self.log_level)

        for log in ["pyflex", "pyadjoint", "pyatoa"]:
            # Set the overall log level
            logger = logging.getLogger(log)
            # Turn off any existing handlers (stream and file for all packages)
            while logger.hasHandlers():
                logger.removeHandler(logger.handlers[0])

            logger.setLevel(self.log_level)
            logger.addHandler(filehandler)
            # logger.addHandler(memhandler)

        return logger, memhandler