"""
Currently under development 27.09.20

Multiprocessing functionality to run Pyatoa in parallel within SeisFlows.
"""
import os
import pyatoa
import logging

from copy import deepcopy
from pyasdf import ASDFDataSet
from pyatoa.utils.asdf.clean import clean_dataset
from concurrent.futures import ProcessPoolExecutor


class Output(dict):
    """
    Dictionary with accessible attributes, used for holding processing stats.
    """
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


class Pyaflowa:
    """
    A multiprocessing class for parallel integration of the Pyatoa workflow into
    larger workflow tools. Allows multiprocessing of events.

    ..note::
        Multithreading of individual events on a per-station basis was partially
        tested, and yielded no increase in speed. It may help with external
        data collection since these processes need to wait for webservice
        queries, but generally the overhead of threading made it slower than
        simply running this task in serial.
    """
    def __init__(self, data, figures):
        """
        Initialize the flow by establishing the directory structures present
        in the SeisFlows workflow
        """
        self.data = data
        self.figures = figures
        self.plot = True
        self.corners = None
        self.fixed = True
        self.config = None

    @staticmethod
    def create_log_handler(fid):
        """
        Create a separate log file for each multiprocessed event.

        :type fid: str
        :param fid: the name of the outputted log file
        :rtype: logging.Logger
        :return: an individualized logging handler
        """
        for log in ["pyatoa", "pyflex", "pyadjoint"]:
            # Set the overall log level
            logger = logging.getLogger(log)
            logger.setLevel("DEBUG")
            # Propogate logging to individual log files
            handler = logging.FileHandler(fid)
            # Maintain the same look as the standard console log messages
            logfmt = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
            formatter = logging.Formatter(logfmt,  datefmt="%Y-%m-%d %H:%M:%S")
            handler.setFormatter(formatter)
            logger.addHandler(handler)

        return logger

    def process_event_single(self, cwd):
        """
        A template processing function to create/open an ASDFDataSet, process
        all stations related to the event, and write the associated adjoint
        sources to disk.
        """
        # Get the event id from the current working directory
        event_id = os.path.basename(cwd)

        # Copy in the Config object and resassign the event id
        config = deepcopy(self.config)
        config.event_id = event_id

        # Create the individualized logger
        logger = self.create_log_handler(fid=f"{config.event_id}.log")

        # Dict-like object used to keep track of processing outputs
        output = Output(misfit=0, nwin=0, stations=0, processed=0, exceptions=0,
                        plot_fids=[])

        # Read in the STATIONS file from the current working directory
        inv = pyatoa.read_station(os.path.join(cwd, "DATA", "STATIONS"))

        # Open the dataset as a context manager and process all events
        with ASDFDataSet(os.path.join(self.data, f"{config.event_id}")) as ds:
            clean_dataset(ds, iteration=config.iteration,
                          step_count=config.step_count
                          )
            config.write(write_to=ds)
            mgmt = pyatoa.Manager(ds=ds, config=config)
            for net in inv:
                for sta in net:
                    # Initialization of Manager and logging
                    mgmt.reset()
                    output.stations += 1
                    logger.info(
                        f"\n{'=' * 80}\n\n{net.code}.{sta.code}\n\n{'=' * 80}"
                    )
                    # Data gathering chunk; if fail, do not continue w/ workflow
                    try:
                        mgmt.gather(code=f"{net.code}.{sta.code}.*.HH?")
                    except pyatoa.ManagerError as e:
                        logger.warning(e)
                        continue
                    # Data processing chunk; if fail, continue to allow plotting
                    try:
                        mgmt.flow(fix_windows=self.fixed)
                        mgmt.write_adjsrcs(path=os.path.join(cwd, "traces",
                                                             "adj")
                                           )
                        output.misfit += mgmt.stats.misfit
                        output.nwin += mgmt.statse.nwin
                        output.processed += 1
                    except pyatoa.ManagerError as e:
                        logger.warning(e)
                        pass
                    except Exception as e:
                        # These are uncontrolled exceptions and should be noted
                        logger.warning(e, exc_info=True)
                        output.exceptions += 1
                    # Plotting chunk; fid is e.g. path/i01s00_NZ_BFZ.pdf
                    if self.plot:
                        plot_fid = "_".join([config.iter_tag, config.step_tag,
                                             net.code, sta.code + ".pdf"]
                                            )
                        save = os.path.join(self.figures, config.event_id,
                                            plot_fid)
                        mgmt.plot(corners=self.corners, show=False, save=save)
                        output.plot_fids.append(save)

        # Finalization chunk
        if output.plot_fids:
            # Merge all output pdfs into a single pdf, delete originals
            plot_fid = "_".join([config.iter_tag, config.step_tag,
                                 config.event_id + ".pdf"])
            save = os.path.join(self.figures, plot_fid)
            pyatoa.utils.images.merge_pdfs(fids=sorted(output.plot_fids),
                                           fid_out=save)
            for fid in output.plot_fids:
                os.remove(fid)

        # Record summary information at the end of log file
        pyatoa.logger.info(f"\n{'=' * 80}\n\nSUMMARY\n\n{'=' * 80}\n"
                           f"SOURCE NAME: {config.event_id}\n"
                           f"STATIONS: {output.processed} / {output.stations}\n"
                           f"WINDOWS: {output.nwin}\n"
                           f"RAW MISFIT: {output.misfit:.2f}\n"
                           f"UNEXPECTED ERRORS: {output.exceptions}"
                           )

        return

    def multiprocess_events(self, source_names):
        """
        Use multiprocessing to run event processing functionality in parallel.
        Max workers is intentionally left blank so that it can be automatically
        determined by the number of processors.
        """
        config = self.create_config()

        with ProcessPoolExecutor() as executor:
            for source, result in zip(source_names,
                                      executor.map(self.process_event_single,
                                                   cwd)):

                pass

        return




