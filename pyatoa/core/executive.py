#!/usr/bin/env python3
"""
A class and functionality to direct multiple Managers to parallelize data
processing in Pyatoa using the in-built concurrent.futures package
"""
import os
from concurrent.futures import ProcessPoolExecutor
from pyasdf import ASDFDataSet
from pyatoa import logger, Manager, Config

logger.setLevel("CRITICAL")


class Executive:
    """
    The Executive is hierarchically above Pyatoa's core class, the Manager. 
    It sets up a simple framework to organize and parallelize misfit 
    quantification.
    """
    def __init__(self, event_ids, station_codes, config, max_stations=4,
                 max_events=1, max_workers=4):
        """
        The Directory needs some key information before it can run processing
        """
        self.config = config

        self.station_codes = station_codes
        self.event_ids = event_ids

        self.max_stations = max_stations
        self.max_events = max_events
        self.max_workers = max_workers

        # Define a rudimentary path structure to keep main dir. light
        self.cwd = os.getcwd()
        self.datasets = os.path.join(self.cwd, "datasets")
        self.figures = os.path.join(self.cwd, "figures")

        for path in [self.datasets, self.figures]:
            if not os.path.exists(path):
                os.mkdir(path)


        self.check()

    def check(self):
        """
        Parameter checking
        """
        assert(self.max_stations * self.max_events <= self.max_workers), (
            "max_stations * max_events > max_workers, please reduce " 
            "max_events or max_stations, or increase allowable max_workers"
        )
        # Ensure entries are lists
        if not isinstance(self.station_codes, list):
            self.station_codes = [self.station_codes]
        if not isinstance(self.event_ids, list):
            self.event_ids = [self.event_ids]

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
        """
        station_misfits = {}
        # Workaround for having it be pretty difficult to pass multiple
        # arguments into an executor, just pass them in together and let the
        # underlying function break them apart again
        codes = [f"{event_id}-{_}" for _ in self.station_codes]
        with ProcessPoolExecutor(max_workers=self.max_stations) as executor:
            for code, misfit in zip(codes, executor.map(self.process_station,
                                                        codes)):
                # Tag dictionary entries by station code
                station_misfits[code.split("-")[1]] = misfit

        return station_misfits

    def process_station(self, event_id_and_station_code):
        """
        Parallel process multiple Managers simultaneously
        """
        event_id, station_code = event_id_and_station_code.split("-")
        ds = ASDFDataSet(os.path.join(self.datasets, f"{event_id}.h5"))

        config = self.config.copy()
        config.event_id = event_id

        mgmt = Manager(config=config, ds=ds)
        mgmt.gather(code=station_code)
        mgmt.flow()
        mgmt.plot(choice="both", show=False, save=os.path.join(
            self.figures, f"{event_id}_{station_code}.png")
                  )

        return mgmt.stats.misfit

