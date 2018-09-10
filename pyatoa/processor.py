#!/usr/bin/env python3
"""
Invoke Pyadjoint and Pyflex on observed and synthetic data to generate misfit
windows and adjoint sources
"""
import copy
import warnings

import pyflex
from obspy.signal.filter import envelope


from .config import Config
from .data_gather import Gatherer

class Processer():
    """
    wrapper to contain all the adjoint source creation functionalities
    """
    def __init__(self,config,ds):
        self.cfg = copy.deepcopy(config)
        self.ds = ds
        self.gatherer = None

    def gather_data(self,station_code):
        """
        call on data_gatherer to populate data space
        """
        if self.gatherer is None:
            gatherer = Gatherer(cfg=self.cfg,ds=self.ds)
            self.gatherer = gatherer
        self.gatherer.gather_all(station_code)

    def run_pyflex(self):
        """
        call pyflex module to calculate best fitting misfit windows
        """
        pf_config = pyflex.Config(
            min_period=self.cfg.min_period, max_period=self.cfg.max_period,
            stalta_waterlevel=self.cfg.pyflex_config[0],
            tshift_acceptance_level=self.cfg.pyflex_config[1],
            dlna_acceptance_level=self.cfg.pyflex_config[2],
            cc_acceptance_level=self.cfg.pyflex_config[3],
            c_0=self.cfg.pyflex_config[4], c_1=self.cfg.pyflex_config[5],
            c_2=self.cfg.pyflex_config[6], c_3a=self.cfg.pyflex_config[7],
            c_3b=self.cfg.pyflex_config[8], c_4a=self.cfg.pyflex_config[9],
            c_4b=self.cfg.pyflex_config[10]
            )
        pf_event = pyflex.Event(
            latitude=self.gatherer.event.origins[0].latitude,
            longitude=self.gatherer.event.origins[0].longitude,
            depth_in_m=self.gatherer.event.origins[0].depth,
            origin_time=self.gatherer.event.origins[0].time
            )
        pf_station = pyflex.Station(
            latitude=self.gatherer.inv[0][0].latitude,
            longitude=self.gatherer.inv[0][0].longitude
            )

        empties = 0
        windows, staltas = {}, {}
        for COMP in self.cfg.component_list:
            # TODO: add log statement here
            window = pyflex.select_windows(
                observed=self.gatherer.st_obs.select(component=COMP),
                synthetic=self.gatherer.st_syn.select(component=COMP),
                config=pf_config, event=pf_event, station=pf_station,
                )
            stalta = pyflex.stalta.sta_lta(
                data=envelope(
                    self.gatherer.st_syn.select(component=COMP)[0].data),
                dt=self.gatherer.st_syn.select(component=COMP)[0].stats.delta,
                min_period=self.cfg.min_period
                )
            staltas[COMP] = stalta

            # TODO: add log statement to print out how many windows found
            print("{} window(s)".format(len(window)))
            if not window:
                empties += 1
                continue
            windows[COMP] = window

        # if all components show empty windows, raise the alarm
        if empties == len(self.cfg.component_list):
            warnings.warn("Empty windows", UserWarning)


        if PD["verbose"]: print("Saving windows to PyASDF dataset")
            for comp in windows.keys():
                for i, window in enumerate(windows[comp]):
                    internalpath = "{net}/{sta}_{comp}_{i}_{m}".format(
                        evid=PD["event_id"],
                        net=PD["network"],
                        sta=PD["station"],
                        m=PD["model"],
                        comp=comp,
                        i=i)
                    # auxiliary data requires a data object, even though we only
                    # want the window parameter dictionary. to save on space
                    winnDixie = create_window_dictionary(window)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        PD["dataset"].add_auxiliary_data(data=np.array([True]),
                                                         data_type="MisfitWindows",
                                                         path=internalpath,
                                                         parameters=winnDixie)

        return windows, staltas, PD