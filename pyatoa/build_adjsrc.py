#!/usr/bin/env python3
"""
Invoke Pyadjoint and Pyflex on observed and synthetic data to generate misfit
windows and adjoint sources
"""
import copy

from .config import Config

class BobTheBuilder():
    """
    wrapper to contain all the adjoint source creation functionalities
    """
    def __init__(self,config,observed=None,synthetic=None,event=None):
        self.cfg = copy.deepcopy(config)
        self.observed = observed
        self.synthetic = synthetic
        self.event = event

    def run_pyflex(self):
        """
        use pyflex to grab windows

        """
        
        CD = choose_config("pyflex",PD)
        config = pyflex.Config(min_period=PD["bounds"][0],
                               max_period=PD["bounds"][1],
                               stalta_waterlevel=CD[0],
                               tshift_acceptance_level=CD[1],
                               dlna_acceptance_level=CD[2],
                               cc_acceptance_level=CD[3],
                               c_0=CD[4],c_1=CD[5],c_2=CD[6],c_3a=CD[7],
                               c_3b=CD[8],c_4a=CD[9],c_4b=CD[10])

        pf_event = pyflex.Event(latitude=event.origins[0].latitude,
                                longitude=event.origins[0].longitude,
                                depth_in_m=event.origins[0].depth,
                                origin_time=event.origins[0].time)

        pf_station = pyflex.Station(latitude=inv[0][0].latitude,
                                    longitude=inv[0][0].longitude)

        # iterate windows by component and place into dictionary output
        # create stalta data from envelopes of synthetic data
        windows,staltas = {},{}
        empties = 0
        for comp in PD["comp_list"]:
            print(comp,end='... ')

            obs,syn = breakout_stream(st.select(component=comp))
            window = pyflex.select_windows(observed=obs,
                                            synthetic=syn,
                                            config=config,
                                            event=pf_event,
                                            station=pf_station,
                                            plot=False)

            # calculate stalta
            syn_envelope = envelope(syn[0].data)
            stalta = pyflex.stalta.sta_lta(data=syn_envelope,
                                           dt=syn[0].stats.delta,
                                           min_period=PD["bounds"][0])
            staltas[comp] = stalta

            # check if pyflex is returning empty windows
            print("{} window(s)".format(len(window)))
            if not window:
                empties+=1
                continue
            windows[comp] = window

        # if all components show empty windows, raise the alarm
        if empties == len(PD["comp_list"]):
            raise Exception("Empty windows")

        # append STA/LTA water level to Par. dict. for plotting
        PD["stalta_wl"] = CD[0]

        # save windows into pyasdf file with stalta as the data and window-
        # parameter dictionaries as external information. dictionary needs
        # to be modified to work in pyasdf format
        if PD["dataset"]:
            if PD["verbose"]:print("Saving windows to PyASDF dataset")
            for comp in windows.keys():
                for i,window in enumerate(windows[comp]):
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

    def populate(self):
        """

        """

