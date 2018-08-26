#!/usr/bin/env python3
"""
Configuration for pyatoa runs
"""

class Config():
    def __init__(self,model_number,event_id,stations,min_period=10,
                 max_period=30,rotate_to_rtz=True,unit_output='DISP',
                 pyflex_config='UAF',adj_src_type='multitaper_misfit',
                 plot_waveform=True,plot_map=True,plot_faults_on_map=False,
                 show_plots=False,save_plots=False,save_pyasdf=False,
                 save_adj_src_separate=False,verbose=False,log=False):
        """
        Configuration file for pyatoa

        Allows the user to control the parameters of the packages called within
        pyatoa, as well as control where the outputs (i.e. pyasdf and plots) are
        sent after processing occurs

        :type model_number: str
        :param model_number: iteration number
        :type event_id: str
        :param event_id: unique event identifier
        :type station: str
        :param station:
        :param minimum_filter_period:
        :param maximum_filter_period:
        :param rotate_to_rtz:
        :param unit_output:
        :param pyflex_config:
        :param adjoint_src_type:
        :param plot_waveform:
        :param plot_map:
        :param plot_faults_on_map:
        :param show_plots:
        :param save_plots:
        :param save_pyasdf:
        :param save_adj_src_separate:
        :param verbose:
        :param log:
        """
        self.model_number = model_number
        self.event_id = event_id
        self.stations = stations
        self.min_period = min_period
        self.max_period = max_period
        self.rotate_to_rtz = rotate_to_rtz
        if unit_output.upper() not in ["DISP","VEL","ACC"]:
            raise PyatoaError("unit output must be 'DISP','VEL' or 'ACC'")
        self.pyflex_config = pyflex_config
        self.adj_src_type = adj_src_type
        self.plot_waveform = plot_waveform
        self.plot_map = plot_map
        self.plot_faults_on_map = plot_faults_on_map
        self.show_plots = show_plots
        self.save_plots = save_plots
        self.save_pyasdf = save_pyasdf
        self.save_adj_src_separate = save_adj_src_separate
        self.verbose = verbose
        self.log = log