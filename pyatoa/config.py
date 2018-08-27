#!/usr/bin/env python3
"""
Configuration for pyatoa runs
"""

class Config():
    def __init__(self,model_number,event_id,stations,min_period=10,
                 max_period=30,rotate_to_rtz=True,unit_output='DISP',pre_p=20,
                 post_p=120,pyflex_config='UAF',
                 adj_src_type='multitaper_misfit',plot_waveform=True,
                 plot_map=True,plot_faults_on_map=False,
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
        self.pre_p = pre_p
        self.post_p = post_p
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
        self.internal_paths = {"synthetic":synthetic_path,"pyasdf":pyasdf_path,
            "quakeml":quakeml_path,"station_list":station_list_path,
            "auxiliary_mseeds":aux_mseed_path}

    def set_internal_pathing(self):
        """
        if working with internally located data, set the pathnames here
        """
        self.internal_files = {}

    def _generate_component_list(self):
        """
        """

    def _generate_pyflex_config(self):
        """
        easy way to quickly generate a pyflex config file
        """
        import pyflex
        if not isinstance(self.pyflex_config,list):
            cfgdict = {"default":[.08,15.,1.,.8,.7,4.,0.,1.,2.,3.,10.],
                       "UAF":[.18,4.,1.5,.71,.7,2.,0.,3.,2.,2.5,12.],
                       "NZ":[]}
            CD = cfgdict[pyflex_config]

        config = pyflex.Config(min_period=PD["bounds"][0],
                               max_period=PD["bounds"][1],
                               stalta_waterlevel=CD[0],
                               tshift_acceptance_level=CD[1],
                               dlna_acceptance_level=CD[2],
                               cc_acceptance_level=CD[3],
                               c_0=CD[4],c_1=CD[5],c_2=CD[6],c_3a=CD[7],
                               c_3b=CD[8],c_4a=CD[9],c_4b=CD[10])
        


    def _write_to_pyasdf(self,ds):
