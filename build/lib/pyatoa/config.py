#!/usr/bin/env python3
"""
Configuration for pyatoa runs
"""
class ConfigException(Exception):
    """
    catch-all exception for config checking
    """
    pass

class Config():
    def __init__(self,model_number,event_id,min_period=10,
                 max_period=30,rotate_to_rtz=True,unit_output='DISP',
                 pyflex_config='UAF',adj_src_type='multitaper_misfit',
                 save_pyasdf=False,verbose=False,log=False):
        """
        Configuration file for pyatoa

        Allows the user to control the parameters of the packages called within
        pyatoa, as well as control where the outputs (i.e. pyasdf and plots) are
        sent after processing occurs

        :type model_number: str
        :param model_number: iteration number
        :type event_id: str
        :param event_id: unique event identifier
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

        example call for config
        config = Config(model_number=0,event_id='2014p240655',
            stations=['NZ.BFZ.10.HH?','NZ.OPRZ.10.HH?','XX.RD01..HH?'],
            min_period=10,max_period=30,rotate_to_rtz=True,unit_output="DISP",
            pyflex_config='UAF',adj_src_type='multitaper_misfit',
            save_pyasdf=True,verbose=True,log=True,
            path_to_waveforms='/Users/chowbr/Documents/subduction/seismic/',
            path_to_pyasdf='/Users/chowbr/Documents/subduction/data/ADJTOMO/PYASDF'

        """
        self.model_number = 'm{:0>2}'.format(model_number)
        self.event_id = event_id
        self.min_period = float(min_period)
        self.max_period = float(max_period)
        self.filter_corners = float(filter_corners)
        self.rotate_to_rtz = rotate_to_rtz
        self.unit_output = unit_output.upper()
        self.pyflex_config = pyflex_config
        self.adj_src_type = adj_src_type
        self.save_pyasdf = save_pyasdf
        self.verbose = verbose
        self.log = log
        self.path_to_waveforms = path_to_waveform
        self._generate_component_list()
        self._check_config_parameters()

    def _generate_component_list(self):
        """
        create a small list for easy access to orthogonal components
        """
        if self.rotate_to_rtz:
            self.component_list = ['Z','R','T']
        else:
            self.component_list = ['Z','N','E']

    def _check_config_parameters(self):
        """
        just make sure that some of the configuration parameters are set proper
        :return:
        """
        if self.pyflex_config not in ['UAF', 'default', 'NZ']:
            raise ConfigException(
                "pyflex config must be 'UAF','default' or 'NZ")
        if self.unit_output not in ['DISP', 'VEL', 'ACC']:
            raise ConfigException("unit output must be 'DISP','VEL' or 'ACC'")

    def write_to_pyasdf(self,ds):
        """
        save the config values as a dictionary in the pyasdf data format
        for easy lookback
        """
        from obspy import UTCDateTime
        par_dict = {"time":UTCDateTime(),
                    "model_number":self.model_number,
                    "event_id":self.event_id,
                    "min_period":self.min_period,
                    "max_period":self.max_period,
                    "filter_corners":self.filter_corners,
                    "rotate_to_rtz":self.rotate_to_rtz
                    "unit_output":self.unit_output,
                    "pyflex_config":self.pyflex_config,
                    "adj_src_type":self.adj_src_type
                    }
        ds.add_auxiliary_data(data_type="Configs", data=np.array([True]),
                              path='{}', parameters=par_dict)


