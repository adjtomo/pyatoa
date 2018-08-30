#!/usr/bin/env python3
"""
Configuration for pyatoa runs
"""
import os
import warnings


class Config():
    def __init__(self, model_number=None, event_id=None, min_period=10,
                 max_period=30, filter_corners=4, rotate_to_rtz=True,
                 unit_output='DISP', pyflex_config='UAF',
                 adj_src_type='multitaper_misfit', paths_to_waveforms=None,
                 path_to_pyasdf=None,
                 verbose=False, log=False):
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
            min_period=10,max_period=30,rotate_to_rtz=True,unit_output="DISP",
            pyflex_config='UAF',adj_src_type='multitaper_misfit',verbose=True,
            log=True,
            paths_to_waveforms=['/Users/chowbr/Documents/subduction/seismic/'],
            paths_to_synthetics=['/Users/chowbr/Documents/subduction/seismic/SPECFEM3D']
            )
        """
        if model_number is not None:
            self.model_number = 'm{:0>2}'.format(model_number)
        self.event_id = event_id
        self.min_period = float(min_period)
        self.max_period = float(max_period)
        self.filter_corners = float(filter_corners)
        self.rotate_to_rtz = rotate_to_rtz
        self.unit_output = unit_output.upper()
        self.pyflex_config = pyflex_config
        self._generate_pyflex_config()
        self.adj_src_type = adj_src_type
        self.verbose = verbose
        self.log = log
        self.paths = {"waveforms":paths_to_waveforms,
                      "pyasdf":path_to_pyasdf
                      }
        self._generate_component_list()

    def __str__(self):
        """
        string representation of class Config for print statements
        :return:
        """
        template = ("Config class\n"
                    "\tEvent ID:              {ei}\n"
                    "\tMinimum Filter Period: {f1}\n"
                    "\tMaximum Filter Period: {f2}\n"
                    "\tFilter Corners:        {fc}\n"
                    "\tRotate to RTZ:         {rr}\n"
                    "\tUnit Output:           {uo}\n"
                    "\tPyflex Config:         {pc}\n"
                    "\tAdjoint Source Type:   {at}\n"
                    "\tPaths to waveforms:    {p1}\n"
                    "\tPaths to synthetics:   {p2}\n"
                    "\tPath to PyASDF:        {p3}\n"
                    )
        return template.format(ei=self.event_id, f1=self.min_period,
                               f2=self.max_period, fc=self.filter_corners,
                               rr=self.rotate_to_rtz, uo=self.unit_output,
                               pc=self.pyflex_config, at=self.adj_src_type,
                               p1=self.paths['waveforms'],
                               p2=self.paths['synthetics'],
                               p3=self.paths['pyasdf']
                               )

    def _generate_component_list(self):
        """
        create a small list for easy access to orthogonal components
        """
        if self.rotate_to_rtz:
            self.component_list = ['Z', 'R', 'T']
        else:
            self.component_list = ['Z', 'N', 'E']

    def _generate_pyflex_config(self):
        """
        allow for dictionary lookup of pyflex config, or manual set
        :return:
        """
        cfg_dict = {"default": [.08, 15., 1., .8, .7, 4., 0., 1., 2., 3., 10.],
                    "UAF": [.18, 4., 1.5, .71, .7, 2., 0., 3., 2., 2.5, 12.]
                    }
        if type(self.pyflex_config) == str:
            self.pyflex_config = cfg_dict[self.pyflex_config]
        elif type(self.pyflex_config) == list and len(self.pyflex_config) != 11:
            warnings.warn("given 'pyflex_config' is not correct length",
                          UserWarning)
        elif type(self.pyflex_config) != list:
            warnings.warn("'pyflex_config' must be type 'str' or 'list",
                          UserWarning)

    def _check_config_parameters(self):
        """
        just make sure that some of the configuration parameters are set proper
        :return:
        """
        if self.paths['waveforms'] is not None:
            for path_ in self.paths.values():
                for p in path_:
                    if not os.path.exists(p):
                        warnings.warn("path '{}' does not exist".format(P2W),
                                      UserWarning)
        if self.unit_output not in ['DISP', 'VEL', 'ACC']:
            warnings.warn("'unit_output' must be 'DISP','VEL' or 'ACC'",
                          UserWarning)

    def write_to_txt_file(self,filename):
        """
        write out conf file to text file
        :param filename:
        :return:
        """

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
                    "rotate_to_rtz":self.rotate_to_rtz,
                    "unit_output":self.unit_output,
                    "pyflex_config":self.pyflex_config,
                    "adj_src_type":self.adj_src_type
                    }
        ds.add_auxiliary_data(data_type="Configs", data=np.array([True]),
                              path='{}', parameters=par_dict)


