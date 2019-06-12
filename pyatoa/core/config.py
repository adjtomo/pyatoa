#!/usr/bin/env python3
"""
Configuration object for Pyatoa.

Fed into the processor class for workflow management, and also used for
information sharing between objects and functions.
"""
import os
import warnings


class Config:
    """
    Configuration class that controls functionalities inside pyatoa
    """
    def __init__(self, model_number=None, event_id=None, min_period=10,
                 max_period=30, filter_corners=4, rotate_to_rtz=True,
                 unit_output='DISP', pyflex_config='default',
                 adj_src_type='multitaper_misfit', startpad=20, endpad=500,
                 paths_to_waveforms=[], paths_to_synthetics=[],
                 paths_to_responses=[]):
        """
        Allows the user to control the parameters of the packages called within
        pyatoa, as well as control where the outputs (i.e. pyasdf and plots) are
        sent after processing occurs

        :type model_number: int
        :param model_number: model iteration number for annotations and tags
        :type event_id: str
        :param event_id: unique event identifier for data gathering, annotations
        :type: min_period: float
        :param min_period: minimum bandpass filter period
        :type max_period: float
        :param max_period: maximum bandpass filter period
        :type filter_corners: int
        :param filter_corners: filter steepness for obspy filter
        :type rotate_to_rtz: bool
        :param rotate_to_rtz: components from NEZ to RTZ
        :type unit_output: str
        :param unit_output: units of stream, to be fed into preprocessor for
            instrument response removal. Available: 'DISP', 'VEL', 'ACC'
        :type pyflex_config: (list of floats) or str
        :param pyflex_config: values to be fed into the pyflex config object.
            Can give dictionary key for presets: 'default', and 'UAF'
            or give a manual entry list of floats with the following format:
            i  Standard Tuning Parameters:
            0: water level for STA/LTA (short term average/long term average)
            1: time lag acceptance level
            2: amplitude ratio acceptance level (dlna)
            3: normalized cross correlation acceptance level
            i  Fine Tuning Parameters
            4: c_0 = for rejection of internal minima
            5: c_1 = for rejection of short windows
            6: c_2 = for rejection of un-prominent windows
            7: c_3a = for rejection of multiple distinct arrivals
            8: c_3b = for rejection of multiple distinct arrivals
            9: c_4a = for curtailing windows w/ emergent starts and/or codas
            10:c_4b = for curtailing windows w/ emergent starts and/or codas
        :type adj_src_type: str
        :param adj_src_type: method of misfit quantification specified by
            Pyadjoint.
            Available: 'waveform', 'cc_traveltime_misfit', 'multitaper_misfit'
            (http://krischer.github.io/pyadjoint/adjoint_sources/index.html)
        :type raw_sampling_rate: float
        :param raw_sampling_rate: waveform streams are saved in raw format
            for later access, but fetched data is usually very high sampling
            rate (e.g. 100Hz for HH? stations). This tag sets the raw sampling
            rate of the saved data, useful for saving space.
        :type paths_to_waveforms: list of str
        :param paths_to_waveforms: any absolute paths for Pyatoa to search for
            waveforms in. If path does not exist, it will automatically be
            skipped. Allows for work on multiple machines, by giving multiple
            paths for the same set of data, without needing to change config.
            Waveforms must be saved in a specific directory structure with a
            specific naming scheme
        :type paths_to_responses: list of str
        :param paths_to_responses: any absolute paths for Pyatoa to search for
            responses in.


        Example call:
        config = Config(model_number=0,event_id='2014p240655',
            min_period=10,max_period=30,rotate_to_rtz=True,unit_output="DISP",
            pyflex_config='UAF',adj_src_type='multitaper_misfit',verbose=True,
            log=True,
            paths_to_waveforms=['/Users/chowbr/Documents/subduction/seismic',
                '/seis/prj/fwi/bchow/seismic','/geonet/seismic'],
            )
        """
        if model_number is not None:
            if isinstance(model_number, str):
                # model_number = "0"
                if not model_number[0] == "m":
                    self.model_number = 'm{:0>2}'.format(model_number)
                # model_number = "m00"
                else:
                    self.model_number = model_number
            # model_number = 0
            elif isinstance(model_number, int):
                self.model_number = 'm{:0>2}'.format(model_number)
        else:
            self.model_number = model_number
        self.event_id = event_id
        self.min_period = float(min_period)
        self.max_period = float(max_period)
        self.filter_corners = float(filter_corners)
        self.rotate_to_rtz = rotate_to_rtz
        self.unit_output = unit_output.upper()
        self.pyflex_config = pyflex_config
        self._generate_pyflex_config()
        self.adj_src_type = adj_src_type
        self.paths = {"waveforms": paths_to_waveforms,
                      "synthetics": paths_to_synthetics,
                      "responses": paths_to_responses,
                      }
        self.startpad = startpad
        self.endpad = endpad
        self._generate_component_list()

    def __str__(self):
        """
        string representation of class Config for print statements
        :return:
        """
        return ("CONFIG\n"
                "\tModel Number:          {model_number}\n"
                "\tEvent ID:              {event_id}\n"
                "\tMinimum Filter Period: {min_period}\n"
                "\tMaximum Filter Period: {max_period}\n"
                "\tFilter Corners:        {corner}\n"
                "\tRotate to RTZ:         {rotate}\n"
                "\tUnit Output:           {output}\n"
                "\tPyflex Config:         {pyflex}\n"
                "\tAdjoint Source Type:   {adjoint}\n"
                "\tPaths to waveforms:    {paths_to_wavs}\n"
                "\tPaths to synthetics:   {paths_to_syns}\n"
                "\tPaths to responses:    {paths_to_resp}"
                ).format(model_number=self.model_number,
                         event_id=self.event_id, min_period=self.min_period,
                         max_period=self.max_period, corner=self.filter_corners,
                         rotate=self.rotate_to_rtz, output=self.unit_output,
                         pyflex=self.pyflex_config, adjoint=self.adj_src_type,
                         paths_to_wavs=self.paths['waveforms'],
                         paths_to_syns=self.paths['synthetics'],
                         paths_to_resp=self.paths['responses']
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
        Allow for dictionary lookup of pyflex config, or manual set
        Acceptance levels can also be set as arrays, defining a time-dependent
        acceptance level. This must be the same length as your waveform

        ++PYFLEX (Maggi et al. 2009)
        i  Standard Tuning Parameters:
        0: water level for STA/LTA (short term average/long term average)
        1: time lag acceptance level
        2: amplitude ratio acceptance level (dlna)
        3: normalized cross correlation acceptance level
        i  Fine Tuning Parameters
        4: c_0 = for rejection of internal minima
        5: c_1 = for rejection of short windows
        6: c_2 = for rejection of un-prominent windows
        7: c_3a = for rejection of multiple distinct arrivals
        8: c_3b = for rejection of multiple distinct arrivals
        9: c_4a = for curtailing windows w/ emergent starts and/or codas
        10:c_4b = for curtailing windows w/ emergent starts and/or codas
        """
        cfg_dict = {"default": [.08, 15., 1., .8, .7, 4., 0., 1., 2., 3., 10.],
                    "UAF": [.18, 4., 1.5, .71, .7, 2., 0., 3., 2., 2.5, 12.],
                    "NZ": [.15, 15., 1.75, .7, .7, 2., 0., 3., 2., 2.5, 12.],
                    }
        if isinstance(self.pyflex_config, str):
            self.pyflex_config = cfg_dict[self.pyflex_config]
        elif isinstance(self.pyflex_config, list) and \
                len(self.pyflex_config) != 11:
            warnings.warn("given 'pyflex_config' is not correct length",
                          UserWarning)
        elif not isinstance(self.pyflex_config, list):
            warnings.warn("'pyflex_config' must be type 'str' or 'list",
                          UserWarning)

    def _check_config_parameters(self):
        """
        just make sure that some of the configuration parameters are set proper
        :return:
        """
        if len(self.paths):
            for path_ in self.paths.values():
                for p in path_:
                    if not os.path.exists(p):
                        warnings.warn("path '{}' does not exist".format(p),
                                      UserWarning)
        if self.unit_output not in ['DISP', 'VEL', 'ACC']:
            warnings.warn("'unit_output' must be 'DISP','VEL' or 'ACC'",
                          UserWarning)

    def set_pyflex_config(self, config):
        """
        convenience function for pyflex testing to quickly change between
        different pyflex configs
        :param config:
        :return:
        """
        self.pyflex_config = config
        self._generate_pyflex_config()

    def write_to_txt_file(self, filename):
        """
        write out conf file to text file
        :param filename:
        :return:
        """
        raise NotImplementedError

    def write_to_asdf(self, ds):
        """
        save the config values as a dictionary in the pyasdf data format
        for easy lookback
        """
        import numpy as np
        from obspy import UTCDateTime
        par_dict = {"creation_time": str(UTCDateTime()),
                    "model_number": self.model_number,
                    "event_id": self.event_id,
                    "min_period": self.min_period,
                    "max_period": self.max_period,
                    "filter_corners": self.filter_corners,
                    "rotate_to_rtz": self.rotate_to_rtz,
                    "unit_output": self.unit_output,
                    "pyflex_config": self.pyflex_config,
                    "adj_src_type": self.adj_src_type,
                    }

        ds.add_auxiliary_data(data_type="Configs", data=np.array([True]),
                              path="{}".format(self.model_number),
                              parameters=par_dict)

