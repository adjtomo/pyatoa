#!/usr/bin/env python3
"""
Configuration object for Pyatoa.

Fed into the processor class for workflow management, and also used for
information sharing between objects and functions.
"""
import pyatoa.utils.configurations.external as extcfg


class Config:
    """
    Configuration class that controls functionalities inside pyatoa
    """
    def __init__(self, model_number=None, event_id=None, min_period=10,
                 max_period=30, filter_corners=4, rotate_to_rtz=False,
                 unit_output='DISP', pyflex_config='default',
                 adj_src_type='multitaper_misfit', start_pad=20, end_pad=500,
                 zero_pad=20, synthetic_unit="DISP", observed_tag='observed',
                 synthetic_tag='synthetic_{model_num}',
                 map_corners={'lat_min': -42.5007, 'lat_max': -36.9488,
                              'lon_min': 172.9998, 'lon_max': 179.5077},
                 paths={'synthetics': [], 'waveforms': [], 'responses': [],
                        'auxiliary_data': []}):
        """
        Allows the user to control the parameters of the packages called within
        pyatoa, as well as control where the outputs (i.e. pyasdf and plots) are
        sent after processing occurs

        Reasonable default values set on initation

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
        :type pyflex_config: str
        :param pyflex_config: name of pyflex preset config to use
        :type adj_src_type: str
        :param adj_src_type: method of misfit quantification for Pyadjoint
        :type start_pad: int
        :param start_pad: seconds before event origintime to grab waveform data
            for use by data gathering class
        :type end_pad: int
        :param end_pad: seconds after event origintime to grab waveform data
        :type zero_pad: int
        :type zero_pad: seconds to zero-pad data front and back, used by the
            preprocess functions, useful for very small source-receiver
            distances where there may not be much time from origin time
            to first arrival
        :type synthetic_unit: str
        :param synthetic_unit: units of Specfem synthetics, 'DISP', 'VEL', 'ACC'
        :type observed_tag: str
        :param observed_tag: Tag to use for asdf dataset to label and search
            for obspy streams of observation data. Defaults 'observed'
        :type synthetic_tag: str
        :param synthetic_tag: Tag to use for asdf dataset to label and search
            for obspy streams of synthetic data. Defaults 'synthetic_{model_num}
            Tag must be formatted before use
        :type paths: dict of str
        :param paths: any absolute paths for Pyatoa to search for
            waveforms in. If path does not exist, it will automatically be
            skipped. Allows for work on multiple machines, by giving multiple
            paths for the same set of data, without needing to change config.
            Waveforms must be saved in a specific directory structure with a
            specific naming scheme



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
            # Format the model number to the way Pyatoa expects it
            if isinstance(model_number, str):
                # If e.g. model_number = "0"
                if not model_number[0] == "m":
                    self.model_number = 'm{:0>2}'.format(model_number)
                # If e.g. model_number = "m00"
                else:
                    self.model_number = model_number
            # If e.g. model_number = 0
            elif isinstance(model_number, int):
                self.model_number = 'm{:0>2}'.format(model_number)
        else:
            self.model_number = None

        self.event_id = event_id
        self.min_period = float(min_period)
        self.max_period = float(max_period)
        self.filter_corners = float(filter_corners)
        self.rotate_to_rtz = rotate_to_rtz
        self.unit_output = unit_output.upper()
        self.synthetic_unit = synthetic_unit.upper()
        self.observed_tag = observed_tag
        self.synthetic_tag = synthetic_tag.format(model_num=self.model_number)
        self.pyflex_config = (pyflex_config, None)
        self.pyadjoint_config = (adj_src_type, None)
        self.map_corners = map_corners
        self.paths = paths
        self.zero_pad = int(zero_pad)
        self.start_pad = int(start_pad)
        self.end_pad = int(end_pad)

        # Run internal functions to check the Config object
        self.component_list = ['Z', 'N', 'E']
        self._check_config()

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
                         pyflex=self.pyflex_config[0],
                         adjoint=self.pyadjoint_config[0],
                         paths_to_wavs=self.paths['waveforms'],
                         paths_to_syns=self.paths['synthetics'],
                         paths_to_resp=self.paths['responses']
                         )

    def _check_config(self):
        """
        Just make sure that some of the configuration parameters are set proper
        """
        # Check period range is acceptable
        assert(self.min_period < self.max_period)

        # Check that the map corners is a dict and contains proper keys
        acceptable_keys = ['lat_min', 'lat_max', 'lon_min', 'lon_max']
        assert(isinstance(self.map_corners, dict)), "map_corners should be dict"
        for key in self.map_corners.keys():
            assert(key in acceptable_keys), "key should be in {}".format(
                acceptable_keys)

        # Check if unit output properly set
        acceptable_units = ['DISP', 'VEL', 'ACC']
        assert(self.unit_output in acceptable_units), \
            "unit_output should be in {}".format(acceptable_units)

        assert(self.synthetic_unit in acceptable_units), \
            "synthetic_unit should be in {}".format(acceptable_units)

        assert(self.pyflex_config[0] in extcfg.pyflex_configs().keys()), \
            "pyflex_config should be in {}".format(
                extcfg.pyflex_config().keys())

        # Check that paths are in the proper format
        acceptable_keys = [
            'synthetics', 'waveforms', 'responses', 'auxiliary_data']
        assert(isinstance(self.paths, dict)), "paths should be a dict"
        for key in self.paths.keys():
            assert(key in acceptable_keys), \
                "path keys can only be in {}".format(acceptable_keys)

        # Rotate component list if necessary
        if self.rotate_to_rtz:
            self.component_list = ['Z', 'R', 'T']

        # If all the assertions pass, set a few behind-the-scenes settings
        # Set Pyflex config as a tuple, (name, pyflex.Config)
        self.pyflex_config = (self.pyflex_config[0], 
                              extcfg.set_pyflex_config(self)
                              )
        # Set Pyadjoint Config as a tuple, (adj source type, pyadjoint.Config)
        self.pyadjoint_config = (
            self.pyadjoint_config[0],
            extcfg.get_pyadjoint_config(choice=self.pyadjoint_config[0],
                                        min_period=self.min_period,
                                        max_period=self.max_period)
                                 )

    def write_to_txt(self, filename="./pyatoa_config.txt"):
        """
        write out config file to text file
        :type filename: str
        :param filename: filename to save config to
        """
        with open(filename.format(model_num=self.model_number), "w") as f:
            f.write("PYATOA CONFIGURATION FILE\n\n")
            f.write("Model Number:            {model_number}\n"
                    "Event ID:                {event_id}\n\n"
                    "PROCESSING\n"
                    "\tMinimum Filter Period: {min_period}\n"
                    "\tMaximum Filter Period: {max_period}\n"
                    "\tZero Pad:              {zero_pad}s\n"
                    "\tStart Pad:             {start_pad}s\n"
                    "\tEnd Pad:               {end_pad}s\n"
                    "\tFilter Corners:        {corner}\n"
                    "\tRotate to RTZ:         {rotate}\n"
                    "\tUnit Output:           {output}\n"
                    "\tSynthetic Unit:        {synunit}\n"
                    "ASDF Dataset\n"
                    "\tObserved Tag:          {obstag}\n"
                    "\tSynthetic Tag:         {syntag}\n"
                    "AUX. CONFIGS\n"
                    "\tPyflex Config:         {pyflex}\n"
                    "\tAdjoint Source Type:   {adjoint}\n"
                    "MISC\n"
                    "\tMap Corners            {mapcorners}\n"
                    "\tPaths to waveforms:    {paths_to_wavs}\n"
                    "\tPaths to synthetics:   {paths_to_syns}\n"
                    "\tPaths to responses:    {paths_to_resp}".format(
                        model_number=self.model_number, event_id=self.event_id,
                        min_period=self.min_period, max_period=self.max_period,
                        zero_pad=self.zero_pad, start_pad=self.start_pad,
                        end_pad=self.end_pad, corner=self.filter_corners,
                        rotate=self.rotate_to_rtz, output=self.unit_output,
                        synunit=self.synthetic_unit, obstag=self.observed_tag,
                        syntag=self.synthetic_tag, pyflex=self.pyflex_config[0],
                        adjoint=self.pyadjoint_config[0],
                        mapcorners=self.map_corners,
                        paths_to_wavs=self.paths['waveforms'],
                        paths_to_syns=self.paths['synthetics'],
                        paths_to_resp=self.paths['responses'])
            )

    def write_to_asdf(self, ds):
        """
        Save the config values as a dictionary in the pyasdf data format
        for easy lookback
        """
        # Lazy imports because this function isn't always called
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
                    "pyflex_config": self.pyflex_config[0],
                    "pyadjoint_config": self.pyadjoint_config[0]
                    }

        ds.add_auxiliary_data(data_type="Configs", data=np.array([True]),
                              path="{}".format(self.model_number),
                              parameters=par_dict)

