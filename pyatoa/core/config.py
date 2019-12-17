#!/usr/bin/env python3
"""
Configuration object for Pyatoa.

The config class is the main interaction object between the User and workflow.
Fed into the processor class for workflow management, and also used for
information sharing between objects and functions.
"""
import yaml
from pyatoa.plugins.pyflex_config import set_pyflex_config
from pyatoa.plugins.pyadjoint_config import set_pyadjoint_config


class Config:
    """
    Configuration class that controls functionalities inside Pyatoa
    """
    def __init__(self, yaml_fid=None, model_number=None, event_id=None,
                 min_period=10, max_period=30, filter_corners=4,
                 rotate_to_rtz=False, unit_output='DISP', pyflex_map='default',
                 component_list=None, adj_src_type='cc_traveltime_misfit',
                 start_pad=20, end_pad=500, zero_pad=0, synthetic_unit="DISP",
                 observed_tag='observed', synthetic_tag='synthetic_{model_num}',
                 synthetics_only=False, window_amplitude_ratio=0.,
                 map_corners=None, cfgpaths=None):
        """
        Allows the user to control the parameters of the workflow, including:
        setting the Config objects for Pyflex and Pyadjoint, and the paths for
        input data, and output figures and results from Pyflex and Pyadjoint.
        Reasonable default values set on initation

        :type yaml_fid: str
        :param yaml_fid: id for .yaml file if config is to be loaded externally
        :type model_number: int or str
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
        :type pyflex_map: str
        :param pyflex_map: name to map to pyflex preset config
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
        :type synthetics_only: bool
        :param synthetics_only: If the user is doing a synthetic-synthetic
            example, e.g. in a checkerboard test, this will tell the internal
            fetcher to search for observation data in the 'waveforms' path
            in the same manner that it searches for synthetic data. Also changes
            tags on waveform plots so it's obvious that the 'data' is synthetic.
        :type observed_tag: str
        :param observed_tag: Tag to use for asdf dataset to label and search
            for obspy streams of observation data. Defaults 'observed'
        :type synthetic_tag: str
        :param synthetic_tag: Tag to use for asdf dataset to label and search
            for obspy streams of synthetic data. Default 'synthetic_{model_num}'
            Tag must be formatted before use.
        :type cfgpaths: dict of str
        :param cfgpaths: any absolute paths for Pyatoa to search for
            waveforms in. If path does not exist, it will automatically be
            skipped. Allows for work on multiple machines, by giving multiple
            paths for the same set of data, without needing to change config.
            Waveforms must be saved in a specific directory structure with a
            specific naming scheme
        """
        # Load config parameters from .yaml
        if yaml_fid:
            self._read_yaml(yaml_fid)
            self._check()
            return

        if model_number is not None:
            # Format the model number to the way Pyatoa expects it
            if isinstance(model_number, str):
                # If e.g. model_number = "0"
                if not model_number[0] == "m":
                    self.model_number = f"m{model_number:0>2}"
                # If e.g. model_number = "m00"
                else:
                    self.model_number = model_number
            # If e.g. model_number = 0
            elif isinstance(model_number, int):
                self.model_number = f"m{model_number:0>2}"
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
        if model_number:
            self.synthetic_tag = synthetic_tag.format(
                model_num=self.model_number)
        else:
            self.synthetic_tag = "synthetic"
        self.pyflex_map = pyflex_map
        self.adj_src_type = adj_src_type
        self.map_corners = map_corners
        self.synthetics_only = synthetics_only
        self.window_amplitude_ratio = window_amplitude_ratio
        self.zero_pad = int(zero_pad)
        self.start_pad = int(start_pad)
        self.end_pad = int(end_pad)
        self.component_list = component_list or ['Z', 'N', 'E'] 

        # Make sure User provided paths are list objects as they will be looped
        # on during the workflow
        if cfgpaths:
            for key in cfgpaths:
                if not isinstance(cfgpaths[key], list):
                    cfgpaths[key] = [cfgpaths[key]]
            self.cfgpaths = cfgpaths
        else:
            self.cfgpaths = {"waveforms": [], "synthetics": [], "responses": []}

        # Run internal functions to check the Config object
        self.pyflex_config = None
        self.pyadjoint_config = None
        self._check()

    def __str__(self):
        """
        String representation of class Config for print statements
        """
        str_out = "CONFIG\n"
        for key, item in vars(self).items():
            str_out += f"\t{key+':':<25}{item}\n"
        return str_out

    def _check(self):
        """
        A check to make sure that the configuration parameters are set properly
        to avoid any problems throughout the workflow.
        """
        # Check period range is acceptable
        assert(self.min_period < self.max_period)

        # Check that the map corners is a dict and contains proper keys
        # else, set to default map corners for New Zealand North Island
        if self.map_corners is not None:
            assert(isinstance(self.map_corners, dict)), \
                "map_corners must be a dictionary object"
            acceptable_keys = ['lat_min', 'lat_max', 'lon_min', 'lon_max']
            for key in self.map_corners.keys():
                assert(key in acceptable_keys), "key should be in {}".format(
                    acceptable_keys)
        else:
            self.map_corners = {'lat_min': -42.5007, 'lat_max': -36.9488,
                                'lon_min': 172.9998, 'lon_max': 179.5077
                                }

        # Check if unit output properly set
        acceptable_units = ['DISP', 'VEL', 'ACC']
        assert(self.unit_output in acceptable_units), \
            "unit_output should be in {acceptable_units}"

        assert(self.synthetic_unit in acceptable_units), \
            "synthetic_unit should be in {acceptable_units}"

        # Check that paths are in the proper format
        required_keys = ['synthetics', 'waveforms', 'responses']
        assert(isinstance(self.cfgpaths, dict)), "paths should be a dict"
        for key in self.cfgpaths.keys():
            assert(key in required_keys), \
                "path keys can only be in {required_keys}"
        # Make sure that all the required keys are given in the dictionary
        for key in required_keys:
            if key not in self.cfgpaths.keys():
                self.cfgpaths[key] = []

        # Rotate component list if necessary
        if self.rotate_to_rtz:
            self.component_list = ['Z', 'R', 'T']

        # Check that the amplitude ratio is a reasonable number
        if self.window_amplitude_ratio > 0:
            assert(self.window_amplitude_ratio < 1), \
                "window amplitude ratio should be < 1"

        # Set Pyflex confict through wrapper function
        self.pyflex_config = set_pyflex_config(choice=self.pyflex_map,
                                               min_period=self.min_period,
                                               max_period=self.max_period
                                               )

        # Set Pyadjoint Config
        self.pyadjoint_config = set_pyadjoint_config(
            choice=self.adj_src_type, min_period=self.min_period,
            max_period=self.max_period
        )

    def write(self, write_to, fmt=None):
        """
        Wrapper for write functions

        :type fmt: str
        :param fmt: format to save parameters to,
            (available: 'yaml', 'ascii', 'asdf')
        :type write_to: str or pyasdf.ASDFDataSet
        :param write_to: filename to save config to, or dataset to save to
        """
        # If no format given, try to guess the format
        acceptable_formats = ["yaml", "asdf", "ascii"]
        if fmt not in acceptable_formats:
            from pyasdf.asdf_data_set import ASDFDataSet

            if isinstance(write_to, str):
                if ("yaml" or "yml") in write_to:
                    fmt = "yaml"
                elif ("txt" or "ascii") in write_to:
                    fmt = "ascii"
                else:
                    print(f"format must be given in {acceptable_formats}")
                    return
            elif isinstance(write_to, ASDFDataSet):
                fmt = "asdf"
            else:
                print(f"format must be given in {acceptable_formats}")
                return

        if fmt.lower() == "ascii":
            self._write_ascii(write_to)
        elif fmt.lower() == "yaml":
            self._write_yaml(write_to)
        elif fmt.lower() == "asdf":
            self._write_asdf(write_to)

    def read(self, read_from, path=None, fmt=None):
        """
        Wrapper for read functions

        :type read_from: str or pyasdf.ASDFDataSet
        :param read_from: filename to read config from, or ds to read from
        :type path: str
        :param path: if fmt='asdf', path to the config in the aux data
        :type fmt: str
        :param fmt: file format to read parameters from, will be guessed but
            can also be explicitely set (available: 'yaml', 'ascii', 'asdf')
        """
        # If no format given, try to guess the format
        acceptable_formats = ["yaml", "asdf"]
        if fmt not in acceptable_formats:
            from pyasdf.asdf_data_set import ASDFDataSet

            if isinstance(read_from, str):
                if ("yaml" or "yml") in read_from:
                    fmt = "yaml"
                elif ("txt" or "ascii") in read_from:
                    fmt = "ascii"
                else:
                    print(f"format must be given in {acceptable_formats}")
                    return
            elif isinstance(read_from, ASDFDataSet):
                fmt = "asdf"
            else:
                print(f"format must be given in {acceptable_formats}")
                return

        if fmt.lower() == "yaml":
            self._read_yaml(read_from)
        elif fmt.lower() == "asdf":
            assert(path is not None), "path must be defined"
            self._read_asdf(read_from, path=path)

    def _write_yaml(self, filename):
        """
        Write config parameters to a yaml file, retain order

        :type filename: str
        :param filename: filename to save yaml file
        """
        # Ensure file ending
        if filename[-5:] != ".yaml":
            filename += ".yaml"
        with open(filename, "w") as f:
            yaml.dump(vars(self), f, default_flow_style=False, sort_keys=False)

    def _write_asdf(self, ds):
        """
        Save the config values as a dictionary in the pyasdf data format
        for easy lookback

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to save the config file to
        """
        # Lazy imports because this function isn't always called
        from numpy import array
        from obspy import UTCDateTime

        # Add/standardize some variables before passing to dataset
        attrs = vars(self)

        # Auxiliary data doesn't like Nonetype objects
        for key, item in attrs.items():
            if item is None:
                attrs[key] = ''

        attrs["creation_time"] = str(UTCDateTime())

        # Auxiliary data can't take dictionaries so convert to lists, variables
        attrs["map_corners"] = [self.map_corners['lat_min'],
                                self.map_corners['lat_max'],
                                self.map_corners['lon_min'],
                                self.map_corners['lon_max']
                                ]

        attrs["cfgpaths_waveforms"] = self.cfgpaths["waveforms"]
        attrs["cfgpaths_synthetics"] = self.cfgpaths["synthetics"]
        attrs["cfgpaths_responses"] = self.cfgpaths["responses"]

        # remove the Config objects because pyASDF doesn't recognize them
        del attrs["pyflex_config"]
        del attrs["pyadjoint_config"]
        del attrs["cfgpaths"]

        ds.add_auxiliary_data(data_type="Configs", data=array([True]),
                              path=self.model_number or "default",
                              parameters=attrs
                              )

    def _write_ascii(self, filename):
        """
        Write the config parameters to an ascii file

        :type filename: str
        :param filename: filename to write the ascii file to
        """
        attrs = vars(self)
        with open(filename, "w") as f:
            f.write("PYATOA CONFIGURATION FILE\n")
            for key_a, item_a in attrs.items():
                # Excludes writing the Pyflex and Pyadjoint Config classes, but
                # instead writes the parameters of those Configs separately
                try:
                    attrs_b = vars(item_a)
                    f.write(f"{key_a}\n")
                    for key_b, item_b in attrs_b.items():
                        f.write(f"\t{key_b}: {item_b}\n")
                except TypeError:
                    f.write(f"{key_a}: {item_a}\n")

    def _read_yaml(self, filename):
        """
        Read config parameters from a yaml file, parse to attributes

        :type filename: str
        :param filename: filename to save yaml file
        """
        with open(filename, "r") as f:
            attrs = yaml.load(f, Loader=yaml.Loader)
        for key, item in attrs.items():
            setattr(self, key, item)

    def _read_asdf(self, ds, path):
        """
        Read and set config parameters from an ASDF Dataset, assumes that all
        necessary parameters are located in the auxiliary data subgroup of the
        dataset, which will be the case if the write_to_asdf() function was used
        Assumes some things about the structure of the auxiliary data.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset with config parameter to read
        :type path: str
        :param path: model number e.g. 'm00' or 'default'
        """
        cfgin = ds.auxiliary_data.Configs[path].parameters
        cfgpaths = {}
        for key, item in cfgin.items():
            if "cfgpaths" in key:
                cfgpaths[key.split('_')[1]] = item.any() or []
            elif key == "map_corners":
                map_corners = {'lat_min': item[0].item(),
                               'lat_max': item[1].item(),
                               'lon_min': item[2].item(),
                               'lon_max': item[3].item()
                               }
                setattr(self, key, map_corners)
            else:
                # Convert numpy objects into native python objects to avoid
                # any confusion when reading from ASDF format
                try:
                    setattr(self, key, item.item())
                except ValueError:
                    setattr(self, key, item.tolist())
                except AttributeError:
                    setattr(self, key, item)
        setattr(self, "cfgpaths", cfgpaths)

        self._check()

