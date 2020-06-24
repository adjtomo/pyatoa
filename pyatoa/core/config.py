#!/usr/bin/env python3
"""
Configuration class that controls User-set parameters within the package.
Contains non-class functions for setting Config objects of Pyflex and Pyadjoint.
"""
import yaml
from pyatoa import logger
from pyatoa.utils.form import format_model_number, format_step_count
from pyatoa.core.pyaflowa import pyaflowa_kwargs
from pyatoa.plugins.pyflex_presets import pyflex_presets

from pyflex import Config as PyflexConfig
from pyadjoint import Config as PyadjointConfig


class Config:
    """
    The Config class is the main interaction object between the User and
    workflow. Fed into the Manager class for workflow management, and also used
    for information sharing between objects and functions. Has the ability to
    read from and write to external files in various formats, or explicitely
    defined through scripts or interactive shells.
    """
    def __init__(self, yaml_fid=None, ds=None, path=None, model=None, step=None,
                 event_id=None, min_period=10, max_period=30, filter_corners=2,
                 client=None, rotate_to_rtz=False, unit_output="DISP",
                 pyflex_preset="default", component_list=None,
                 adj_src_type="cc_traveltime_misfit", start_pad=20, end_pad=500,
                 synthetic_unit="DISP", observed_tag="observed",
                 synthetic_tag="synthetic", synthetics_only=False,
                 win_amp_ratio=0., cfgpaths=None, save_to_ds=True, **kwargs):
        """
        Allows the user to control the parameters of the workflow, including:
        setting the Config objects for Pyflex and Pyadjoint, and the paths for
        input data, and output figures and results from Pyflex and Pyadjoint.
        Reasonable default values set on initation

        Kwargs are passed to Pyflex and Pyadjoint config objects so those can be
        set by the User through this Config object

        :type yaml_fid: str
        :param yaml_fid: id for .yaml file if config is to be loaded externally
        :type model: int
        :param model: model number, will be formatted for use in tags 
        :type step: int
        :param step: step count, will be formatted for use in tags
        :type event_id: str
        :param event_id: unique event identifier for data gathering, annotations
        :type: min_period: float
        :param min_period: minimum bandpass filter period
        :type max_period: float
        :param max_period: maximum bandpass filter period
        :type filter_corners: int
        :param filter_corners: filter steepness for obspy filter
        :type client: str
        :param client: ObsPy FDSN Client to be used for data gathering.
        :type rotate_to_rtz: bool
        :param rotate_to_rtz: components from NEZ to RTZ
        :type unit_output: str
        :param unit_output: units of stream, to be fed into preprocessor for
            instrument response removal. Available: 'DISP', 'VEL', 'ACC'
        :type pyflex_preset: str
        :param pyflex_preset: name to map to pyflex preset config
        :type adj_src_type: str
        :param adj_src_type: method of misfit quantification for Pyadjoint
        :type start_pad: int
        :param start_pad: seconds before event origintime to grab waveform data
            for use by data gathering class
        :type end_pad: int
        :param end_pad: seconds after event origintime to grab waveform data
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
        :type save_to_ds: bool
        :param save_to_ds: allow toggling saving to the dataset when new data
            is gathered/collected. This is useful, e.g. if a dataset that
            contains data is passed to the Manager, but you don't want to
            overwrite the data inside while you do some temporary processing.
        """
        self.model = model
        self.step = step
        self.event_id = event_id
        self.min_period = float(min_period)
        self.max_period = float(max_period)
        self.filter_corners = float(filter_corners)
        self.client = client
        self.rotate_to_rtz = rotate_to_rtz
        self.unit_output = unit_output.upper()
        self.synthetic_unit = synthetic_unit.upper()
        self.observed_tag = observed_tag

        self.synthetic_tag = synthetic_tag
        if self.model:
            # Tag based on model number and step count, e.g. synthetic_m00s00
            self.synthetic_tag += (f"_{self.model_number}"
                                   f"{self.step_count or ''}"
                                   )

        self.pyflex_preset = pyflex_preset
        self.adj_src_type = adj_src_type
        self.synthetics_only = synthetics_only
        self.win_amp_ratio = win_amp_ratio
        self.start_pad = int(start_pad)
        self.end_pad = int(end_pad)
        self.component_list = component_list or ["Z", "N", "E"]
        self.save_to_ds = save_to_ds

        # These are filled in with actual Config objects by _check()
        self.pyflex_config = None
        self.pyadjoint_config = None

        # Make sure User provided paths are list objects
        if cfgpaths:
            for key in cfgpaths:
                if not isinstance(cfgpaths[key], list):
                    cfgpaths[key] = [cfgpaths[key]]
            self.cfgpaths = cfgpaths
        else:
            self.cfgpaths = {"waveforms": [], "synthetics": [], "responses": []}

        # Overwrite config parameters from .yaml if given
        if yaml_fid:
            kwargs = self._read_yaml(yaml_fid)
        elif ds:
            assert(path is not None), "'path' is required to load from dataset"
            self._read_asdf(ds, path=path)

        # Run internal sanity checks
        self._check(**kwargs)

    def __str__(self):
        """
        String representation of class Config for print statements.
        Separate into similar labels for easier reading.
        """
        # Model and step need to be formatted before printing
        str_out = ("Config\n"
                   f"    {'model:':<25}{self.model_number}\n"
                   f"    {'step:':<25}{self.step_count}\n"
                   f"    {'event:':<25}{self.event_id}\n"
                   )
        # Format the remainder of the keys identically
        key_dict = {"Gather": ["client", "start_pad", "end_pad", "save_to_ds"],
                    "Process": ["min_period", "max_period", "filter_corners",
                                "unit_output", "synthetic_unit",
                                "rotate_to_rtz", "win_amp_ratio",
                                "synthetics_only"],
                    "Labels": ["component_list", "observed_tag",
                               "synthetic_tag", "cfgpaths"],
                    "External": ["pyflex_preset", "adj_src_type",
                                 "pyflex_config", "pyadjoint_config"
                                 ]
                    }
        for key, items in key_dict.items():
            str_out += f"{key.upper()}\n"
            for item in items:
                str_out += f"    {item+':':<25}{getattr(self, item)}\n"
        return str_out

    def __repr__(self):
        """Simply call string representation"""
        return self.__str__()

    @property
    def model_number(self):
        """string formatted version of model, e.g. 'm00'"""
        if self.model is not None:
            return format_model_number(self.model)
        else:
            return None

    @property
    def step_count(self):
        """string formatted version of step, e.g. 's00'"""
        if self.step is not None:
            return format_step_count(self.step)
        else:
            return None

    def _check(self, **kwargs):
        """
        A series of sanity checks to make sure that the configuration parameters
        are set properly to avoid any problems throughout the workflow.
        """
        # Check period range is acceptable
        assert(self.min_period < self.max_period), \
            "min_period must be less than max_period"

        # Check if unit output properly set, dictated by ObsPy units
        acceptable_units = ['DISP', 'VEL', 'ACC']
        assert(self.unit_output in acceptable_units), \
            f"unit_output should be in {acceptable_units}"

        assert(self.synthetic_unit in acceptable_units), \
            f"synthetic_unit should be in {acceptable_units}"

        # Check that paths are in the proper format, dictated by Pyatoa
        required_keys = ['synthetics', 'waveforms', 'responses']
        assert(isinstance(self.cfgpaths, dict)), "paths should be a dict"
        for key in self.cfgpaths.keys():
            assert(key in required_keys), \
                f"path keys can only be in {required_keys}"
        # Make sure that all the required keys are given in the dictionary
        for key in required_keys:
            if key not in self.cfgpaths.keys():
                self.cfgpaths[key] = []

        # Rotate component list if necessary
        if self.rotate_to_rtz:
            logger.debug("Components changed ZNE -> ZRT")
            for comp in ["N", "E"]:
                assert(comp not in self.component_list), \
                        f"rotated component list cannot include '{comp}'"

        # Check that the amplitude ratio is a reasonable number
        if self.win_amp_ratio > 0:
            assert(self.win_amp_ratio < 1), \
                "window amplitude ratio should be < 1"

        # Make sure adjoint source type is formatted properly
        self.adj_src_type = format_adj_src_type(self.adj_src_type)

        # Set Pyflex confict through wrapper function
        self.pyflex_config, unused_kwargs_pf = set_pyflex_config(
            choice=self.pyflex_preset, min_period=self.min_period,
            max_period=self.max_period, **kwargs
        )

        # Set Pyadjoint Config
        self.pyadjoint_config, unused_kwargs_pa = set_pyadjoint_config(
            min_period=self.min_period, max_period=self.max_period, **kwargs
        )

        # Check for unnused kwargs
        unused_kwargs = []
        for kwarg in unused_kwargs_pf:
            if kwarg in unused_kwargs_pa and kwarg not in pyaflowa_kwargs:
                unused_kwargs.append(kwarg)
        if unused_kwargs:
            raise ValueError(f"{unused_kwargs} are not keyword arguments in "
                             f"Pyatoa, Pyflex or Pyadjoint.")

    @staticmethod
    def _check_io_format(fid, fmt=None):
        """
        A simple check before reading or writing the config to determine what
        file format to use.

        :type fmt: str
        :param fmt: format specified by the User
        :rtype: str
        :return: format string to be understood by the calling function
        """
        acceptable_formats = ["yaml", "asdf", "ascii"]
        if fmt not in acceptable_formats:
            # If no format given, try to guess the format based on file ending
            from pyasdf.asdf_data_set import ASDFDataSet
            if isinstance(fid, str):
                if ("yaml" or "yml") in fid:
                    return "yaml"
                elif ("txt" or "ascii") in fid:
                    return "ascii"
                else:
                    raise TypeError(
                        "format must be given in {acceptable_formats}")
            elif isinstance(fid, ASDFDataSet):
                return "asdf"
            else:
                raise TypeError("file must be given in {acceptable_formats}")
        else:
            return fmt

    def write(self, write_to, fmt=None):
        """
        Wrapper for write functions

        :type fmt: str
        :param fmt: format to save parameters to,
            (available: 'yaml', 'ascii', 'asdf')
        :type write_to: str or pyasdf.ASDFDataSet
        :param write_to: filename to save config to, or dataset to save to
        """
        fmt = self._check_io_format(write_to, fmt)

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
        fmt = self._check_io_format(read_from, fmt)

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
        from copy import deepcopy

        # Add/standardize some variables before passing to dataset
        # deep copy to ensure that we aren't editing the Config parameters
        attrs = vars(deepcopy(self))

        # Auxiliary data doesn't like NoneType objects
        for key, item in attrs.items():
            if item is None:
                attrs[key] = ''

        attrs["creation_time"] = str(UTCDateTime())
        attrs["cfgpaths_waveforms"] = self.cfgpaths["waveforms"]
        attrs["cfgpaths_synthetics"] = self.cfgpaths["synthetics"]
        attrs["cfgpaths_responses"] = self.cfgpaths["responses"]

        # remove Config objects because pyASDF won't recognize dicts or objects
        # these will be reinstated by _check() if/when read back in
        del attrs["pyflex_config"]
        del attrs["pyadjoint_config"]
        del attrs["cfgpaths"]

        # Figure out how to tag the data in the dataset
        if self.model and self.step:
            # model/step/window_tag
            path = f"{self.model}/{self.step}"
        elif self.model:
            path = self.model
        else:
            path = "default"

        ds.add_auxiliary_data(data_type="Configs", data=array([True]),
                              path=path, parameters=attrs
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
        Read config parameters from a yaml file, parse to attributes.
        Any non-standard parameters will be passed through as kwargs.

        :type filename: str
        :param filename: filename to save yaml file
        :rtype: dict
        :return: key word arguments that do not belong to Pyatoa are passed back
            as a dictionary object, these are expected to be arguments that are
            to be used in Pyflex and Pyadjoint configs
        """
        with open(filename, "r") as f:
            attrs = yaml.load(f, Loader=yaml.Loader)
        # Check if we're reading from a Seisflows yaml file
        if 'PYATOA' in attrs.keys():
            attr_list = attrs['PYATOA'].items()
        else:
            attr_list = attrs.items()
        
        kwargs = {}
        for key, item in attr_list:
            if hasattr(self, key.lower()):
                # Special case: ensure cfgpaths don't overwrite, but append
                if key == "cfgpaths":
                    for cfgkey, cfgitem in self.cfgpaths.items():
                        item[cfgkey] += cfgitem
                setattr(self, key.lower(), item)
            else:
                kwargs[key.lower()] = item

        return kwargs

    def _read_asdf(self, ds, path):
        """
        Read and set config parameters from an ASDF Dataset, assumes that all
        necessary parameters are located in the auxiliary data subgroup of the
        dataset, which will be the case if the write_to_asdf() function was used
        Assumes some things about the structure of the auxiliary data.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset with config parameter to read
        :type path: str
        :param path: model number e.g. 'm00' or 'default', or 'm00/s00'
        """
        # Check if nested paths are provided
        splitpath = path.split("/")
        if len(splitpath) > 1:
            cfgin = ds.auxiliary_data.Configs
            for p in splitpath:
                cfgin = cfgin[p]
            cfgin = cfgin.parameters
        else:
            cfgin = ds.auxiliary_data.Configs[path].parameters

        cfgpaths = {}
        for key, item in cfgin.items():
            if "cfgpaths" in key:
                cfgpaths[key.split('_')[1]] = item.any() or []
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


def set_pyflex_config(min_period, max_period, choice=None, **kwargs):
    """
    Overwriting the default Pyflex parameters with User-defined criteria

    :type choice: str or dict
    :param choice: name of map to choose the Pyflex config options, if None,
        default values are used. Kwargs can still overload default values.
        Also dicts can be passed in as User-defined preset
    :type min_period: float
    :param min_period: min period of the data
    :type max_period: float
    :param max_period: max period of the data
    :rtype: pyflex.Config
    :return: the pyflex Config option to use when running Pyflex
    """
    # Instantiate the pyflex Config object
    pfconfig = PyflexConfig(min_period=min_period, max_period=max_period)

    # Set preset configuration parameters based on hard-coded presets
    if isinstance(choice, str):
        if choice in pyflex_presets.keys():
            preset = pyflex_presets[choice]
            for key, item in preset.items():
                setattr(pfconfig, key, item)
        else:
            raise KeyError(
                f"'{choice}' does not match any available presets for Pyflex. "
                f"Presets include {list(pyflex_presets.keys())}"
            )
    # Allow dictionary object to be passed in as a preset
    elif isinstance(choice, dict):
        for key, item in choice.items():
            setattr(pfconfig, key, item)

    # Kwargs can also be passed from the pyatoa.Config object to avoid having to
    # define pre-set values. Kwargs will override preset values
    unused_kwargs = []
    for key, item in kwargs.items():
        if hasattr(pfconfig, key):
            setattr(pfconfig, key, item)
        else:
            unused_kwargs.append(key)

    return pfconfig, unused_kwargs


def set_pyadjoint_config(min_period, max_period, **kwargs):
    """
    Set the Pyadjoint config based on Pyatoa Config parameters.
    Kwargs can be fed to the Pyadjoint Config object. Returns unnused kwargs.

    Config parameters can be found at:
    https://github.com/krischer/pyadjoint/blob/master/src/pyadjoint/config.py

    :type min_period: float
    :param min_period: min period of the data
    :type max_period: float
    :param max_period: max period of the data
    :rtype cfgout: pyadjoint.Config
    :return cfgout: properly set pyadjoint configuration object
    """
    paconfig = PyadjointConfig(min_period=min_period,
                               max_period=max_period
                               )
    unused_kwargs = []
    for key, item in kwargs.items():
        if hasattr(paconfig, key):
            setattr(paconfig, key, item)
        else:
            unused_kwargs.append(key)

    return paconfig, unused_kwargs


def format_adj_src_type(choice):
    """
    Pyadjoint requires that a standard adjoint source type is given in the
    calculate_adjoint_source() function. This function acts as simple dictionary
    input to provide the correct input for that function. Allows for various
    spellings and variations of the same name.

    Throws ValueError if choice isn't in the provided lists.

    :type choice: str
    :param choice: pyatoa.Config.adj_src_type
    :rtype: str
    :return: pyadjoint adj_src_type
    """
    if choice in ["cc", "cc_traveltime_misfit", "cross_correlation"]:
        adj_src_type = "cc_traveltime_misfit"
    elif choice in ["mt", "mtm", "multitaper_misfit", "multitaper"]:
        adj_src_type = "multitaper_misfit"
    elif choice in ["wav", "wave", "waveform", "w"]:
        adj_src_type = "waveform"
    else:
        raise ValueError(f"'{choice}' does not match available adjoint source "
                         f"types, must be 'cc', 'mt', or 'wav'")
    return adj_src_type
