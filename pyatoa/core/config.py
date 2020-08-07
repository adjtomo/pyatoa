#!/usr/bin/env python3
"""
Configuration class that controls User-set parameters within the package.
Contains non-class functions for setting Config objects of Pyflex and Pyadjoint.
"""
import yaml
from pyatoa import logger
from pyatoa.utils.form import format_iter, format_step
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
    def __init__(self, yaml_fid=None, ds=None, path=None, iteration=None, 
                 step_count=None, event_id=None, min_period=10, max_period=30, 
                 filter_corners=2, client=None, rotate_to_rtz=False, 
                 unit_output="DISP", pyflex_preset="default", 
                 component_list=None, adj_src_type="cc_traveltime_misfit", 
                 start_pad=20, end_pad=500, observed_tag="observed", 
                 synthetic_tag=None, synthetics_only=False, win_amp_ratio=0., 
                 cfgpaths=None, save_to_ds=True, **kwargs):
        """
        Allows the user to control the parameters of the workflow, including:
        setting the Config objects for Pyflex and Pyadjoint, and the paths for
        input data, and output figures and results from Pyflex and Pyadjoint.
        Reasonable default values set on initation

        Kwargs are passed to Pyflex and Pyadjoint config objects so those can be
        set by the User through this Config object

        :type yaml_fid: str
        :param yaml_fid: id for .yaml file if config is to be loaded externally
        :type iteration: int
        :param iteration: if running an inversion, the current iteration. Used 
            for internal path naming, as well as interaction with Seisflows via
            Pyaflowa.
        :type step_count: int
        :param step: if running an inversion, the current step count in the
            line search, will be used for internal path naming, and interaction
            with Seisflows via Pyaflowa.
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
        self.iteration = iteration
        self.step_count = step_count
        self.event_id = event_id
        self.min_period = float(min_period)
        self.max_period = float(max_period)
        self.filter_corners = float(filter_corners)
        self.client = client
        self.rotate_to_rtz = rotate_to_rtz
        self.unit_output = unit_output.upper()
        self.observed_tag = observed_tag
        
        # Allow manual override of synthetic tag, but keep internal and rely 
        # on calling property for actual value
        self._synthetic_tag = synthetic_tag

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
        str_out = ("CONFIG\n"
                   f"    {'iteration:':<25}{self.iter_tag}\n"
                   f"    {'step_count:':<25}{self.step_tag}\n"
                   f"    {'event_id:':<25}{self.event_id}\n"
                   )
        # Format the remainder of the keys identically
        key_dict = {"Gather": ["client", "start_pad", "end_pad", "save_to_ds"],
                    "Process": ["min_period", "max_period", "filter_corners",
                                "unit_output", "rotate_to_rtz", "win_amp_ratio",
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
    def iter_tag(self):
        """string formatted version of iteration, e.g. 'i00'"""
        if self.iteration is not None:
            return format_iter(self.iteration)
        else:
            return None

    @property
    def step_tag(self):
        """string formatted version of step, e.g. 's00'"""
        if self.step_count is not None:
            return format_step(self.step_count)
        else:
            return None

    @property
    def synthetic_tag(self):
        """tag to be used for synthetic data, uses iteration and step count"""
        if self._synthetic_tag is not None:
            return self._synthetic_tag

        # If no override value given, fall back to default
        tag = self._get_aux_path(default=None, separator='')
        if tag is not None:
            return f"synthetic_{tag}"
        else:
            return "synthetic"

    @property
    def aux_path(self):
        """property to quickly get a bog-standard aux path e.g. i00/s00"""
        return self._get_aux_path()

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
        if unused_kwargs:
            raise ValueError(f"{unused_kwargs} are not keyword arguments in "
                             f"Pyatoa, Pyflex or Pyadjoint.")

    def _get_aux_path(self, default="default", separator="/"):
        """
        Pre-formatted path to be used for tagging and identification in 
        ASDF dataset auxiliary data. Internal function to be called by property
        aux_path.

        :type default: str
        :param default: if no iteration or step information is given, path will
            default to this string. By default it is 'default'.
        :type separator: str
        :param separator: if an iteration and step_count are available, 
            separator will be placed between. Defaults to '/', use '' for no
            separator.
        """
        if (self.iter_tag is not None) and self.step_tag is not None:
            # model/step/window_tag
            path = separator.join([self.iter_tag, self.step_tag])
        elif self.iter_tag is not None:
            path = self.iter_tag
        else:
            path = default

        return path

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
        from os.path import splitext

        # Ensure file ending
        if splitext(filename)[1] != ".yaml":
            filename += ".yaml"
        with open(filename, "w") as f:
            yaml.dump(vars(self), f, default_flow_style=False, sort_keys=False)

    def _write_asdf(self, ds):
        """
        Save the Config values as a parameter dictionary in the ASDF Data set
        Converts types to play nice with ASDF Auxiliary Data.
        Flattens dictionaries and external Config objects for easy storage.

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset to save the config file to
        """
        # Lazy imports because this function isn't always called
        from numpy import array
        from copy import deepcopy

        # Deep copy to ensure that we aren't editing the Config parameters
        attrs = vars(deepcopy(self))
        
        add_attrs = {}
        del_attrs = []
        for key, item in attrs.items():
            if item is None:
                # HDF doesn't support NoneType so convert to string        
                attrs[key] = "None"
            elif isinstance(item, (dict, PyflexConfig, PyadjointConfig)):
                # Flatten dictionaries, add prefix, delete original
                try:
                    # Config objects will need to be converted to dictionaries
                    vars_ = vars(item)
                except TypeError:
                    vars_ = item
                # Prepend a prefix for easier read-back, also convert NoneTypes
                vars_ = {f"{key}_{k}": ('' if i is None else i)
                         for k, i in vars_.items()
                         }
                del_attrs.append(key)
                add_attrs.update(vars_)

        # Update the dictionary after the fact
        for key in del_attrs:
            attrs.pop(key)
        attrs.update(add_attrs)

        ds.add_auxiliary_data(data_type="Configs", data=array([True]),
                              path=self.aux_path, parameters=attrs
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

        # Parameters from flattened dictionaries will need special treatment
        cfgpaths, pyflex_config, pyadjoint_config = {}, {}, {}

        for key, item in cfgin.items():
            # Convert the item into expected native Python objects
            if isinstance(item, str):
                item = None if item == "None" else item
            else:
                try:
                    item = item.item()
                except ValueError:
                    item = item.tolist()

            # Put the item in the correct dictionary
            if "cfgpaths" in key:
                # e.g. cfgpaths_waveforms -> waveforms
                cfgpaths[key.split('_')[1]] = item
            elif "pyflex_config" in key:
                pyflex_config["_".join(key.split('_')[2:])] = item
            elif "pyadjoint_config" in key:
                # e.g. pyadjoint_config_dlna_sigma_min -> dlna_sigma_min
                pyadjoint_config["_".join(key.split('_')[2:])] = item
            else:
                # Normal Config attribute
                setattr(self, key, item)

        # Assign the flattened dictionaries back into nested dictionaries
        setattr(self, "cfgpaths", cfgpaths)

        pyflex_config, _ = set_pyflex_config(**pyflex_config, choice=None)
        setattr(self, "pyflex_config", pyflex_config)

        pyadjoint_config, _ = set_pyadjoint_config(**pyadjoint_config)
        setattr(self, "pyadjoint_config", pyadjoint_config)


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
