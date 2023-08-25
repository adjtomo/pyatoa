#!/usr/bin/env python3
"""
Configuration of User-set parameters within the package.
Contains external functions to set Config objects of Pyflex and Pyadjoint.
"""
import yaml
import numpy as np
from copy import deepcopy
from pyatoa.utils.form import format_iter, format_step

from pyflex import Config as PyflexConfig
from pyadjoint import get_config as get_pyadjoint_config
from pyadjoint import ADJSRC_TYPES


class Config:
    """
    The Config class is the main interaction object between the User and
    workflow. It is used by :class:`Manager <pyatoa.core.manager.Manager>` for
    workflow management, and also for information sharing between Pyatoa objects
    and functions. The Config can be read to and written from external files and
    ASDFDataSets.
    """
    def __init__(self, yaml_fid=None, ds=None, path=None, iteration=None,
                 step_count=None, event_id=None, min_period=10, max_period=100,
                 rotate_to_rtz=False, unit_output="DISP",  component_list=None,
                 adj_src_type="cc_traveltime", observed_tag="observed",
                 synthetic_tag=None, st_obs_type="obs", st_syn_type="syn",
                 win_amp_ratio=0., pyflex_parameters=None,
                 pyadjoint_parameters=None):
        """
        Initiate the Config object either from scratch, or read from external.

        .. note::
            keyword arguments are passed to Pyflex and Pyadjoint config objects
            so that there is only one entry point to all config objects.

        :type yaml_fid: str
        :param yaml_fid: id for .yaml file if config is to be loaded externally
        :type iteration: int
        :param iteration: if running an inversion, the current iteration. Used
            for internal path naming, as well as interaction with Seisflows via
            Pyaflowa.
        :type step_count: int
        :param step_count: if running an inversion, the current step count in the
            line search, will be used for internal path naming, and interaction
            with Seisflows via Pyaflowa.
        :type event_id: str
        :param event_id: unique event identifier for data gathering, annotations
        :type min_period: float
        :param min_period: minimum bandpass filter period
        :type max_period: float
        :param max_period: maximum bandpass filter period
        :type rotate_to_rtz: bool
        :param rotate_to_rtz: components from NEZ to RTZ
        :type unit_output: str
        :param unit_output: units of stream, to be fed into preprocessor for
            instrument response removal. Available: 'DISP', 'VEL', 'ACC'
        :type adj_src_type: str
        :param adj_src_type: method of misfit quantification for Pyadjoint
        :type st_obs_type: str
        :param st_obs_type: Tell Pyatoa how to treat `st_obs`, either
            - 'data': as data, which involves instrument response removal and
                data gathering based on SEED formatted directories
            - 'syn': as syntheitcs, which skips instrument response removal
                and data gathering is based on simpler synthetic dir. structure
            Defaults to 'data'
        :type st_syn_type: str
        :param st_syn_type: Tell Pyatoa how to treat `st_syn`, either
            - 'data': as data, which involves instrument response removal and
                data gathering based on SEED formatted directories
            - 'syn': as syntheitcs, which skips instrument response removal
                and data gathering is based on simpler synthetic dir. structure
            Defaults to 'syn'
        :type observed_tag: str
        :param observed_tag: Tag to use for asdf dataset to label and search
            for obspy streams of observation data. Defaults 'observed'
        :type synthetic_tag: str
        :param synthetic_tag: Tag to use for asdf dataset to label and search
            for obspy streams of synthetic data. Default 'synthetic_{model_num}'
            Tag must be formatted before use.
        :type pyflex_parameters: dict
        :param pyflex_parameters: overwrite for Pyflex parameters defined
            in the Pyflex.Config object. Incorrectly defined argument names
            will raise a TypeError. See Pyflex docs for detailed parameter defs:
            http://adjtomo.github.io/pyflex/#config-object
        :type pyadjoint_parameters: dict
        :param pyadjoint_parameters: overwrite for Pyadjoint parameters defined
            in the Pyadjoint.Config object for the given `adj_src_type`.
            Incorrectly defined argument names will raise a TypeError. See
            Pyadjoint docs for detailed parameter definitions:
            https://adjtomo.github.io/pyadjoint/
        :raises TypeError: If incorrect arguments provided to the underlying
            Pyflex or Pyadjoint Config objects.
        """
        self.iteration = iteration
        self.step_count = step_count
        self.event_id = event_id
        self.min_period = min_period
        self.max_period = max_period
        self.rotate_to_rtz = rotate_to_rtz
        self.unit_output = unit_output.upper()
        self.observed_tag = observed_tag

        # Allow manual override of synthetic tag, but keep internal and rely
        # on calling property for actual value
        self._synthetic_tag = synthetic_tag

        self.adj_src_type = adj_src_type
        self.st_obs_type = st_obs_type
        self.st_syn_type = st_syn_type
        self.win_amp_ratio = win_amp_ratio
        self.component_list = component_list

        # To be filled in by reading or with default parameters
        self.pyflex_config = None
        self.pyadjoint_config = None

        # If reading from a YAML file or from a dataset, do not set the external
        # Configs (pyflex and pyadjoint) because these will be read in verbatim
        if ds or yaml_fid:
            if ds:
                assert(path is not None), "'path' required to load from dataset"
                self._read_asdf(ds, path=path)
            elif yaml_fid:
                self._read_yaml(yaml_fid)
        # If initiating normally, need to set external Configs based on map
        # names and keyword arguments
        else:
            # Set Pyflex and Pyadjoint Config objects as attributes
            pyflex_parameters = pyflex_parameters or {}
            self.pyflex_config = PyflexConfig(min_period=min_period,
                                              max_period=max_period,
                                              **pyflex_parameters)
            pyadjoint_parameters = pyadjoint_parameters or {}
            # Double difference flag will be set by the adjoint source type
            self.pyadjoint_config = get_pyadjoint_config(
                adjsrc_type=adj_src_type, min_period=min_period,
                max_period=max_period, **pyadjoint_parameters
            )

        # Run internal sanity checks
        self._check()

    def __str__(self):
        """
        String representation of the class for print statements.
        It separates information into similar bins for readability.
        """
        # Model and step need to be formatted before printing
        str_out = ("CONFIG\n"
                   f"    {'iteration:':<25}{self.iter_tag}\n"
                   f"    {'step_count:':<25}{self.step_tag}\n"
                   f"    {'event_id:':<25}{self.event_id}\n"
                   )
        # Format the remainder of the keys identically
        key_dict = {"Process": ["min_period", "max_period",  "unit_output",
                                "rotate_to_rtz", "win_amp_ratio", "st_obs_type",
                                "st_syn_type"],
                    "Labels": ["component_list", "observed_tag",
                               "synthetic_tag"],
                    "External": ["adj_src_type", "pyflex_config",
                                 "pyadjoint_config"
                                 ]
                    }
        for key, items in key_dict.items():
            str_out += f"{key.upper()}\n"
            for item in items:
                str_out += f"    {item+':':<25}{getattr(self, item)}\n"
        return str_out

    def __repr__(self):
        """Simple call string representation"""
        return self.__str__()

    @property
    def pfcfg(self):
        """simple dictionary print of pyflex config object"""
        return vars(self.pyflex_config)

    @property
    def pacfg(self):
        """simple dictionary print of pyflex config object"""
        return vars(self.pyadjoint_config)

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
    def eval_tag(self):
        """string formatted version of iter and step, e.g. 'i01s00'"""
        return f"{self.iter_tag}{self.step_tag}"

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

    def _check(self):
        """
        A series of sanity checks to make sure that the configuration parameters
        are set properly to avoid any problems throughout the workflow. Should
        normally be run after any parameters are changed to make sure that they
        are acceptable.
        """
        if self.iteration is not None:
            assert(self.iteration >= 1), "Iterations must start at 1"

        if self.step_count is not None:
            assert(self.step_count >= 0), "Step count must start from 0"

        # Check period range is acceptable
        if self.min_period and self.max_period:
            assert(self.min_period < self.max_period), \
                "min_period must be less than max_period"

        # Check if unit output properly set, dictated by ObsPy units
        acceptable_units = ['DISP', 'VEL', 'ACC']
        assert(self.unit_output in acceptable_units), \
            f"unit_output should be in {acceptable_units}"

        # Set the component list. Rotate component list if necessary
        if self.rotate_to_rtz:
            if not self.component_list:
                self.component_list = ["R", "T", "Z"]
            else:
                for comp in ["N", "E"]:
                    assert(comp not in self.component_list), \
                            f"rotated component list cannot include '{comp}'"
        else:
            if not self.component_list:
                self.component_list = ["E", "N", "Z"]

        # Check that the amplitude ratio is a reasonable number
        if self.win_amp_ratio > 0:
            assert(self.win_amp_ratio < 1), \
                "window amplitude ratio should be < 1"

        assert(self.adj_src_type in ADJSRC_TYPES), \
            f"Pyadjoint `adj_src_type` must be in {ADJSRC_TYPES}"

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
        file format to use. Currently accepted file formats are yaml, asdf and
        ascii.

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

    def copy(self):
        """
        Simply convenience function to return a deep copy of the Config
        """
        return deepcopy(self)

    def write(self, write_to, fmt=None):
        """
        Wrapper for underlying low-level write functions

        :type fmt: str
        :param fmt: format to save parameters to. Available:

            * yaml: Write all parameters to a .yaml file which can be read later
            * ascii: Write parameters to a simple ascii file, not very smart and
              yaml is prefereable in most cases
            * asdf: Save the Config into an ASDFDataSet under the auxiliary
              data attribute
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
        Wrapper for underlying low-level read functions

        :type read_from: str or pyasdf.asdf_data_set.ASDFDataSet
        :param read_from: filename to read config from, or ds to read from
        :type path: str
        :param path: if fmt='asdf', path to the config in the aux data
        :type fmt: str
        :param fmt: file format to read parameters from, will be guessed but
            can also be explicitely set (available: 'yaml', 'ascii', 'asdf')
        """
        fmt = self._check_io_format(read_from, fmt)

        if fmt.lower() == "yaml":
            try:
                self._read_yaml(read_from)
            except ValueError as e:
                print(f"Unknown yaml format for file {read_from}, {e}")
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

        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset to save the config file to
        """
        # Deep copy to ensure that we aren't editing the Config parameters
        attrs = vars(deepcopy(self))

        add_attrs = {}
        del_attrs = []
        for key, item in attrs.items():
            if item is None:
                # HDF doesn't support NoneType so convert to string
                attrs[key] = "None"
            elif isinstance(item, dict) or ("config" in key):
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

        ds.add_auxiliary_data(data_type="Configs", data=np.array([True]),
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

        :type filename: str
        :param filename: filename to save yaml file
        :rtype: dict
        :return: key word arguments that do not belong to Pyatoa are passed back
            as a dictionary object, these are expected to be arguments that are
            to be used in Pyflex and Pyadjoint configs
        :raises ValueError: if unrecognized kwargs are found in the yaml file
        """
        with open(filename, "r") as f:
            attrs = yaml.load(f, Loader=yaml.Loader)

        unused_kwargs = {}
        for key, item in attrs.items():
            if hasattr(self, key.lower()):
                setattr(self, key.lower(), item)
            else:
                unused_kwargs[key.lower()] = item

        if unused_kwargs:
            raise ValueError(f"{list(unused_kwargs)} are not recognized "
                             "keyword arguments for a Config yaml file. Maybe "
                             "you meant to use the parameter 'seisflows_yaml'"
                             )

    def _read_asdf(self, ds, path):
        """
        Read and set config parameters from an ASDF Dataset, assumes that all
        necessary parameters are located in the auxiliary data subgroup of the
        dataset, which will be the case if the write_to_asdf() function was used
        Assumes some things about the structure of the auxiliary data.

        :type ds: pyasdf.asdf_data_set.ASDFDataSet
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
        pyflex_config, pyadjoint_config = {}, {}

        for key, item in cfgin.items():
            # Convert the item into expected native Python objects
            if isinstance(item, str):
                item = None if (item == "None" or item == "") else item
            else:
                try:
                    item = item.item()
                except ValueError:
                    item = item.tolist()

            # Put the item in the correct dictionary
            if "pyflex_config" in key:
                # Ensure that empties are set to NoneType
                pyflex_config["_".join(key.split('_')[2:])] = item
            elif "pyadjoint_config" in key:
                # e.g. pyadjoint_config_dlna_sigma_min -> dlna_sigma_min
                pyadjoint_config["_".join(key.split('_')[2:])] = item
            else:
                # Normal Config attribute
                setattr(self, key, item)

        # Set Pyflex and Pyadjoint Config objects as attributes
        self.pyflex_config = PyflexConfig(**pyflex_config)
        # Double difference is stored but not required
        pyadjoint_config.pop("double_difference")
        self.pyadjoint_config = get_pyadjoint_config(**pyadjoint_config)

