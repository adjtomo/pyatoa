:py:mod:`pyatoa.core.config`
============================

.. py:module:: pyatoa.core.config

.. autoapi-nested-parse::

   Configuration of User-set parameters within the package.
   Contains external functions to set Config objects of Pyflex and Pyadjoint.

   To Do 02.04.21: Allow config to set pyflex and pyadjoint config as a function,
       currently it's obscured behind some private functions



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.core.config.Config



Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.core.config.set_pyflex_config
   pyatoa.core.config.set_pyadjoint_config



.. py:class:: Config(yaml_fid=None, ds=None, path=None, iteration=None, step_count=None, event_id=None, min_period=10, max_period=100, rotate_to_rtz=False, unit_output='DISP', component_list=None, pyflex_preset='default', adj_src_type='cc_traveltime', observed_tag='observed', synthetic_tag=None, st_obs_type='obs', st_syn_type='syn', win_amp_ratio=0.0, paths=None, save_to_ds=True, **kwargs)

   The Config class is the main interaction object between the User and
   workflow. It is used by :class:`Manager <pyatoa.core.manager.Manager>` for
   workflow management, and also for information sharing between Pyatoa objects
   and functions. The Config can be read to and written from external files and
   ASDFDataSets.

   .. py:property:: pfcfg

      simple dictionary print of pyflex config object

   .. py:property:: pacfg

      simple dictionary print of pyflex config object

   .. py:property:: iter_tag

      string formatted version of iteration, e.g. 'i00'

   .. py:property:: step_tag

      string formatted version of step, e.g. 's00'

   .. py:property:: eval_tag

      string formatted version of iter and step, e.g. 'i01s00'

   .. py:property:: synthetic_tag

      tag to be used for synthetic data, uses iteration and step count

   .. py:property:: aux_path

      property to quickly get a bog-standard aux path e.g. i00/s00

   .. py:method:: __str__()

      String representation of the class for print statements.
      It separates information into similar bins for readability.


   .. py:method:: __repr__()

      Simple call string representation


   .. py:method:: _check()

      A series of sanity checks to make sure that the configuration parameters
      are set properly to avoid any problems throughout the workflow. Should
      normally be run after any parameters are changed to make sure that they
      are acceptable.


   .. py:method:: _set_external_configs(check_unused=False, **kwargs)

      Set the Pyflex and Pyadjoint Config parameters using kwargs provided
      to the init function. Allows the user to request unused kwargs be
      returned with an Error statement.

      :type check_unnused: bool
      :param check_unnused: check if kwargs passed to Pyadjoint and Pyflex
          do not match config parameters for either, which may mean
          misspelled parameter names, or just kwargs that were not meant for
          either
      :raises ValueError: if check_unnused is True and unnused kwargs found


   .. py:method:: _get_aux_path(default='default', separator='/')

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


   .. py:method:: _check_io_format(fid, fmt=None)
      :staticmethod:

      A simple check before reading or writing the config to determine what
      file format to use. Currently accepted file formats are yaml, asdf and
      ascii.

      :type fmt: str
      :param fmt: format specified by the User
      :rtype: str
      :return: format string to be understood by the calling function


   .. py:method:: copy()

      Simply convenience function to return a deep copy of the Config


   .. py:method:: write(write_to, fmt=None)

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


   .. py:method:: read(read_from, path=None, fmt=None)

      Wrapper for underlying low-level read functions

      :type read_from: str or pyasdf.asdf_data_set.ASDFDataSet
      :param read_from: filename to read config from, or ds to read from
      :type path: str
      :param path: if fmt='asdf', path to the config in the aux data
      :type fmt: str
      :param fmt: file format to read parameters from, will be guessed but
          can also be explicitely set (available: 'yaml', 'ascii', 'asdf')


   .. py:method:: _write_yaml(filename)

      Write config parameters to a yaml file, retain order

      :type filename: str
      :param filename: filename to save yaml file


   .. py:method:: _write_asdf(ds)

      Save the Config values as a parameter dictionary in the ASDF Data set
      Converts types to play nice with ASDF Auxiliary Data.
      Flattens dictionaries and external Config objects for easy storage.

      :type ds: pyasdf.asdf_data_set.ASDFDataSet
      :param ds: dataset to save the config file to


   .. py:method:: _write_ascii(filename)

      Write the config parameters to an ascii file

      :type filename: str
      :param filename: filename to write the ascii file to


   .. py:method:: _read_yaml(filename)

      Read config parameters from a yaml file, parse to attributes.

      :type filename: str
      :param filename: filename to save yaml file
      :rtype: dict
      :return: key word arguments that do not belong to Pyatoa are passed back
          as a dictionary object, these are expected to be arguments that are
          to be used in Pyflex and Pyadjoint configs
      :raises ValueError: if unrecognized kwargs are found in the yaml file


   .. py:method:: _read_asdf(ds, path)

      Read and set config parameters from an ASDF Dataset, assumes that all
      necessary parameters are located in the auxiliary data subgroup of the
      dataset, which will be the case if the write_to_asdf() function was used
      Assumes some things about the structure of the auxiliary data.

      :type ds: pyasdf.asdf_data_set.ASDFDataSet
      :param ds: dataset with config parameter to read
      :type path: str
      :param path: model number e.g. 'm00' or 'default', or 'm00/s00'



.. py:function:: set_pyflex_config(min_period, max_period, choice=None, **kwargs)

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


.. py:function:: set_pyadjoint_config(adjsrc_type, min_period, max_period, **kwargs)

   Set the Pyadjoint config based on Pyatoa Config parameters.
   Kwargs can be fed to the Pyadjoint Config object. Returns unnused kwargs.

   Config parameters can be found at:
   http://adjtomo.github.io/pyadjoint/autoapi/pyadjoint/config/index.html

   :type min_period: float
   :param min_period: min period of the data
   :type max_period: float
   :param max_period: max period of the data
   :rtype cfgout: pyadjoint.Config
   :return cfgout: properly set pyadjoint configuration object


