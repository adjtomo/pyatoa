:py:mod:`pyatoa`
================

.. py:module:: pyatoa


Subpackages
-----------
.. toctree::
   :titlesonly:
   :maxdepth: 3

   core/index.rst
   plugins/index.rst
   scripts/index.rst
   tests/index.rst
   utils/index.rst
   visuals/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.Config
   pyatoa.Manager
   pyatoa.Executive
   pyatoa.Gatherer
   pyatoa.Inspector




Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.logger
   pyatoa.propagate
   pyatoa.ch
   pyatoa.FORMAT
   pyatoa.formatter


.. py:data:: logger
   

   

.. py:data:: propagate
   :annotation: = 0

   

.. py:data:: ch
   

   

.. py:data:: FORMAT
   :annotation: = [%(asctime)s] - %(name)s - %(levelname)s: %(message)s

   

.. py:data:: formatter
   

   

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



.. py:class:: Manager(config=None, ds=None, event=None, st_obs=None, st_syn=None, inv=None, windows=None, staltas=None, adjsrcs=None, gcd=None, baz=None, gatherer=None)

   Pyatoas core workflow object.

   Manager is the central workflow control object. It calls on mid and
   low level classes to gather data, standardize and preprocess stream objects,
   generate misfit windows, and calculate adjoint sources. Has a variety of
   internal sanity checks to ensure that the workflow stays on the rails.

   .. py:property:: st

      Simplified call to return all streams available, observed and synthetic

   .. py:method:: __str__()

      Print statement shows available data detailing workflow


   .. py:method:: __repr__()

      Return repr(self).


   .. py:method:: check()

      (Re)check the stats of the workflow and data within the Manager.

      Rechecks conditions whenever called, incase something has gone awry
      mid-workflow. Stats should only be set by this function.


   .. py:method:: reset()

      Restart workflow by deleting all collected data in the Manager, but
      retain dataset, event, config, and gatherer so a new station can be
      processed with the same configuration as the previous workflow.


   .. py:method:: write(ds=None)

      Write the data collected inside Manager to an ASDFDataSet,

      :type ds: pyasdf.asdf_data_set.ASDFDataSet or None
      :param ds: write to a given ASDFDataSet. If None, will look for
          internal attribute `self.ds` to write to. Allows overwriting to
          new datasets


   .. py:method:: write_adjsrcs(path='./', write_blanks=True)

      Write internally stored adjoint source traces into SPECFEM defined
      two-column ascii files. Filenames are based on what is expected by
      Specfem, that is: 'NN.SSS.CCC.adj'

      ..note::
          By default writes adjoint sources for ALL components if one
          component has an adjoint source. If an adjoint sourced doesn't exist
          for a given component, it will be written with zeros. This is to
          satisfy SPECFEM3D requirements.

      :type path: str
      :param path: path to save the
      :type write_blanks: bool
      :param write_blanks: write zeroed out adjoint sources for components
          with no adjoint sources to meet the requirements of SPECFEM3D.
          defaults to True


   .. py:method:: load(code=None, path=None, ds=None, synthetic_tag=None, observed_tag=None, config=True, windows=False, adjsrcs=False)

      Populate the manager using a previously populated ASDFDataSet.
      Useful for re-instantiating an existing workflow that has already
      gathered data and saved it to an ASDFDataSet.

      .. note::
          mgmt.load() will return example data with no dataset

      .. warning::
          Loading any floating point values may result in rounding errors.
          Be careful to round off floating points to the correct place before
          using in future work.

      :type code: str
      :param code: SEED conv. code, e.g. NZ.BFZ.10.HHZ
      :type path: str
      :param path: if no Config object is given during init, the User
          can specify the config path here to load data from the dataset.
          This skips the need to initiate a separate Config object.
      :type ds: None or pyasdf.asdf_data_set.ASDFDataSet
      :param ds: dataset can be given to load from, will not set the ds
      :type synthetic_tag: str
      :param synthetic_tag: waveform tag of the synthetic data in the dataset
          e.g. 'synthetic_m00s00'. If None given, will use `config` attribute.
      :type observed_tag: str
      :param observed_tag: waveform tag of the observed data in the dataset
          e.g. 'observed'. If None given, will use `config` attribute.
      :type config: bool
      :param config: load config from the dataset, defaults to True but
          can be set False if Config should be instantiated by the User
      :type windows: bool
      :param windows: load misfit windows from the dataset, defaults to False
      :type adjsrcs: bool
      :param adjsrcs: load adjoint sources from the dataset, defaults to False


   .. py:method:: flow(**kwargs)

      A convenience function to run the full workflow with a single command.
      Does not include gathering. Takes kwargs related to all underlying
      functions.
      .. code:: python
          mgmt = Manager()
          mgmt.flow() == mgmt.standardize().preprocess().window().measure()
      :raises ManagerError: for any controlled exceptions


   .. py:method:: flow_multiband(periods, plot=False, **kwargs)

      Run the full workflow for a number of distinct period bands, returning
      a final set of adjoint sources generated as a summation of adjoint
      sources from each of these period bands.

      .. rubric::
          manager.flow_multiband(periods=[(1, 5), (10, 30), (40, 100)])

      :type periods: list of tuples
      :param periods: a list of tuples that define multiple period bands to
          generate windows and adjoint sources for. Overwrites the Config's
          internal `min_period` and `max_period` parameters. The final
          adjoint source will be a summation of all adjoint sources generated.
      :type plot: str
      :param plot: name of figure if given, will plot waveform and map for
          each period band and append period band to figure name `plot`
      :rtype: tuple of dict
      :return: (windows, adjoint_sources), returns all the collected
          measurements from each of the period bands
      :raises ManagerError: for any controlled exceptions


   .. py:method:: gather(code=None, choice=None, event_id=None, **kwargs)

      Gather station dataless and waveform data using the Gatherer class.
      In order collect observed waveforms, dataless, and finally synthetics.

      For valid kwargs see methods in :doc:`core.gatherer`

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type choice: list
      :param choice: allows user to gather individual bits of data, rather
          than gathering all. Allowed: 'inv', 'st_obs', 'st_syn'
      :raises ManagerError: if any part of the gathering fails.

      :keyword bool try_fm: Try to retrieve and append focal mechanism information to the
                            Event object.
      :keyword str prefix: Prefix for event id when searching for event information,
                           can be used to search ordered files e.g., CMTSOLUTION_001
      :keyword str suffix: Suffix for event id when searching for event information
      :keyword str station_level: The level of the station metadata if retrieved using the ObsPy
                                  Client. Defaults to 'response'
      :keyword str resp_dir_template: Directory structure template to search for response files.
                                      By default follows the SEED convention:
                                      'path/to/RESPONSE/{sta}.{net}/'
      :keyword str resp_fid_template: Response file naming template to search for station dataless.
                                      By default, follows the SEED convention
                                      'RESP.{net}.{sta}.{loc}.{cha}'
      :keyword str obs_dir_template: directory structure to search for observation data. Follows the
                                     SEED convention: 'path/to/obs_data/{year}/{net}/{sta}/{cha}'
      :keyword str obs_fid_template: File naming template to search for observation data. Follows the
                                     SEED convention: '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
      :keyword str syn_cfgpath: Config.cfgpaths key to search for synthetic data. Defaults to
                                'synthetics', but for the may need to be set to 'waveforms' in
                                certain use-cases, e.g. synthetics-synthetic inversions.
      :keyword str syn_unit: Optional argument to specify the letter used to identify the
                             units of the synthetic data: For Specfem3D: ["d", "v", "a", "?"]
                             'd' for displacement, 'v' for velocity,  'a' for acceleration.
                             Wildcards okay. Defaults to '?'
      :keyword str syn_dir_template: Directory structure template to search for synthetic waveforms.
                                     Defaults to empty string
      :keyword str syn_fid_template: The naming template of synthetic waveforms defaults to:
                                     "{net}.{sta}.*{cmp}.sem{syn_unit}"


   .. py:method:: standardize(force=False, standardize_to='syn')

      Standardize the observed and synthetic traces in place.
      Ensures Streams have the same starttime, endtime, sampling rate, npts.

      :type force: bool
      :param force: allow the User to force the function to run even if checks
          say that the two Streams are already standardized
      :type standardize_to: str
      :param standardize_to: allows User to set which Stream conforms to which
          by default the Observed traces should conform to the Synthetic ones
          because exports to Specfem should be controlled by the Synthetic
          sampling rate, npts, etc.


   .. py:method:: preprocess(which='both', overwrite=None, **kwargs)

      Preprocess observed and synthetic waveforms in place.
      Default preprocessing tasks: Remove response (observed), rotate, filter,
      convolve with source time function (synthetic).

      .. note::
          Default preprocessing can be overwritten using a
          user-defined function that takes Manager and choice as inputs
          and outputs an ObsPy Stream object.

      .. note::
          Documented kwargs only apply to default preprocessing.

      :type which: str
      :param which: "obs", "syn" or "both" to choose which stream to process
          defaults to both
      :type overwrite: function
      :param overwrite: If a function is provided, it will overwrite the
          standard preprocessing function. All arguments that are given
          to the standard preprocessing function will be passed as kwargs to
          the new function. This allows for customized preprocessing

      :keyword int water_level: water level for response removal
      :keyword float taper_percentage: amount to taper ends of waveform
      :keyword bool remove_response: remove instrument response using the Manager's inventory object.
                                     Defaults to True
      :keyword bool apply_filter: filter the waveforms using the Config's min_period and
                                  max_period parameters. Defaults to True
      :keyword bool convolve_with_stf: Convolve synthetic data with a Gaussian source time function if
                                       a half duration is provided.
      :keyword bool rotate_to_rtz: Use the `rotate_baz` variable to rotate streams
                                   from ZNE components to RTZ


   .. py:method:: window(fix_windows=False, iteration=None, step_count=None, force=False, save=True)

      Evaluate misfit windows using Pyflex. Save windows to ASDFDataSet.
      Allows previously defined windows to be retrieved from ASDFDataSet.

      .. note::
          * Windows are stored as dictionaries of pyflex.Window objects.
          * All windows are saved into the ASDFDataSet, even if retrieved.
          * STA/LTA information is collected and stored internally.

      :type fix_windows: bool
      :param fix_windows: do not pick new windows, but load windows from the
          given dataset from 'iteration' and 'step_count'
      :type iteration: int or str
      :param iteration: if 'fix_windows' is True, look for windows in this
          iteration. If None, will check the latest iteration/step_count
          in the given dataset
      :type step_count: int or str
      :param step_count: if 'fix_windows' is True, look for windows in this
          step_count. If None, will check the latest iteration/step_count
          in the given dataset
      :type force: bool
      :param force: ignore flag checks and run function, useful if e.g.
          external preprocessing is used that doesn't meet flag criteria
      :type save: bool
      :param save: save the gathered windows to an ASDF Dataset


   .. py:method:: retrieve_windows(iteration, step_count, return_previous)

      Mid-level window selection function that retrieves windows from a
      PyASDF Dataset, recalculates window criteria, and attaches window
      information to Manager. No access to rejected window information.

      :type iteration: int or str
      :param iteration: retrieve windows from the given iteration
      :type step_count: int or str
      :param step_count: retrieve windows from the given step count
          in the given dataset
      :type return_previous: bool
      :param return_previous: if True: return windows from the previous
          step count in relation to the given iteration/step_count.
          if False: return windows from the given iteration/step_count


   .. py:method:: select_windows_plus()

      Mid-level custom window selection function that calls Pyflex select
      windows, but includes additional window suppression functionality.
      Includes custom Pyflex addition of outputting rejected windows, which
      will be used internally for plotting.

      .. note::
          Pyflex will throw a ValueError if the arrival of the P-wave
          is too close to the initial portion of the waveform, considered the
          'noise' section. This happens for short source-receiver distances
          (< 100km).

          This error becomes a PyflexError if no event/station attributes
          are provided to the WindowSelector

          We could potentially deal with this by zero-padding the
          waveforms, and running select_windows() again, but for now we just
          raise a ManagerError and allow processing to continue


   .. py:method:: measure(force=False, save=True)

      Measure misfit and calculate adjoint sources using PyAdjoint.

      Method for caluculating misfit set in Config, Pyadjoint expects
      standardized traces with the same spectral content, so this function
      will not run unless these flags are passed.

      Returns a dictionary of adjoint sources based on component.
      Saves resultant dictionary to a pyasdf dataset if given.

      .. note::
          Pyadjoint returns an unscaled misfit value for an entire set of
          windows. To return a "total misfit" value as defined by
          Tape (2010) Eq. 6, the total summed misfit will need to be scaled by
          the number of misfit windows chosen in Manager.window().

      :type force: bool
      :param force: ignore flag checks and run function, useful if e.g.
          external preprocessing is used that doesn't meet flag criteria
      :type save: bool
      :param save: save adjoint sources to ASDFDataSet


   .. py:method:: save_windows(ds=None, force=False)

      Convenience function to save collected misfit windows into an
      ASDFDataSet with some preliminary checks

      Auxiliary data tag is hardcoded as 'MisfitWindows'

      :type ds: pyasdf.ASDFDataSet
      :param ds: allow replacement of the internal `ds` dataset. If None,
          will try to write to internal `ds`
      :type force: bool
      :param force: force saving windows even if Config says don't do it.
          This is used by write() to bypass the default 'dont save' behavior


   .. py:method:: save_adjsrcs(ds=None, force=False)

      Convenience function to save collected adjoint sources into an
      ASDFDataSet with some preliminary checks

      Auxiliary data tag is hardcoded as 'AdjointSources'

      :type ds: pyasdf.ASDFDataSet
      :param ds: allow replacement of the internal `ds` dataset. If None,
          will try to write to internal `ds`
      :type force: bool
      :param force: force saving windows even if Config says don't do it.
          This is used by write() to bypass the default 'dont save' behavior


   .. py:method:: _format_windows()

      .. note::
          In `pyadjoint.calculate_adjoint_source`, the window needs to be a
          list of lists, with each list containing the
          [left_window, right_window]; each window argument should be given in
          units of time (seconds). This is not in the PyAdjoint docs.

      :rtype: dict of list of lists
      :return: dictionary with key related to individual components,
          and corresponding to a list of lists containing window start and end


   .. py:method:: plot(choice='both', save=None, show=True, corners=None, figsize=None, dpi=100, **kwargs)

      Plot observed and synthetics waveforms, misfit windows, STA/LTA and
      adjoint sources for all available components. Append information
      about misfit, windows and window selection. Also as subplot create a
      source receiver map which contains annotated information detailing
      src-rcv relationship like distance and BAz. Options to plot either or.

      For valid key word arguments see `visuals.manager_plotter` and
      `visuals.map_maker`

      :type show: bool
      :param show: show the plot once generated, defaults to False
      :type save: str
      :param save: absolute filepath and filename if figure should be saved
      :param corners: {lat_min, lat_max, lon_min, lon_max}
          corners to cut the map to, otherwise a global map is provided
      :type choice: str
      :param choice: choice for what to plot:
          * 'wav': plot waveform figure only
          * 'map': plot a source-receiver map only
          * 'both' (default): plot waveform and source-receiver map together
      :type figsize: tuple
      :param figsize: optional size of the figure, set by plot()
      :type dpi: int
      :param dpi: optional dots per inch (resolution) of figure



.. py:exception:: ManagerError

   Bases: :py:obj:`Exception`

   A class-wide custom exception raised when functions fail gracefully


.. py:class:: Executive(event_ids, station_codes, config, max_stations=4, max_events=1, cat='+', log_level='DEBUG', cwd=None, datasets=None, figures=None, logs=None, adjsrcs=None, ds_fid_template=None)

   The Executive is hierarchically above Pyatoa's core class, the Manager.
   It sets up a simple framework to organize and parallelize misfit
   quantification.

   .. py:property:: codes

      Define a set of event-station codes that are used to traverse through
      all possible source receiver combinations.

      .. note::
          Workaround for having it be pretty difficult to pass multiple
          arguments into an executor. Just pass a list of strings that is
          split by the parallel processes

   .. py:method:: check()

      Parameter checking


   .. py:method:: process()

      Process all events concurrently


   .. py:method:: process_event(event_id)

      Process all given stations concurrently for a single event

      :type event_id: str
      :param event_id: one value from the Executor.events list specifying
          a given event to process


   .. py:method:: process_station(event_id_and_station_code)

      Parallel process multiple Managers simultaneously, which is the biggest
      time sync. IO is done in serial to get around BlockingIO

      .. note::
          Very broad exceptions to keep process running smoothly, you will
          need to check log messages individually to figure out if and where
          things did not work

      .. note::
          Employs a workaround to inability to parallel write to HDF5 files
          BlockingIOError by doing the processing first, and then waiting
          for each process to finish writing before accessing.

      :type event_id_and_station_code: str
      :param event_id_and_station_code: a string concatenation of a given
          event id and station code, which will be used to process a single
          source receiver pair


   .. py:method:: _check_rank(event_id_and_station_code)

      Poor man's method for determining the processor rank for a given event.
      Used so that processes that happen only once (e.g., writing config) are
      done consistently by one process

      :type event_id_and_station_code: str
      :param event_id_and_station_code: a string concatenation of a given
          event id and station code, which will be used to process a single
          source receiver pair
      :rtype: int
      :return: rank index in Executive.codes based on event and station


   .. py:method:: _generate_logger(log_path)

      Create a log file for each source. No stream handler, only file output
      Also create a memory handler to dump all log messages at once, rather
      than as they happen, allowing multiple stations to write to the same
      file sequentially

      :type log_path: str
      :param log_path: path and filename to save log file



.. py:class:: Gatherer(config, ds=None, origintime=None)

   A mid-level data gathering class used to get data internally and externally.
   All saving to ASDFDataSet taken care of by the Gatherer class.

   .. py:method:: gather_event(event_id=None, **kwargs)

      Gather an ObsPy Event object by searching disk
      Event info need only be retrieved once per Pyatoa workflow.

      :type event_id: str
      :param event_id: a unique event idenfitier to search and tag event info
      :rtype: obspy.core.event.Event
      :return: event retrieved either via internal or external methods
      :raises GathererNoDataException: if no event information is found.


   .. py:method:: gather_station(code, **kwargs)

      Gather StationXML information from disk

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.inventory.Inventory
      :return: inventory containing relevant network and stations


   .. py:method:: gather_observed(code, **kwargs)

      Gather observed waveforms from disk as ObsPy streams

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.stream.Stream
      :return: stream object containing relevant waveforms


   .. py:method:: gather_synthetic(code, **kwargs)

      Gather synthetic waveforms as ObsPy streams.

      Only possible to check ASDFDataSet and local filesystem, cannot gather
      synthetics from webservice.

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.stream.Stream
      :return: stream object containing relevant waveforms
      :raises GathererNoDataException: if no synthetic data is found


   .. py:method:: fetch_event_from_dataset()

      Return Event information from ASDFDataSet.

      .. note::
          Assumes that the ASDF Dataset will only contain one event, which is
          dictated by the structure of Pyatoa.

      :rtype event: obspy.core.event.Event
      :return event: event object
      :raises AttributeError: if no event attribute found in ASDFDataSet
      :raises IndexError: if event attribute found but no events


   .. py:method:: fetch_inv_from_dataset(code)

      Return StationXML from ASDFDataSet based on station code.

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.inventory.network.Network
      :return: network containing relevant station information
      :raises KeyError: if no matching StationXML found


   .. py:method:: fetch_waveform_from_dataset(code, tag)

      Return waveforms as Stream objects from ASDFDataSet.

      .. note:
          * Allows for wildcard selection of component (? or *)
          * Selects by component because synthetic channel naming may differ
          from observation channels.
          * Component is assumed to be the last index in the channel,
          following SEED convention.

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type tag: str
      :param tag: internal asdf tag labelling waveforms
      :rtype: obspy.core.stream.Stream
      :return: waveform contained in a stream, or None if no matching value


   .. py:method:: fetch_event_by_dir(event_id, prefix='', suffix='', format_=None, **kwargs)

      Fetch event information via directory structure on disk. Developed to
      parse CMTSOLUTION and QUAKEML files, but theoretically accepts any
      format that the ObsPy read_events() function will accept.

      Will search through all paths given until a matching source file found.

      .. note::
          This function will search for the following path
          /path/to/event_dir/{prefix}{event_id}{suffix}

          so, if e.g., searching for a CMTSOLUTION file in the current dir:
          ./CMTSOLUTION_{event_id}

          Wildcards are okay but the function will return the first match

      :type event_id: str
      :param event_id: Unique event identifier to search source file by.
          e.g., a New Zealand earthquake ID '2018p130600'. A prefix or suffix
          will be tacked onto this
      :rtype event: obspy.core.event.Event or None
      :return event: event object if found, else None.
      :type prefix: str
      :param prefix Prefix to prepend to event id for file name searching.
          Wildcards are okay.
      :type suffix: str
      :param suffix: Suffix to append to event id for file name searching.
          Wildcards are okay.
      :type format_: str or NoneType
      :param format_: Expected format of the file to read, e.g., 'QUAKEML',
          passed to ObsPy read_events. NoneType means read_events() will guess


   .. py:method:: fetch_inv_by_dir(code, resp_dir_template='{sta}.{net}', resp_fid_template='RESP.{net}.{sta}.{loc}.{cha}', **kwargs)

      Fetch station dataless via directory structure on disk.
      Will search through all paths given until StationXML found.

      .. note::
          Default path naming follows SEED convention, that is:
          path/to/dataless/{NET}.{STA}/RESP.{NET}.{STA}.{LOC}.{CHA}
          e.g. path/to/dataless/NZ.BFZ/RESP.NZ.BFZ.10.HHZ

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type resp_dir_template: str
      :param resp_dir_template: Directory structure template to search for
          response files. By default follows the SEED convention:
          'path/to/RESPONSE/{sta}.{net}/'
      :type resp_fid_template: str
      :param resp_fid_template: Response file naming template to search for
          station dataless. By default, follows the SEED convention:
          'RESP.{net}.{sta}.{loc}.{cha}'
      :rtype inv: obspy.core.inventory.Inventory or None
      :return inv: inventory containing relevant network and stations


   .. py:method:: fetch_observed_by_dir(code, obs_dir_template='{year}/{net}/{sta}/{cha}', obs_fid_template='{net}.{sta}.{loc}.{cha}.{year}.{jday:0>3}', **kwargs)

      Fetch observation waveforms via directory structure on disk.

      .. note::
          Default waveform directory structure assumed to follow SEED
          convention. That is:
          path/to/data/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/{FID}
          e.g. path/to/data/2017/NZ/OPRZ/HHZ.D/NZ.OPRZ.10.HHZ.D

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type obs_dir_template: str
      :param obs_dir_template: directory structure to search for observation
          data. Follows the SEED convention:
          'path/to/obs_data/{year}/{net}/{sta}/{cha}'
      :type obs_fid_template: str
      :param obs_fid_template: File naming template to search for observation
          data. Follows the SEED convention:
          '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
      :rtype stream: obspy.core.stream.Stream or None
      :return stream: stream object containing relevant waveforms, else None


   .. py:method:: fetch_synthetic_by_dir(code, syn_cfgpath='synthetics', syn_unit='?', syn_dir_template='', syn_fid_template='{net}.{sta}.*{cmp}.sem{dva}*', **kwargs)

      Fetch synthetic waveforms from Specfem3D via directory structure on
      disk, if necessary convert native ASCII format to Stream object.

      .. note::
          By default, synthetics will be searched for with the following path
          config.paths[syn_cfgpath]/syn_dir_template/syn_fid_template.format()

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type syn_cfgpath: str
      :param syn_cfgpath: Config.paths key to search for synthetic data.
          Defaults to 'synthetics', but for the may need to be set to
          'waveforms' in certain use-cases.
      :type syn_unit: str
      :param syn_unit: Optional argument to specify the letter used to
          identify the units of the synthetic data: For Specfem3D:
          ["d", "v", "a", "?"] 'd' for displacement, 'v' for velocity,
          'a' for acceleration. Wildcards okay. Defaults to '?'
      :type syn_dir_template: str
      :param syn_dir_template: Directory structure template to search for
          synthetic waveforms. Defaults to empty string
      :type syn_fid_template: str
      :param syn_fid_template: The naming template of synthetic waveforms
          defaults to "{net}.{sta}.*{cmp}.sem{syn_unit}"
      :rtype stream: obspy.core.stream.Stream or None
      :return stream: stream object containing relevant waveforms


   .. py:method:: save_waveforms_to_dataset(st, tag)

      Save waveformsm to the ASDFDataSet with a simple check for existence
      of dataset and save parameter. Passes if waveforms already exist while
      ignoring the PyASDF warning that gets thrown if waveforms exist.

      :type st: obspy.core.stream.Stream
      :param st: Stream object to be saved into the dataset
      :type tag: str
      :param tag: unique identifier to save the waveforms under



.. py:class:: Inspector(tag='default', verbose=True)

   Bases: :py:obj:`pyatoa.visuals.insp_plot.InspectorPlotter`

   This plugin object will collect information from a Pyatoa run folder and
   allow the User to easily understand statistical information or generate
   statistical plots to help understand a seismic inversion.

   Inherits plotting capabilities from InspectorPlotter class to reduce clutter

   .. py:property:: keys

      Shorthand to access the keys of the Windows dataframe

   .. py:property:: events

      Return an array of all event ids

   .. py:property:: stations

      Return an array of all stations

   .. py:property:: networks

      Return an array of all stations

   .. py:property:: netsta

      Return a Dataframe containing unique network-station idents

   .. py:property:: srcrcv

      Return a dataframe with source-receiver information, dists and baz

   .. py:property:: pairs

      Determine the number of unique source-receiver pairs

   .. py:property:: iterations

      Return an array of all iteration

   .. py:property:: steps

      Returns a pandas. Series of iteration with values listing steps

   .. py:property:: models

      Return a dict of model numbers related to a unique iteration/step

   .. py:property:: initial_model

      Return tuple of the iteration and step count corresponding M00

   .. py:property:: final_model

      Return tuple of iteration and step count for final accepted model

   .. py:property:: good_models

      Return models that are only status 0 or 1 (initial or success)

   .. py:property:: restarts

      Try to guess the indices of restarts for convergence plot based on
      misfit increase in adjacent good models as well as discontinous misfit
      values for the final line search model and subsequent initial model.
      Not guaranteed to catch everything so may require manual review using
      the convergence() function

   .. py:property:: evaluations

      Returns the number of iterations, or the sum of all step counts

   .. py:property:: mags

      Return a dictionary of event magnitudes

   .. py:property:: times

      Return a dictionary of event origin times

   .. py:property:: depths

      Return a dictionary of event depths in units of meters

   .. py:method:: _get_str()

      Get the string representation once and save as internal attribute


   .. py:method:: __str__()

      Return a list of all variables and functions available for quick ref


   .. py:method:: __repr__()

      Return repr(self).


   .. py:method:: _try_print(a)

      Try-except catch for property print statements


   .. py:method:: _get_srcrcv_from_dataset(ds)

      Get source and receiver information from dataset, this includes
      latitude and longitude values for both, and event information including
      magnitude, origin time, id, etc.

      Returns Dataframes for sources and receivers iff they are not already
      contained in the class dataframes, to avoid duplicates.

      Returns empty DataFrames if no unique info was found.

      :type ds: pyasdf.ASDFDataSet
      :param ds: dataset to query for distances
      :rtype source: pandas.core.frame.DataFrame
      :return source: single row Dataframe containing event info from dataset
      :rtype receivers: multiindexed dataframe containing unique station info


   .. py:method:: _get_windows_from_dataset(ds)

      Get window and misfit information from dataset auxiliary data
      Model and Step information should match between the two
      auxiliary data objects MisfitWindows and AdjointSources

      TODO: break this into _get_windows_from_dataset and
            _get_adjsrcs_from_dataset?

      :type ds: pyasdf.ASDFDataSet
      :param ds: dataset to query for misfit:
      :rtype: pandas.DataFrame
      :return: a dataframe object containing information per misfit window


   .. py:method:: _parse_nonetype_eval(iteration, step_count)

      Whenever a user does not choose an iteration or step count, e.g., in
      plotting functions, this function defines default values based on the
      initial model (if neither given), or the last step count for a given
      iteration (if only iteration is given). Only step count is not allowed

      :type iteration: str
      :param iteration: chosen iteration, formatted as e.g., 'i01'
      :type step_count: str
      :param step_count: chosen step count, formatted as e.g., 's00'
      :rtype: tuple of str
      :return: (iteration, step_count) default values for the iteration
          and step_count


   .. py:method:: discover(path='./', ignore_symlinks=True)

      Allow the Inspector to scour through a path and find relevant files,
      appending them to the internal structure as necessary.

      :type path: str
      :param path: path to the pyasdf.asdf_data_set.ASDFDataSets that were
          outputted by the Seisflows workflow
      :type ignore_symlinks: bool
      :param ignore_symlinks: skip over symlinked HDF5 files when discovering


   .. py:method:: append(dsfid, srcrcv=True, windows=True)

      Simple function to parse information from a
      pyasdf.asdf_data_setASDFDataSet file and append it to the currect
      collection of information.

      :type dsfid: str
      :param dsfid: fid of the dataset
      :type srcrcv: bool
      :param srcrcv: gather source-receiver information
      :type windows: bool
      :param windows: gather window information


   .. py:method:: extend(windows)

      Extend the current Inspector data frames with the windows from another
      Inspector. This is useful for when an inversion has been run in legs, so
      two individual inspectors constitute a single inversion.

      .. note::
          The current inspector is considered leg A, and the argument
          'windows' is considered leg B. Leg B will have its iteration numbers
          changed to reflect this

      .. warning::
          This will only work if all the events and stations are the same.
          That is, only two identical inversion scenarios can be used.

      :type windows: pandas.core.data_frame.DataFrame or list of DataFrames
      :param windows: Windows from a separate inspector object that will be
          used to extend the current Inspector. Can also be provided as a list
          of DataFrames to extend multiple times.


   .. py:method:: save(path='./', fmt='csv', tag=None)

      Save the downloaded attributes into JSON files for easier re-loading.

      .. note::
          fmt == 'hdf' requires 'pytables' to be installed in the environment

      :type tag: str
      :param tag: tag to use to save files, defaults to the class tag
          but allows for the option of overwriting that
      :type path: str
      :param path: optional path to save to, defaults to cwd
      :type fmt: str
      :param fmt: format of the files to write, default csv


   .. py:method:: write(**kwargs)

      Same as Inspector.save(), but I kept writing .write()


   .. py:method:: read(path='./', fmt=None, tag=None)

      Load previously saved attributes to avoid re-processing data.

      :type tag: str
      :param tag: tag to use to look for files, defaults to the class tag
          but allows for the option of overwriting that
      :type path: str
      :param path: optional path to file, defaults to cwd
      :type fmt: str
      :param fmt: format of the files to read, default csv


   .. py:method:: copy()

      Return a deep copy of the Inspector


   .. py:method:: reset()

      Simple function to wipe out all the internal attributes


   .. py:method:: isolate(iteration=None, step_count=None, event=None, network=None, station=None, channel=None, component=None, keys=None, exclude=None, unique_key=None)

      Returns a new dataframe that is grouped by a given index if variable is
      None, defaults to returning all available values

      :type event: str
      :param event: event id e.g. '2018p130600' (optional)
      :type iteration: str
      :param iteration: iteration e.g. 'i00' (optional)
      :type step_count: str
      :param step_count: step count e.g. 's00' (optional)
      :type station: str
      :param station: station name e.g. 'BKZ' (optional)
      :type network: str
      :param network: network name e.g. 'NZ' (optional)
      :type channel: str
      :param channel: channel name e.g. 'HHE' (optional)
      :type component: str
      :param component: component name e.g. 'Z' (optional)
      :type unique_key: str
      :param unique_key: isolates model, event and station information,
          alongside a single info key, such as dlnA.
          Useful for looking at one variable without have to write out long
          lists to 'exclude' or 'keys'
      :type keys: list
      :param keys: list of keys to retain in returned dataset, 'exclude'
          will override this variable, best to use them separately
      :type exclude: list
      :param exclude: list of keys to remove from returned dataset
      :rtype: pandas.DataFrame
      :return: DataFrame with selected rows based on selected column values


   .. py:method:: nwin(level='step')

      Find the cumulative length of misfit windows for a given iter/step,
      or the number of misfit windows for a given iter/step.

      .. note::
          Neat trick to select just by station:
          insp.windows(level='station').query("station == 'BFZ'")

      :type level: str
      :param level: Level to get number of windows by. Default is 'step'

          * step: to get the total window length and number of windows for the
            given step count.
          * station: to get this on a per-station basis,
            useful for identifying sta quality.
      :rtype: pandas.DataFrame
      :return: a DataFrame with indices corresponding to iter, step,
          columns listing the number of windows (nwin) and the cumulative
          length of windows in seconds (length_s)


   .. py:method:: misfit(level='step', reset=False)

      Sum the total misfit for a given iteration based on the individual
      misfits for each misfit window, and the number of sources used.
      Calculated misfits are stored internally to avoid needing to recalculate
      each time this function is called

      .. note::
          To get per-station misfit on a per-step basis
              df = insp.misfits(level="station").query("station == 'TOZ'")
              df.groupby(['iteration', 'step']).sum()

      :type level: str
      :param level:  Default is 'step'
          'station': unscaled misfit on a per-station basis
          'step': to get total misfit for a given step count.
          'event': to get this on a per-event misfit.
      :type reset: bool
      :param reset: reset internally stored attribute and re-calculate misfit
      :rtype: dict
      :return: total misfit for each iteration in the class


   .. py:method:: stats(level='event', choice='mean', key=None, iteration=None, step_count=None)

      Calculate the per-level statistical values for DataFrame

      :type level: str
      :param level: get statistical values per 'event' or 'station'
      :type choice: str
      :param choice: Pandas function, 'mean', 'std', 'var', etc.
      :type key: windows column header, e.g. 'cc_shift_in_seconds'
      :type iteration: str
      :param iteration: filter for a given iteration
      :type step_count: str
      :param step_count: filter for a given step count
      :rtype: pandas.DataFrame
      :return: DataFrame containing the `choice` of stats for given options


   .. py:method:: minmax(iteration=None, step_count=None, keys=None, quantities=None, pprint=True)

      Calculate and print the min/max values for a whole slew of parameters
      for a given iteration and step count. Useful for understanding the
      worst/ best case scenarios and their relation to the average.

      :type iteration: str
      :param iteration: filter for a given iteration
      :type step_count: str
      :param step_count: filter for a given step count
      :type keys: list of str
      :param keys: keys to calculate minmax values for, must be a subset of
          Inspector.windows.keys()
      :type quantities: list of str
      :param quantities: quantities to get values for, e.g. min, max, median,
          must be an attribute of pandas.core.series.Series
      :type pprint: bool
      :param pprint: pretty print the resulting values
      :rtype: dict
      :return: dictionary containing the minmax stats


   .. py:method:: compare(iteration_a=None, step_count_a=None, iteration_b=None, step_count_b=None)

      Compare the misfit and number of windows on an event by event basis
      between two evaluations. Provides absolute values as well as
      differences. Final dataframe is sorted by the difference in misfit,
      showing the most and least improved events.

      :type iteration_a: str
      :param iteration_a: initial iteration to use in comparison
      :type step_count_a: str
      :param step_count_a: initial step count to use in comparison
      :type iteration_b: str
      :param iteration_b: final iteration to use in comparison
      :type step_count_b: str
      :param step_count_b: final step count to use in comparison
      :rtype: pandas.core.data_frame.DataFrame
      :return: a sorted data frame containing the difference of misfit and
          number of windows between final and initial


   .. py:method:: compare_windows(iteration_a=None, step_count_a=None, iteration_b=None, step_count_b=None)

      Compare individual, matching misfit windows between two evaluations.

      .. note::
          This will only work/make sense if the windows were fixed between
          the two evaluations, such that they share the exact same window
          selections.

      :type iteration_a: str
      :param iteration_a: initial iteration to use in comparison
      :type step_count_a: str
      :param step_count_a: initial step count to use in comparison
      :type iteration_b: str
      :param iteration_b: final iteration to use in comparison
      :type step_count_b: str
      :param step_count_b: final step count to use in comparison
      :rtype: pandas.core.data_frame.DataFrame
      :return: a data frame containing differences of windowing paramenters
          between final and initial models


   .. py:method:: filter_sources(lat_min=None, lat_max=None, lon_min=None, lon_max=None, depth_min=None, depth_max=None, mag_min=None, mag_max=None, min_start=None, max_start=None)

      Go through misfits and windows and remove events that fall outside
      a certain bounding box. Return sources that fall within the box.
      Bounds are inclusive of given values.

      :type lat_min: float
      :param lat_min: minimum latitude in degrees
      :type lat_max: float
      :param lat_max: maximum latitude in degrees
      :type lon_min: float
      :param lon_min: minimum longitude in degrees
      :type lon_max: float
      :param lon_max: maximum longitude in degrees
      :type depth_min: float
      :param depth_min: minimum depth of event in km, depth is positive
      :type depth_max: float
      :param depth_max: maximum depth of event in km, depth is positive
      :type mag_min: float
      :param mag_min: minimum magnitude
      :type mag_max: float
      :param mag_max: maximum magnitude
      :type min_start: obspy.UTCDateTime()
      :param min_start: minimum origintime of event
      :type max_start: obspy.UTCDateTime()
      :param max_start: maximum origintime of event


   .. py:method:: get_models()

      Return a sorted list of misfits which correspond to accepted models,
      label discards of the line search, and differentiate the final accepted
      line search evaluation from the previous iteration and the initial
      evaluation of the current iteration.

      .. note::
          State and status is given as:
          0 == INITIAL function evaluation for the model;
          1 == SUCCESS -ful function evaluation for the model;
          -1 == DISCARD trial step from line search.

      :rtype: pandas.core.data_frame.DataFrame
      :return: a dataframe containing model numbers, their corresponding
          iteration, step count and misfit value, and the status of the
          function evaluation.


   .. py:method:: get_srcrcv()

      Retrieve information regarding source-receiver pairs including distance,
      backazimuth and theoretical traveltimes for a 1D Earth model.

      :rtype: pandas.core.frame.DataFrame
      :return: separate dataframe with distance and backazimuth columns, that
          may be used as a lookup table


   .. py:method:: get_unique_models(float_precision=3)

      Find all accepted models (status 0 or 1) that have a unique misfit
      value. Because some forward evaluations are repeats of the previous
      line search evaluation, they will effectively be the same evaluation so
      they can be removed

      :type float_precision: int
      :param float_precision: identical misfit values will differ after some
          decimal place. this value determines which decimal place to
          truncate the values for comparison



