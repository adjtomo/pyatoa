:py:mod:`pyatoa.core.manager`
=============================

.. py:module:: pyatoa.core.manager

.. autoapi-nested-parse::

   A class to control workflow and temporarily store and manipulate data



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.core.manager.ManagerStats
   pyatoa.core.manager.Manager




.. py:exception:: ManagerError

   Bases: :py:obj:`Exception`

   A class-wide custom exception raised when functions fail gracefully


.. py:class:: ManagerStats

   Bases: :py:obj:`dict`

   A simple dictionary that can get and set keys as attributes and has a
   cleaner looking print statement, used for storing internal statistics
   in the Manager class

   .. py:method:: __setattr__(key, value)

      Implement setattr(self, name, value).


   .. py:method:: __getattr__(key)


   .. py:method:: __str__()

      Return str(self).


   .. py:method:: reset()

      Convenience function to reset stats to None



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



