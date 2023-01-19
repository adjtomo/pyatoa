:py:mod:`pyatoa.core.inspector`
===============================

.. py:module:: pyatoa.core.inspector

.. autoapi-nested-parse::

   A class to aggregate time windows, source-receiver information and misfit
   using Pandas.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.core.inspector.Inspector




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



