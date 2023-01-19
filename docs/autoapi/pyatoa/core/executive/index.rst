:py:mod:`pyatoa.core.executive`
===============================

.. py:module:: pyatoa.core.executive

.. autoapi-nested-parse::

   A class and functionality to direct multiple Managers to parallelize data
   processing in Pyatoa using the in-built concurrent.futures package

   .. rubric::
       >>> from pyatoa import Executive
       >>> exc = Executive(event_ids=..., station_codes=..., config=...)
       >>> misfits = exc.process()

   .. note::
       I was getting the following error randomly on my runs, I think it had to
       do with the fact that I was requesting all of my cores to do the job (16).
       Dropped the max number of jobs to 8 and things ran fine. Will investigate

       concurrent.futures.process.BrokenProcessPool:
       A process in the process pool was terminated abruptly while the future
       was running or pending.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.core.executive.Executive




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



