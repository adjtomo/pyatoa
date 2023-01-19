:py:mod:`pyatoa.utils.form`
===========================

.. py:module:: pyatoa.utils.form

.. autoapi-nested-parse::

   Formatting functionality

   Pyatoa relies on data structure being ordered and consistent throughout all the
   various bits of data required, as well as a few standardized string formatters
   to keep everything playing nice. Functions here will aid in reshaping data
   into the correct formats.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.form.format_iter
   pyatoa.utils.form.format_step
   pyatoa.utils.form.format_event_name
   pyatoa.utils.form.convert_stations_to_seed



.. py:function:: format_iter(iteration)

   Format the iteration to be used in internal path naming and labelling.
   Standard is to format with a leading 'i' character followed by two digits.
   Inputs can be strings or integers. Assumes iterations won't go above 99.

   :type iteration: str or int
   :param iteration: input model number, e.g. 0, '0', 'i0', 'i00'
   :rtype: str
   :return: formatted model number, e.g. 'i00', or None if no matching format


.. py:function:: format_step(count)

   Same as for iteration but step count is formatted with a leading 's'

   :type count: str or int
   :param count: input model number, e.g. 0, '0', 's0', 's00'
   :rtype: str
   :return: formatted model number, e.g. 's00', or None if no matching format


.. py:function:: format_event_name(ds_or_event)

   Formalize the definition of Event ID in Pyatoa

   :type ds_or_event: pyasdf.ASDFDataSet or obspy.core.event.Event or str
   :param ds_or_event: get dataset event name from the event resource_id
   :rtype: str
   :return: the event name to be used for naming schema in the workflow


.. py:function:: convert_stations_to_seed(stations_file='./STATIONS', path_to_save_seed='./seed', **kwargs)

   A convenience function to format a SPECFEM file into SeisFlows3 ready files.
   Convert a SPECFEM STATIONS file into a directory of SEED formatted
   StationXML files which are REQUIRED by a Pyatoa + SeisFlows3 workflow.
   Kwargs are passed to pyatoa.utils.write.write_inv_seed()
   See above function for available options on writing SEED StationXML files

   :type stations_file: str
   :param stations_file: path to the STATIONS file defined in SPECFEM format
   :type path_to_save_seed: str
   :param path_to_save_seed: path to save the output SEED formatted
       StationXML files


