:py:mod:`pyatoa.utils.write`
============================

.. py:module:: pyatoa.utils.write

.. autoapi-nested-parse::

   For writing various output files used by Pyatoa, Specfem and Seisflows



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.write.write_inv_seed
   pyatoa.utils.write.write_misfit
   pyatoa.utils.write.write_stations_adjoint
   pyatoa.utils.write.write_adj_src_to_ascii



.. py:function:: write_inv_seed(inv, path='./', dir_structure='{sta}.{net}', file_template='RESP.{net}.{sta}.{loc}.{cha}', components='ZNE', channel_code='HX{comp}', **kwargs)

   Pyatoa requires stations to be discoverable in SEED format, i.e., in a data
   center repository. This structure dictates that each component of each
   station has its own individual StationXML file, saved in a specific
   directory structure with a unique file naming schema.

   This utility is useful for creating the necessary StationXML files for
   temporary or synthetic stations which are not discoverable via
   FDSN or through datacenters.

   .. note::
       kwargs are passed to obspy.core.inventory.channel.Channel

   :type path: str
   :param path: location to save StationXML files to
   :type dir_structure: str
   :param dir_structure: template for directory structure, likely should not
       need to change from the default
   :type file_template: str
   :param file_template: template for file naming, likely should not change
       from default template
   :type channels: str
   :param channels: OPTIONAL, if inventory does not contain channels (e.g.,
       if read from a SPECFEM STATIONS file), channels will be generated here.
   :type channel_code: str
   :param channel_code: Explicitely defined default channel values for
   generating channels on the fly when none are provided by the inventory


.. py:function:: write_misfit(ds, iteration, step_count=None, path='./', fidout=None)

   This function writes a text file containing event misfit.
   This misfit value corresponds to F_S^T of Eq 6. Tape et al. (2010)

   e.g. path/to/misfits/{iteration}/{event_id}

   These files will then need to be read by: seisflows.workflow.write_misfit()

   :type ds: pyasdf.ASDFDataSet
   :param ds: processed dataset, assumed to contain auxiliary_data.Statistics
   :type iteration: str or int
   :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
   :type step_count: str or int
   :param step_count: step count e.g. "s00". Will be formatted so int ok.
   :type path: str
   :param path: output path to write the misfit. fid will be the event name
   :type fidout: str
   :param fidout: allow user defined filename, otherwise default to name of ds
       note: if given, var 'pathout' is not used, this must be a full path


.. py:function:: write_stations_adjoint(ds, iteration, specfem_station_file, step_count=None, pathout=None)

   Generate the STATIONS_ADJOINT file for Specfem input by reading in the
   STATIONS file and cross-checking which adjoint sources are available in the
   Pyasdf dataset.

   :type ds: pyasdf.ASDFDataSet
   :param ds: dataset containing AdjointSources auxiliary data
   :type iteration: str or int
   :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
   :type step_count: str or int
   :param step_count: step count e.g. "s00". Will be formatted so int ok.
       If NoneType, final step of the iteration will be chosen automatically.
   :type specfem_station_file: str
   :param specfem_station_file: path/to/specfem/DATA/STATIONS
   :type pathout: str
   :param pathout: path to save file 'STATIONS_ADJOINT'


.. py:function:: write_adj_src_to_ascii(ds, iteration, step_count=None, pathout=None, comp_list='ZNE')

   Take AdjointSource auxiliary data from a Pyasdf dataset and write out
   the adjoint sources into ascii files with proper formatting, for input
   into PyASDF.

   .. note::
       Specfem dictates that if a station is given as an adjoint source,
       all components must be present, even if some components don't have
       any misfit windows. This function writes blank adjoint sources
       (an array of 0's) to satisfy this requirement.

   :type ds: pyasdf.ASDFDataSet
   :param ds: dataset containing adjoint sources
   :type iteration: str or int
   :param iteration: iteration number, e.g. "i00". Will be formatted so int ok.
   :type step_count: str or int
   :param step_count: step count e.g. "s00". Will be formatted so int ok.
           If NoneType, final step of the iteration will be chosen automatically.
   :type pathout: str
   :param pathout: path to write the adjoint sources to
   :type comp_list: str
   :param comp_list: component list to check when writing blank adjoint sources
       defaults to N, E, Z, but can also be e.g. R, T, Z


