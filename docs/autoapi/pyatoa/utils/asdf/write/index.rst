:py:mod:`pyatoa.utils.asdf.write`
=================================

.. py:module:: pyatoa.utils.asdf.write

.. autoapi-nested-parse::

   Functions to convert ASDFDataSet into individual data files that can be
   stored in a directory structure. Not very smart, just dumps all data into a
   given path.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.asdf.write.write_all
   pyatoa.utils.asdf.write.write_events
   pyatoa.utils.asdf.write.write_stations
   pyatoa.utils.asdf.write.write_waveforms
   pyatoa.utils.asdf.write.write_windows
   pyatoa.utils.asdf.write.write_adjoint_sources



.. py:function:: write_all(ds, path='./')

   Convenience function to dump everything inside a dataset

   :type ds: pyasdf.ASDFDataSet
   :param ds: Dataset containing info the write
   :type path: str
   :param path: path to save data to


.. py:function:: write_events(ds, path='./')

   Write Event object as a QuakeML file.

   :type ds: pyasdf.ASDFDataSet
   :param ds: Dataset containing info the write
   :type path: str
   :param path: path to save data to


.. py:function:: write_stations(ds, path='./')

   Write Stations dataless files as STATIONXML files

   :type ds: pyasdf.ASDFDataSet
   :param ds: Dataset containing info the write
   :type path: str
   :param path: path to save data to


.. py:function:: write_waveforms(ds, path='./', station=None, tag=None, format='MSEED')

   Write waveforms as MSEED files

   :type ds: pyasdf.ASDFDataSet
   :param ds: Dataset containing info the write
   :type path: str
   :param path: path to save data to


.. py:function:: write_windows(ds, path='./')

   Write MisfitWindows as .JSON files

   :type ds: pyasdf.ASDFDataSet
   :param ds: Dataset containing info the write
   :type path: str
   :param path: path to save data to


.. py:function:: write_adjoint_sources(ds, path='./')

   Write AdjointSources as ASCII files in directories corresponding to model
   number and step count.

   :type ds: pyasdf.ASDFDataSet
   :param ds: Dataset containing info the write
   :type path: str
   :param path: path to save data to


