:py:mod:`pyatoa.utils.asdf.add`
===============================

.. py:module:: pyatoa.utils.asdf.add

.. autoapi-nested-parse::

   ASDF Datasets can be given auxiliary data to supplement the existing waveform,
   event and station information contained. The functions contained in this script
   add new auxiliary data structures to existing ASDF datasets



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.asdf.add.add_misfit_windows
   pyatoa.utils.asdf.add.add_adjoint_sources



.. py:function:: add_misfit_windows(windows, ds, path)

   Write Pyflex misfit windows into the auxiliary data of an ASDFDataSet

   :type windows: dict of list of pyflex.Window
   :param windows: dictionary of lists of window objects with keys
       corresponding to components related to each window
   :type ds: pyasdf.ASDFDataSet
   :param ds: ASDF data set to save windows to
   :type path: str
   :param path: internal pathing to save location of auxiliary data


.. py:function:: add_adjoint_sources(adjsrcs, ds, path, time_offset)

   Writes the adjoint source to an ASDF file.

   .. warning::
       It is inherently assumed SPECFEM will be using the adjoint source

   .. note::
       Modified from Pyadjoint source code:
       pyadjoint.adjoint_source.write_to_asdf()

   .. note::
       SPECFEM requires one additional parameter: the temporal offset of the
       first sample in seconds. This will have been set by the parameter
       USER_T0 in the constants.h.in file of SPECFEM's setup directory

   :type adjsrcs: list of pyadjoint.asdf_data_set.ASDFDataSet
   :param adjsrcs: adjoint source to save
   :type ds: pyasdf.asdf_data_set.ASDFDataSet
   :type path: str
   :param path: internal pathing for save location in the auxiliary data attr.
   :type ds: pyasdf.ASDFDataSet
   :param ds: The ASDF data structure read in using pyasdf.
   :type time_offset: float
   :param time_offset: The temporal offset of the first sample in seconds.
       This is required if using the adjoint source as input to SPECFEM.


