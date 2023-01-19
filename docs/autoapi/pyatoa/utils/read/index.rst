:py:mod:`pyatoa.utils.read`
===========================

.. py:module:: pyatoa.utils.read

.. autoapi-nested-parse::

   Utilities for reading various file types, mostly from Specfem3D to ObsPy classes
   These are meant to be standalone functions so they may repeat some functionality
   found elsewhere in the package.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.read.read_fortran_binary
   pyatoa.utils.read.read_station_codes



.. py:function:: read_fortran_binary(path)

   Convert a Specfem3D fortran .bin file into a NumPy array,
   Copied verbatim from Seisflows/plugins/solver_io/fortran_binary.py/_read()

   :type path: str
   :param path: path to fortran .bin file
   :rtype: np.array
   :return: fortran binary data as a numpy array


.. py:function:: read_station_codes(path_to_stations, loc='??', cha='*', seed_template='{net}.{sta}.{loc}.{cha}')

   Read the SPECFEM3D STATIONS file and return a list of codes (Pyatoa format)
   that are accepted by the Manager and Pyaflowa classes. Since the STATIONS
   file only provides NET and STA information, the user must provide the
   location and channel information, which can be wildcards.

   :type path_to_stations: str
   :param path_to_stations: full path to the STATIONS file
   :type loc: str
   :param loc: formatting of the location section of the code, defaults to
       '??' two-digit wildcard
   :type cha: str
   :param cha: formatting of the channel section fo the code, defaults to
       'HH?' for wildcard component of a high-gain seismometer. Follows SEED
       convention (see IRIS).
   :type seed_template: str
   :param seed_template: string template to be formatted with some combination
       of 'net', 'sta', 'loc' and 'cha', used for generating station codes
   :rtype: list of str
   :return: list of codes to be used by the Manager or Pyaflowa classes for
       data gathering and processing


