PYATOA
------

``Pyatoa`` is Python’s Adjoint Tomography Operations Assistant, a
Python-based toolbox for misfit assessment in seismic tomography.
``Pyatoa`` is built to tackle the adjoint tomography problem by
addressing data gathering, waveform standardization and preprocessing,
waveform time-windowing, misfit measurement, misfit assesment, and
visualization.

--------------

The source code of can be found here: https://github.com/bch0w/pyatoa

--------------

``Pyatoa`` is built ontop of the following Python packages:

| `Obspy: <https://github.com/obspy/obspy/wiki>`__ for seismic data
  fetching, handling, processing and organization.
| `Pyflex: <https://krischer.github.io/pyflex/>`__ a Python port of
  Flexwin, an automatic time window selection algorithm.
| `Pyadjoint: <http://krischer.github.io/pyadjoint/>`__ a package for
  calculating misfit and creating adjoint sources.
| `PyASDF: <https://seismicdata.github.io/pyasdf/>`__ heirarchical data
  storage for seismic data.
| `Pandas: <https://pandas.pydata.org/>`__ data analysis and
  manipulation tool for misfit assessment.

--------------

Installation
------------

| Pyatoa is still currently under development so a package manager
  install has not been implemented.
| Install will have to be through the Github repository:

.. code:: bash

   $ git clone https://github.com/bch0w/pyatoa.git
   $ cd pyatoa
   $ python setup.py install

--------------

Doc Contents
------------

.. toctree::
   :maxdepth: 2

   config
   gathering
   manager
   storage
   

