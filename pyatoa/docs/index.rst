Python’s
~~~~~~~~

Adjoint Tomography Operations Assistant
---------------------------------------

``Pyatoa`` is a misfit quantification toolbox for the modern
tomographer. It provides abstraction over key Python packages to
facilitate data gathering, preprocessing, misfit analysis and
visualization in a seismic tomography problem. The design of ``Pyatoa``
was inspired by ``Obspy``, translating to a design ethos of rapid
development, through scripting or shell interaction, and useful
abstractions to accomplish tasks in tomography.

The source code of ``Pyatoa`` can be found here:
https://github.com/bch0w/pyatoa

--------------

``Pyatoa`` is built around a handful of Python packages:

| `Obspy: <https://github.com/obspy/obspy/wiki>`__ for seismic data
  fetching, handling, processing and organization.
| `Pyflex: <https://krischer.github.io/pyflex/>`__ a Python port of
  Flexwin, an automatic time window selection algorithm.
| `Pyadjoint: <http://krischer.github.io/pyadjoint/>`__ a package for
  calculating misfit and creating adjoint sources.
| `PyASDF: <https://seismicdata.github.io/pyasdf/>`__ heirarchical data
  storage for seismic data.
| `Matplotlib: <https://matplotlib.org/>`__ 2D plotting library for
  visualization of waveforms, statistics, misfit etc
| `Basemap: <https://matplotlib.org/basemap/>`__ A mapping library for
  source receiver distributions, raypaths, etc. (deprecated, Cartopy in
  the future).

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

Doc Contents:
------------------
.. toctree:: 
    :maxdepth: 2  

    tutorials  

