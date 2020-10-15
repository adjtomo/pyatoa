Pyatoa
------

**Python’s Adjoint Tomography Operations Assistant**

--------------

``Pyatoa`` is a Python-based toolbox for misfit assessment in seismic
tomography. ``Pyatoa`` was designed to facilitate automation within the
adjoint tomography workflow by abstracting the following operations:
data gathering, waveform standardization, preprocessing, time-windowing,
misfit measurement, misfit assesment, and visualization.

To accomplish this, ``Pyatoa`` is built ontop of the following Python
packages.

-  `Obspy: <https://github.com/obspy/obspy/wiki>`__ to tackle seismic
   data handling and processing
-  `Pyflex: <https://krischer.github.io/pyflex/>`__ for automatic time
   window selection of waveform data
-  `Pyadjoint: <http://krischer.github.io/pyadjoint/>`__ to calculate
   waveform misfit and generate adjoint sources
-  `PyASDF: <https://seismicdata.github.io/pyasdf/>`__ for efficient,
   hierarchical, self-describing storage of seismic data
-  `Pandas: <https://pandas.pydata.org/>`__ to simplify data aggregation
   and misfit assessment

--------------

Source code
~~~~~~~~~~~

| ``Pyatoa`` can be found on GitHub: https://github.com/bch0w/pyatoa

--------------

Installation
~~~~~~~~~~~~

Install is currently done through cloning the Github repository.

.. code:: bash

   $ git clone https://github.com/bch0w/pyatoa.git
   $ cd pyatoa
   $ python setup.py install

Package manager installation coming soon!

.. toctree::
   :maxdepth: 3
   :hidden:
  
   intro
   tutorials
   standards
   logging
   pyatoa_api
   
