Pyatoa
------

**[ Python’s Adjoint Tomography Operations Assistant ]**


``Pyatoa`` is a Python-based toolbox for misfit assessment in adjoint
tomography. It began as a loose collection of scripts, and has grown into a 
cohesive set of tools for: 
seismic data gathering, waveform processing, time-windowing, misfit measurement, 
measurement aggregation, and visualization.

This docs page provides an overview of intended use through 
quickstart introductions, in-depth tutorials, and API for each of Pyatoa's core 
classes and supporting utilities.

Source code can be found on Github: https://github.com/bch0w/pyatoa

--------------

Installation
~~~~~~~~~~~~

``Pyatoa`` is in ongoing development so package manager installation is not 
currently available. Install should be accomplished using pip.  

It is recomended that Pyatoa be installed inside a Conda environment to
preserve your root environment. The 'devel' branch provides the latest codebase.


.. code:: bash

   $ conda create -n pyatoa python=3.7
   $ source activate pyatoa
   $ git clone https://github.com/bch0w/pyatoa.git
   $ cd pyatoa
   $ git checkout devel
   $ pip install -r requirements.txt .


Running Tests
~~~~~~~~~~~~~

Tests ensure ``Pyatoa`` runs as expected, these require installation of Pytest.

.. code:: bash

   $ pip install pytest
   $ cd pyatoa/tests
   $ pytest

--------------

Dependencies
~~~~~~~~~~~~

Credit where credit is due, ``Pyatoa`` is a high-level API relying on the following:

-  `Obspy: <https://github.com/obspy/obspy/wiki>`__ for seismic
   data handling, data gathering, and waveform processing
-  `Pyflex: <https://krischer.github.io/pyflex/>`__ for automatic time
   window selection of time-series data
-  `Pyadjoint: <http://krischer.github.io/pyadjoint/>`__ for calculation of
   waveform misfit and generation of adjoint sources
-  `PyASDF: <https://seismicdata.github.io/pyasdf/>`__ for efficient,
   hierarchical, self-describing storage of seismic data
-  `Pandas: <https://pandas.pydata.org/>`__ to simplify data aggregation
   and misfit assessment  

--------------

Cite
~~~~

If you use ``Pyatoa``, consider citing our recent publication:
`Chow et al. (2020) <https://academic.oup.com/gji/article/223/3/1461/5897358>`__.

Chow, B., Kaneko, Y., Tape, C., Modrak, R., & Townend, J. (2020).   
An automated workflow for adjoint tomography — waveform misfits and synthetic inversions for the North Island, New Zealand.  
Geophysical Journal International, 223(3), 1461-1480.


.. toctree::
   :maxdepth: 3
   :hidden:
  
   intro
   tutorials
   standards
   logging
   pyatoa_api
   
