===================================================
Python’s Adjoint Tomography Operations Assistant
===================================================
``Pyatoa`` is a Python-based toolbox meant to facilitate waveform comparisons in 
adjoint tomography. With humble origins as a disjointed collection of scripts, 
it has grown into a cohesive library for misfit quantification with
purpose-built objects to facilitate automation of common tasks in full waveform
tomography.

.. figure:: images/data-synthetic_misfit.png
    :alt: An example waveform figure showing off some of Pyatoa's features

    An example of ``Pyatoa``'s waveform comparison capabilities. Observed 
    (black) and synthetic (red) waveforms compared within
    time windows (orange boxes), culminating in adjoint sources (green)

----------

This docs page provides an overview of Pyatoa through introductory material,
in-depth tutorials, and API for core classes and supporting utilities. Have a 
look at the :doc:`Quick Start </quickstart>` and 
:doc:`Gallery </gallery>` pages to get an introductory overview of the package.

Source code can be found on GitHub: https://github.com/bch0w/pyatoa


--------------

Installation
~~~~~~~~~~~~

``Pyatoa`` is in ongoing development so package manager installation is not
currently available, however pip and conda install procedures are on the to-do
list. For now, install should be accomplished via pip and the GitHub repository.

It is recommended that ``Pyatoa`` be installed inside a Conda environment to
preserve your root environment. The 'devel' branch provides the latest codebase.

.. code:: bash

   $ conda create -n pyatoa python=3.7
   $ conda activate pyatoa
   $ git clone --branch devel https://github.com/bch0w/pyatoa.git
   $ cd pyatoa
   $ pip install .


Running Tests
~~~~~~~~~~~~~

Tests ensure ``Pyatoa`` runs as expected, these require installation of Pytest.
If any changes are made to the source code, please run tests to ensure nothing
is broken.

.. code:: bash

   $ cd pyatoa/tests
   $ pytest

--------------

Dependencies
~~~~~~~~~~~~

Credit where credit is due, ``Pyatoa`` is a high-level purpose-buiolt wrapper 
of the following packages:

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

How to Cite
~~~~~~~~~~~~

If you use ``Pyatoa``, consider citing the related publication:
`Chow et al. (2020) <https://academic.oup.com/gji/article/223/3/1461/5897358>`__.

Chow, B., Kaneko, Y., Tape, C., Modrak, R., & Townend, J. (2020).   
An automated workflow for adjoint tomography — waveform misfits and synthetic inversions for the North Island, New Zealand.  
Geophysical Journal International, 223(3), 1461-1480.


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Introduction

   overview
   quickstart
   gallery

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Tutorials

   datasyn_misfit
   inversion_prep

.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: Core Functionality

   config
   manager
   gatherer
   pyaflowa
   storage
   inspector
   utilities

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Interface

   standards
   logging

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   code_dev_plan
   changelog

