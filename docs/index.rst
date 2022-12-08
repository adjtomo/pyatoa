===================================================
Python’s Adjoint Tomography Operations Assistant
===================================================
`Pyatoa <https://github.com/adjtomo/pyatoa>`__ is an open-source, Python-based
misfit quantification toolbox for full waveform inversion / adjoint tomography.

Primarily intended to calculate misfit between two waveforms, it also hosts
tools for data visualization, data storage, and statistical analysis of misfit
for seismic inversions run with
`SeisFlows <https://github.com/adjtomo/seisflows>`__.

Pyatoa is hosted on `GitHub <https://github.com/adjtomo/pyatoa>`__ as part of
the `adjTomo organization <https://github.com/adjtomo>`__.

.. figure:: images/data-synthetic_misfit.png
    :alt: An example waveform figure showing off some of Pyatoa's features

    An example of Pyatoa's waveform comparison capabilities. Observed 
    (black) and synthetic (red) waveforms compared within
    time windows (orange boxes), culminating in adjoint sources (green)

----------

Quickstart
~~~~~~~~~~~

- Have a look at the `Overview <overview.html>`__ page to learn about Pyatoa.
- The `Getting Started <getting_started.html>`__ page shows you how to run a
  small misfit quantification example.
- The `Gallery <gallery.html>`__ provides a number of figures which illustrate
  how Pyatoa can be used to assess misfit during a seismic inversion.

--------------

Installation
~~~~~~~~~~~~

It is recommended that Pyatoa be installed inside a Conda environment.
The 'devel' branch provides the latest codebase.

.. code:: bash

   git clone --branch devel https://github.com/bch0w/pyatoa.git
   cd pyatoa
   conda env create -f environment.yml
   conda activate pyatoa


Running Tests
`````````````

Tests ensure Pyatoa runs as expected after changes are made to the source code.
You can run tests with Pytest

.. code:: bash

   $ cd pyatoa/tests
   $ pytest

--------------

Cite Pyatoa
~~~~~~~~~~~~

If you use Pyatoa in your own research, please consider citing the related
publication:
`Chow et al. (2020) <https://academic.oup.com/gji/article/223/3/1461/5897358>`__.

    Chow, B., Kaneko, Y., Tape, C., Modrak, R., & Townend, J. (2020).
    *An automated workflow for adjoint tomography — waveform misfits and synthetic
    inversions for the North Island, New Zealand.*
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

.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: Core Functionality

   config
   manager
   gatherer
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

