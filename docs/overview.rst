Overview
==============

Pyatoa is a waveform assessment package designed to compare one set of wiggles
to another set of (similar) wiggles. Core functionalities of the package are:

- Waveform standardization and preprocesing
  (with `Obspy <https://github.com/obspy/obspy/>`__)
- Time windowing (with `Pyflex <https://krischer.github.io/pyflex/>`__) and
  adjoint source generation (with
  `Pyadjoint <http://krischer.github.io/pyadjoint/>`__)
- Hierarchical data storage of waveforms, metadata, and measurements
  (with `PyASDF <https://seismicdata.github.io/pyasdf/>`__)
- Bulk measurement aggregation and analysis
  (with `Pandas <https://pandas.pydata.org/>`__)
- Waveform and measurement plotting
- Built-in interface with the automated workflow tool:
  `SeisFlows <https://github.com/adjtomo/seisflows>`__

**In a nutshell**: Feed Pyatoa two waveforms, metadata and some configuration
parameters, it spits out misfit windows and adjoint sources. Feed it *more*
waveforms and metadata, it spits out more windows and adjoint sources, while
also neatly organizing the input data for later bulk analysis.

Why is Pyatoa necessary?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Q**: If Pyatoa simply wraps a number of already-working dependencies, why
do we need it?

**A**: Pyatoa acts as a replacement for the inevitable scripts a User would
have to write to link the data collection capabilities of ObsPy with the
windowing functions of Pyflex with the adjoint source generation of Pyadjoint.

It also provides internal checks along the way to make sure you're not feeding
it garbage, and provides useful routines for performing and understanding
seismic inversions.

What can't Pyatoa?
~~~~~~~~~~~~~~~~~~

Pyatoa cannot generate sythetic waveforms, interface with HPC systems, or
gather seismic data. These capabilities are supplemented by other packages,
some of which are hosted within the `adjTomo <https://github.com/adjtomo/>`__
ecosystem.


How do (I use) Pyatoa?
~~~~~~~~~~~~~~~~~~~~~~~

Pyatoa was written following the design philosophy of ObsPy i.e., called
through a Python interface or within Python scripts, Jupyter notebooks or other
Python packages.

In other words:

.. code:: python
    from pyatoa import Manager

    mgmt = Manager()


Check out `First Glance <first_glance.html>`__ for a short code-based
introduction to Pyatoa.
    
