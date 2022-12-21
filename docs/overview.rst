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

**In a nutshell**: Provide two similar waveforms and some configuration
parameters, receive misfit windows and adjoint sources. Do this in bulk and
analyze all returned misfit in aggregate.

Why is Pyatoa necessary?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If Pyatoa simply wraps a number of already-working dependencies, why
do we need it?

Pyatoa provides an easy-to-use tool that acts in liue of the inevitable 
collection scripts a User would have to write to link the data collection 
capabilities of ObsPy with the windowing functions of Flexwin/Pyflex with the 
adjoint source generation of Pyadjoint.

It also provides structure for these processes, to ensure that data-synthetic 
misfit is calculated properly. Similarly, useful routines within the package 
provide assistance during seismic inversions.

What doesn't Pyatoa do?
~~~~~~~~~~~~~~~~~~~~~~~

Pyatoa does not: generate sythetic waveforms, interface with HPC systems,
gather seismic data. These capabilities are supplemented by other packages,
some of which are hosted within the `adjTomo <https://github.com/adjtomo/>`__
ecosystem.

Recommended tools to use with Pyatoa:

- Gather seismic data: `PySEP <https://github.com/adjtomo/pysep>`__
- Generate synthetics: `SPECFEM <https://github.com/specfem>`__
- Automate inversions on HPCs: `SeisFlows <https://github.com/adjtomo/seisflows>`__


How do (I use) Pyatoa?
~~~~~~~~~~~~~~~~~~~~~~~

Pyatoa was written following the design philosophy of ObsPy i.e., called
through a Python interface or within Python scripts, Jupyter notebooks or other
Python packages. That is:

.. code:: python

    from pyatoa import Manager

    mgmt = Manager()


Check out `First Glance <first_glance.html>`__ for a short code-based
introduction to Pyatoa.
