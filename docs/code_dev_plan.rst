Code Development Plan
======================

Based on the `Computational Infrastructure for Geodynamics
<https://geodynamics.org/>`__ (CIG) `Software Development Best Practices 
<https://github.com/geodynamics/best_practices/blob/master/
SoftwareDevelopmentBestPractices.md>`__ we propose the following (yearly updated)
Code Development Plan (CDP) which lays out prioritization of new features and 
estimated timetables for their implementation. 

The goals listed here are target dates and are by no means hard deadlines.
Rather, they provide the user-base and community with an idea of where package
development is headed for Pyatoa, and what new features may be expected for
the package. If the community identifies features that would be useful to 
include but are not listed in the CDP, they are free to open up a GitHub issue
or email the developers and we will consider adding to the CDP.

It is our intention to work to meet all of the `Standard Best Practices 
<https://github.com/geodynamics/best_practices/blob/master/
SoftwareDevelopmentBestPractices.md#standard-best-practices>`__
in the CIG Best Practices guideline, and aim for `Target Best Practices 
<https://github.com/geodynamics/best_practices/blob/master/
SoftwareDevelopmentBestPractices.md#target-best-practices>`__ where
possible.

----------------------

CDP 2022 -- 2023 
----------------
[Last updated 3/3/22] This CDP outlines the goals for the year 2022.


Short-Term Priorities [1--2 months]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Package manager installation with Pip and/or Conda (no need for cloning repo) 
- Pre-commit framework (Flake8, Black, pre-commit, codecoverage) 
- Travis Continuous Integration
- Update all dependencies to latest release, address deprecation warnings
- Unit tests for plotting utilities
- Integration tests
- Code version citation (Zenodo)
- FAQs and user knowledge base

Medium-Term Objectives [3--6 months]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Cookbook examples and `ObsPy-like gallery 
  <https://docs.obspy.org/gallery.html>`__ for new user introduction
- Developer documentation for extending code in anticipated ways
- Guidelines for parameter scales/combinations for which code is designed/tested
- Flowchart or diagrammatic illustrations of Pyatoa core objects.
- Remove mapping dependency on Basemap (now deprecated), and either shift to
  PyGMT, Cartopy, Matplotlib, separate repo, or drop altogether
- Pyflex and Pyadjoint dependencies currently point directly to GitHub repos 
  which limits portability. Figure out best way to actively maintain these
  dependencies (options: hard fork to Pyatoa org page and host on Pip? 
  Copy entire codebase directly into Pyatoa?)
- Introduce more flexibility in handling FORCESOLUTION and SOURCE files. 
  Current implementation was built around CMTSOLUTION files only.
- Introduce flexibility in handling different types of temporary data. 
  E.g., flagging data which has already been preprocessed, or had instrument 
  response removed
- Add operator functions (e.g., +=, -=, /=) to core classes where relevant
  (similar to ObsPy Stream() += Stream())
- More emphasis on location codes, which are currently not fully integrated
  into bookkeeping.


Long-Term Goals [> 6 months]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Integrate `Meshfem3D Python-based meshing tool 
  <https://github.com/bch0w/simutils/blob/master/meshing/meshfem/
  prepare_meshfem.py>`__ into pyatoa.core 
- Integrate event declustering scripts into utilities 
  (see Amanda McPherson, Carl Tape)
- Standalone SPECFEM2D kernel creation examples with adjoint source and 
  auxiliary file generation handled by Pyatoa. 
- Improved workflow integration with 
  `SeisFlows3 <https://github.com/bch0w/seisflows3>`__
  (e.g., passing objects in memory and not on disk, actual parallel 
  implementation as opposed to embarrasingly parallel)
- Enhance Pyatoa's simulation domain recognition abilities 
  (e.g., defining domain border based on longest wavelength, station distance)
- Improved parallel capacity for waveform processing (with concurrent.futures)
- Reduce filesize overhead during inversions (e.g., stop storing redundant 
  Station metadata in ASDFDataSets, utilize HDF5 file compression)
- Checkpointing or failsafes to safeguard against HDF5 file corruption

