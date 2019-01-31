.. Pyatoa documentation master file, created by
   sphinx-quickstart on Wed Jan 23 14:02:21 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyatoa 
======
Adjoint Tomography Operational Assistance
_________________________________________
``Pyatoa`` is a workflow management tool designed for the adjoint tomography problem. Its purpose is to provide an easy to use command line or scripting interface for intermediate steps in the tomography problem. These steps include data retrieval, standardization, preprocessing, visualization. 
``Pyatoa`` relies heavily on ObsPy_ for all seismic trace processing and manipulation. It also calls upon Pyflex_ and Pyadjoint_ for misfit quantification and adjoint source creation, and PyASDF_ for I/O.

.. _Pyflex: http://krischer.github.io/pyflex/
.. _Pyadjoint: http://krischer.github.io/pyadjoint/
.. _PyASDF: http://seismicdata.github.io/pyasdf/
.. _ObsPy: https://github.com/obspy/obspy/wiki
.. toctree::
   tutorial 
   :maxdepth: 2
   :caption: Contents:


