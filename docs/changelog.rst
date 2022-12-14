Change Log
==============

Version 0.2.0
~~~~~~~~~~~~~~~
- Renamed 'Quickstart' doc to 'A short example', created a new 'Quickstart' doc which has a short code snippet that creates a figure.

- Revamped documentation, switched to new style of building documentation using only .rst files (rather than building off of Jupyter notebooks directly in RTD, which was high in memory consumption)

- Switched API to sphinx-autoapi as opposed to autodoc (moves load to local side rather than RTD generating API)

- Added new hard requirement for pillow>=8.4.0 for image manipulation

- Inspector now sets aspect ratios on map-like plots based on latitude. Rough equation but makes things like maps and raypath plots scale better

- Added Inspector, Gallery docs pages. Included some docs statements separating Pyatoa and SeisFlows3

- Added test data pertaining to docs. Docs now do not work directly with test data but rather make copies in non-repo'd directories. 

- Added __init__.py to all relevant directories that were missing it before

- Removed 'Pyaflowa' class, tests and related docs. Pyaflowa was the 
  Pyatoa-SeisFlows3 interaction class. All of its functionality has been moved
  directly into the `seisflows.preprocess.pyaflowa.Pyaflowa` class

- Added hard requirements for Cartopy and Proj in requirements.txt as their 
  absence was causing some dependency conflicts

Version 0.3.0
~~~~~~~~~~~~~~

- Removed all FDSN gathering routines from Pyatoa completely to keep the package
  more lightweight. This functionality has moved to PySEP.

- Shifted all SPECFEM-based I/O routines to PySEP (e.g., reading synthetics, 
  cmtsolutions). This functionality has moved to PySEP

- Removed tests and documentation related to the above 

- Added PySEP as a dependency of Pyatoa

- Overhauled documentation