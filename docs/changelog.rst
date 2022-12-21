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

- Removed all FDSN gathering routines from Pyatoa completely to keep the package
  more lightweight. This functionality has moved to PySEP.

- Shifted all SPECFEM-based I/O routines to PySEP (e.g., reading synthetics, 
  cmtsolutions). This functionality has moved to PySEP

- Removed tests and documentation related to the above 

- Added PySEP as a dependency of Pyatoa

- Overhauled documentation to minimize use of Jupyter notebooks and instead
  hardcode code snippets. This leads to cleaner, more manageable code. Also
  attempted to slim down documentation (read: less wordy)

- Added example data reading functions

- Bombed out pyatoa/scripts repository which contained unused or old scripts

- `Config` class dropped seisflows parameter file reading functions as these 
  were not used

- `Config` class dropped 'start_pad' and 'end_pad' parameters which were tied 
  in to the gathering functionality

- Default preprocessing function changed inputs from Manager class to ObsPy
  objects and soe flags to make it more general

- Cleaned up logging messages by shortening overall log messages to make it
  easier for users to determine what is going on (too wordy before)

- Changes install procedure from the old-style 'setup.py' to a Conda-based
  'environment.yml' file. Also introduces the new-style 'pyproject.toml' file
  for Pip/PyPi
