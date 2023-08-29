Changelog
=========

v0.3.0
------

   **Note**: The motivation behind the changes in v0.3.0 were that the
   original data gathering setup used by Pyatoa was very abstract,
   opaque, and unncessarily rigid (e.g.,, building path strings out of
   various components of filenames and internal attributes). The new
   approach to data gathering is to use PySEP to perform all data
   gathering once-and-for-all, including one time tasks like instrument
   removal. The resulting SAC files can then be read in with ObsPy and
   directly fed into the Manager class for misfit quantification. This
   also gives the User much more control over their data gathering
   without getting confused by Pyatoa’s internal data gathering system.

-  Removed ``pyatoa.core.gatherer.Gatherer`` class from package
   entirely, all data gathering capabilities have been migrated to
   PySEP, Pyatoa will now only accept input data as already-defined
   ObsPy objects
-  Removed Gatherer-related tests and documentation from package
-  Removed ``paths`` attribute from ``pyatoa.core.config.Config`` and
   all references to the paths attribute throughout the package as these
   were only accessed by the now removed ``Gatherer`` class
-  Changed Pyflex and Pyadjoint configuration building procedure in
   ``pyatoa.core.config.Config`` as it was previously abstracted behind
   a few unncessary functions. ``Config`` now accepts parameters
   ``pyflex_parameters`` and ``pyadjoint_parameters`` (dictionaries)
   that overwrite default Config parameters in the underlying Config
   objects
-  Changed ``pyatoa.core.manager.Manager.write()`` to
   ``write_to_dataset`` to be clearer in explaning it’s role
-  Exposed the default preprocessing procedures directly in the
   ``Manager.preprocess`` function, rather than having it hidden behind
   a function call to a utility script. Users who want to overwrite the
   preprocessing need only skip the call to preprocess and perform their
   own tasks on the internally defined ``st_obs`` and ``st_syn``
   attributes.
-  Removed ``pyatoa.core.manager.Manager``\ ’s ability to save to
   ASDFDataSet mid workflow (i.e., during window and measure). Manager
   must now use the ``write_to_dataset`` function if it wants to save
   data to an ASDFDataSet
-  Removed the ``pyatoa/plugins`` directory which only contained the
   pyflex preset dictionaries. These were not very flexible, instead
   they have been converted to a docs page for easier accessibility.
-  Created Docs page for Pyflex presets that can be copy-pasted into
   misfit quantification routines
-  Added a ``plt.close('all')`` to the end of the Manager’s plot routine
   as as a final precaution against leaving an excessive number of
   Matplotlib figures open
-  Overhauled ``pyatoa.core.manager.Manager.flow_multiband`` to mimic
   behavior the standard behavior of ``Manager.flow``, that is: return
   internal attributes ``windows`` and ``adjsrcs`` which are
   component-wise dictionaries that each contain Pyflex Windows and
   Pyadjoint AdjointSource objects, respectively. Previously this
   function returned dictionaries of dictionaries which needed to be
   further manipulated, now the function averages all adjoint sources
   from all period bands, and also collects all windows.
-  Adjusted and fixed tests based on all the above changes.

v0.2.2
------

-  Bugfix: Gatherer attempting to access a removed Config parameter
-  Resolve PyPDF2 -> PyPDF dependency deprecation warning
-  Bugfix: Manager.standardize() only resamples if required, otherwise
   small time shifting is introduced (Issue #34)

v0.2.1
------

-  Updated internal call structures to deal with Pyadjoint v0.2.1 API
   changes
-  Changed internal test ASDFDataSet and created a script to generate
   new dataset because the old one had no way of being remade.
-  New Docs + Example + Example data: Processing data with Pyatoa and
   MPI
-  Remove GitHub Pip install links for PySEP, Pyflex and Pyadjoint

v0.2.0
------
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
