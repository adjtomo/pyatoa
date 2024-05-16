# CHANGELOG

## v0.4.0

>__Note__: Mainly improvements to the Inspector class and its ability to
> generate 'reports' that automatically summarize SeisFlows inversion results.


- #41 Made Manager window selection function calls more explicit; removed 
  simple detrend from preprocessing steps which adversely affected
  non-tapered data
- #42 Introduce `Inspector.generate_report()` function; introduce `revalidate`
  feature in Manager.window() that re-calculates window criteria for updated
  waveforms and drops windows if they do not meet the original criteria; 
  improves `comp_wave` function to plot waveform updates through inversion
- #43 Bugfixes `trim_streams` failing silently causing waveform data to not be
  trimmed but returned as if they were
- #44 Further improves `Inspector.generate_report()` functionality

## v0.3.0

>__Note__: The motivation behind the changes in v0.3.0 were that the original 
> data gathering setup used by Pyatoa was very abstract, opaque, and 
> unncessarily rigid (e.g.,, building path strings out of various components of
> filenames and internal attributes). The new approach to data gathering is to
> use PySEP to perform all data gathering once-and-for-all, including one time
> tasks like instrument removal. The resulting SAC files can then be read in 
> with ObsPy and directly fed into the Manager class for misfit quantification.
> This also gives the User much more control over their data gathering without
> getting confused by Pyatoa's internal data gathering system. 

- Removed ``pyatoa.core.gatherer.Gatherer`` class from package entirely, all 
  data gathering capabilities have been migrated to PySEP, Pyatoa will now only 
  accept input data as already-defined ObsPy objects
- Removed Gatherer-related tests and documentation from package
- Removed ``paths`` attribute from ``pyatoa.core.config.Config`` and all 
  references to the paths attribute throughout the package as these were only
  accessed by the now removed ``Gatherer`` class
- Changed Pyflex and Pyadjoint configuration building procedure in
  ``pyatoa.core.config.Config`` as it was previously abstracted behind a few 
  unncessary functions. ``Config`` now accepts parameters ``pyflex_parameters``
  and ``pyadjoint_parameters`` (dictionaries) that overwrite default Config
  parameters in the underlying Config objects
- Changed ``pyatoa.core.manager.Manager.write()`` to ``write_to_dataset`` to be
  clearer in explaning it's role
- Exposed the default preprocessing procedures directly in the
  ``Manager.preprocess`` function, rather than having it hidden behind a 
  function call to a utility script. Users who want to overwrite the  
  preprocessing need only skip the call to preprocess and perform their own
  tasks on the internally defined ``st_obs`` and ``st_syn`` attributes.
- Removed ``pyatoa.core.manager.Manager``'s ability to save to ASDFDataSet mid
  workflow (i.e., during window and measure). Manager must now use the 
  ``write_to_dataset`` function if it wants to save data to an ASDFDataSet
- Removed the ``pyatoa/plugins`` directory which only contained the pyflex
  preset dictionaries. These were not very flexible, instead they have been
  converted to a docs page for easier accessibility.
- Created Docs page for Pyflex presets that can be copy-pasted into misfit 
  quantification routines
- Added a ``plt.close('all')`` to the end of the Manager's plot routine as
  as a final precaution against leaving an excessive number of Matplotlib 
  figures open
- Overhauled ``pyatoa.core.manager.Manager.flow_multiband`` to mimic behavior 
  the standard behavior of ``Manager.flow``, that is: return internal attributes
  ``windows`` and ``adjsrcs`` which are component-wise dictionaries that each
  contain Pyflex Windows and Pyadjoint AdjointSource objects, respectively. 
  Previously this function returned dictionaries of dictionaries which needed 
  to be further manipulated, now the function averages all adjoint sources 
  from all period bands, and also collects all windows.
- Adjusted and fixed tests based on all the above changes.

## v0.2.2

- Bugfix: Gatherer attempting to access a removed Config parameter
- Resolve PyPDF2 -> PyPDF dependency deprecation warning
- Bugfix: Manager.standardize() only resamples if required, otherwise small time shifting is introduced (Issue \#34)

## v0.2.1

- Updated internal call structures to deal with Pyadjoint v0.2.1 API changes
- Changed internal test ASDFDataSet and created a script to generate new dataset
  because the old one had no way of being remade.
- New Docs + Example + Example data: Processing data with Pyatoa and MPI
- Remove GitHub Pip install links for PySEP, Pyflex and Pyadjoint

## v0.2.2

- Bugfix: Gatherer attempting to access a removed Config parameter
- Resolve PyPDF2 -> PyPDF dependency deprecation warning
