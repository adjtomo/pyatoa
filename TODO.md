## For Version 0.2.0

#### Bugs
- [X] FDSNException in gatherer.gather_observed( Uknown Error (timeout))
      Added catch for FDSNException in gather_observed() and return st_obs=None
- [ ] Pyflex value Error is being thrown (pyflex.window_selector() line 427 np.abs(noise).max() zero size array)
- [X] Pyaflowa exit gracefully if no data is gathered for an entire event, currently finalize is throwing uncaught errors
      Pyaflowa counts successful processes, if none, finalize is skipped
- [X] Fix: waveform plots not deleted if theyre not included in the composite pdf
      Moved purge outside loop and set it to delete all files inside the dir with the correct tag
- [X] Pyaflowa waveform composites aren't made if a thrifty inversion is done because it skips s00
      removed the step count requirement in the if statement

#### General
- [ ] Proper docs for all functions. methods, classes 
- [ ] Finish writing tests with PyTest
- [X] ipynb in gitignore
- [ ] remove large data files from test data
- [ ] make model number, step count formatters standard package wide functions
- [X] check unused kwargs in the config in the case of typos
      Added a check call at the end of check() that puts up a user warning if unnused kwargs fall through
- [ ] include __repr__ for all classes

#### Config
- [ ] Include UTM projection into config and propogate into scripts

#### Manager
- [X] Initialize an empty Manager with an empty Config to remove the need to call Config separately
- [ ] Move window by amplitude in Manager.window() into its own function
- [ ] Check if convolve_stf properly performs the time shift
- [ ] Option to save processed streams in the dataset
- [ ] Change preprocess() function to take a Manager object, that way when the `overwrite` parameter is called, the User knows exactly
      how to write their new preprocess function
- [ ] Fixed windowing might encounter some problems because the synthetic trace is changing, so the values of max_cc_, cc_shift and dlnA       are not being re-evaluated. Can this be remedied?
- [ ] Weighting adjoint sources by station proximity?
- [ ] ability to turn off saving for pyasdf dataset

#### ASDF
- [ ] ASDF export all data to directory structure, saving waveforms, inventories, aux data etc into individual 
- [ ] Generate waveform plots, maps, from an ASDF dataset. As in remove the need to create a Manager just to make 
      waveform plots, if a datset already exists
- [X] Processing provenance saved into auxiliary_data?
      Saved processing stats from obspy stream into the Config object for each model/step
- [X] Retain step count information for MisfitWindows and AdjointSources

#### Misc.
- [X] Move all hardcoded stuff into plugins, e.g. fault_coordinates, geonet_moment_tensor
- [X] Move tools and visualize out of utils dir into main dir

#### Visuals
- [X] Depth cross section of a Catalog object
      In pyaflowa artist 
- [ ] Plot output optim only show models, maybe iterations as smaller points, or plot both
- [ ] VTK plotter
- [ ] Remove rcParams and explicitely set all plot attributes in calls
- [x] Windows extend the height of the waveform plot

#### Pyaflowa:
- [X] fix windows within a single iteration
- [ ] source receiver json file with event-station information such as lat/lon, dist and BAz, only once for m00s00
- [ ] ensure end_pad >= PAR.NT * PAR.DT so that observations will be as long as synthetics
- [ ] sanity check parameters in initialize
- [ ] if given parameters (i.e. misspelled) are not used, notify the user somehow, maybe in seisflows
- [ ] auxiliary_data.Statistics: store time shift (min, max, time shift components), format misfit with 'e' notation not 'f'
- [X] auxiliary_data store windows and adjoint sources by step lengths? Or just store s00 because that is the "initial" step
      Waveforms, misfit windows and adjoint sources are now stored by model and step count in Manager
- [ ] allow reading data from previously collected hdf5 files from other inversions, to avoid recollecting obs data. would need to sanitize the datasets so that none of the other data gets in the way
- [X] make rcvs.vtk for all receivers
- [X] make srcs_w_depth.vtk, and perhaps with depth as planar slices
      New src_vtk_from_specfem() function that takes constant x y or z values to make planar slices

