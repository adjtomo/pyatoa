## For Version 0.2.0

#### Features
- [X] Use Pandas in the Inspector class to do the large-scale data analysis required for all the misfit windows, etc.

#### Bugs
- [ ] Pyflex value Error is being thrown (pyflex.window_selector() line 427 np.abs(noise).max() zero size array)
- [X] FDSNException in gatherer.gather_observed( Uknown Error (timeout))  
      *Added catch for FDSNException in gather_observed() and return st_obs=None*
- [X] Pyaflowa exit gracefully if no data is gathered for an entire event, currently finalize is throwing uncaught errors  
      *Pyaflowa counts successful processes, if none, finalize is skipped*
- [X] Fix: waveform plots not deleted if theyre not included in the composite pdf  
      *Moved purge outside loop and set it to delete all files inside the dir with the correct tag*
- [X] Pyaflowa waveform composites aren't made if a thrifty inversion is done because it skips s00  
      *removed the step count requirement in the if statement*
      
#### Questions
- [ ] Fixed windowing might encounter some problems because the synthetic trace is changing, so the values of max_cc_, cc_shift and dlnA are not being re-evaluated. Can this be remedied? Can we add some functionality to Pyflex to reevaluate misfit values based on waveforms?
- [ ] Is weighting adjoint sources by station proximity something that Pyatoa should do, how could it be implemented?
- [ ] Pyflex: Very close source-receiver distances means P-wave arrival is within the first wavelength, meaning no noise amplitude calculations can take place, and windows are not picked even for good waveforms. Can this be reconciled in Pyflex, or do we need to exclude distances <100km e.g.?
- [X] Should we be able to export ASDFDataSets to hardcoded directory structures. This would provide a form of 'backwards compatability' to old styles of tomography, but if the User already can access a Dataset, do they need this functionality?


#### General
- [ ] Big change: Rename 'model' to 'iteration', and start counting iterations from 1 not 0, to match Seisflows.  
      Refer to 'iterations' as 'evaluations', in reference to function evaluations.  
      Use integers to refer to iterations and step counts, only format for print statements, removes a lot of difficult   
      Current implementation is wrong and misleading, as e.g. m00s03 refers to the line search for model 1.
- [ ] Flag/cache system to tells the Manage/Gatherer/Pyaflowa that an FDSNNoDataException has been thrown, so it won't query FDSN in future iterations
- [ ] Change mapping to Cartopy or just drop mapping capabilities completely? Seems to make things more difficult for not much gain. Or have a completely indepenent module to 
     make maps all in one go given CMTSOLUTIONS and STATION file list. Or use ObsPy mapping functionality to replace Manager.srcrcvmap()? Current Basemap calling is a bit 
     rough.
- [ ] Make network calls optional in e.g. Inspector(), search for station name only and maybe if multiple stations have the same name, then require network code. Station should be enough to identify.
- [ ] Documentation, tutorials, examples and API
- [ ] Tests with PyTest
- [X] ipynb in gitignore
- [X] remove large data files from test data
- [X] make model number, step count formatters standard package wide functions  
     *utils.format now has model() and step() functions that standardize these string formatters*
- [X] check unused kwargs in the config in the case of typos  
      *Added a check call at the end of check() that puts up a user warning if unnused kwargs fall through*
- [X] include __repr__ for all classes
     *config, manager, pyaflowa... todo: inspector  

#### Config
- [ ] Include a set() function to prevent incorrect parameter sets or overwriting read-only parameters
- [ ] ~~Include UTM projection into config and propogate into scripts~~
- [X] Change 'model_number' to model

#### Manager
- [ ] Remove window_amplitude_ratio() from Manager, Pyflex already has this with 'check_global_data_quality'
- [ ] Check if convolve_stf properly performs the time shift  
- [ ] Enforce zero padding front and back of waveforms for short source-receiver distance, as windowing ignores these because the stalta starts too early.
- [ ] Log statement stating fid when saving figures
- [ ] Manager.load() should allow loading misfit windows and adjoint sources as well, will need a warning statement to say waveforms are not processed
- [ ] Functions return self to allow chaining
- [ ] Custom error messages, or specific error messages from Manager rather than returning 0 or 1
- [X] Initialize an empty Manager with an empty Config to remove the need to call Config separately
- [X] Move window by amplitude in Manager.window() into its own function  
     *moved into pyatoa.utils.window and import by Manager*
- [ ] ~~Option to save processed streams in the dataset~~  
     *decided to just allow saving preprocessing attributes to config in pyaflowa, don't think this is super important*  
- [X] Change preprocess() function to take a Manager object, that way when the `overwrite` parameter is called, the User knows exactly
      how to write their new preprocess function  
      *standard preproc() fuction now takes manager as an input, which shortens the call in Manager class*
- [X] Ability to turn off saving for pyasdf dataset  
     *added a 'save_to_ds' bool parameter to the config object that toggles saving for waveforms, metadata, windows, adj_srcs
- [X] Manager.load() take model number from Config if available


#### ASDF
- [ ] Save Pyflex/Pyadjoint Config parameters, not just the map name, incase map names change
- [ ] Generate waveform plots, maps, from an ASDF dataset. As in remove the need to create a Manager just to make 
      waveform plots, if a dataset already exists
- [X] Processing provenance saved into auxiliary_data?  
      *Saved processing stats from obspy stream into the Config object for each model/step*
- [X] Retain step count information for MisfitWindows and AdjointSources

#### Misc.
- [X] Move all hardcoded stuff into plugins, e.g. fault_coordinates, geonet_moment_tensor
- [X] Move tools and visualize out of utils dir into main dir


#### Visuals
- [X] Depth cross section of a Catalog object  
      *In pyaflowa artist *  
- [ ] ~~Plot output optim only show models, maybe iterations as smaller points, or plot both~~  
     *this is okay, Inspector can make more detailed plots if we want, just want a quick visualization**  
- [X] VTK plotter  
     *VTK plotter using mayavi has been started, but needs fine tuning, and implementation into the greater workflow*  
- [ ] ~~Remove rcParams and explicitely set all plot attributes in calls~~   
     *let sleeping dogs lie for now, hopefully this is okay but maybe change in future if it causes problems*  
- [x] Windows extend the height of the waveform plot
- [ ] Plot c_0 * stalta_waterlevel for water level rejection

#### Pyaflowa:
- [X] fix windows within a single iteration  
- [ ] ~~source receiver json file with event-station information such as lat/lon, dist and BAz, only once for m00s00~~  
     *This is taken care of by the Inspector class, don't need to perform twice*  
- [X] sanity check parameters in initialize, e.g. ensure end_pad >= PAR.NT * PAR.DT so that observations will be as long as synthetics  
     *Pyaflowa now has a _check() function that checks the parameters passed in from Seisflows*  
- [ ] if given parameters (i.e. misspelled) are not used, notify the user somehow, maybe in seisflows  
- [ ] ~~auxiliary_data.Statistics: store time shift (min, max, time shift components), format misfit with 'e' notation not 'f'~~   
     *Got rid of Statistics as I never really used it, would prefer to use the Inspector class as that has a cleaner approach to  collecting stats data*  
- [X] auxiliary_data store windows and adjoint sources by step lengths? Or just store s00 because that is the "initial" step  
      *Waveforms, misfit windows and adjoint sources are now stored by model and step count in Manager*  
- [ ] ~~allow reading data from previously collected hdf5 files from other inversions, to avoid recollecting obs data. would need to clean up the datasets so that none of the other data gets in the way~~   
     *this doesn't need to be done in Pyaflowa, can probably be turned into a script or a function in Pyatoa?*  
- [X] make rcvs.vtk for all receivers
- [X] make srcs_w_depth.vtk, and perhaps with depth as planar slices  
      *New src_vtk_from_specfem() function that takes constant x y or z values to make planar slices*
- [X] Waveform improvement script should take advantage of Manager.load() rather than accessing waveform tags
- [X] Inspector convergence plot should plot line search models as non-connected points as in Krischer et al. 2018(?)

#### Inspector:
- [ ] Automatically create a list of windows corresponding to largest time shift, or misfit or dlnA

