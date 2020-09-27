## For Version 0.2.0

#### Features
- [ ] Introduce multiprocessing with concurrent.futures or multiprocessing packages.
- [X] Use Pandas in the Inspector class to do the large-scale data analysis required for all the misfit windows, etc.

#### Bugs
- [X] Pyflex value Error is being thrown (pyflex.window_selector() line 427 np.abs(noise).max() zero size array). Very close source-receiver distances means P-wave arrival is within the first wavelength, meaning no noise amplitude calculations can take place, and windows are not picked even for good waveforms. Can this be reconciled in Pyflex, or do we need to exclude distances <100km e.g.?
*Pyatoa now catches these as known exceptions*
- [ ] If Manager.gather() is called after Manager.load() the gatherer attribute has no origintime. Also the synthetic tag is incorrectly set. These need to be properly propogated in the load() command.
- [ ] Gathered obs data can sometimes have multiple (>1) traces per component. Need to cap this to 1 per component and also let the user know that this is happening.
- [ ] If windows are fixed, Pyflex no longer controls the maximum criteria within the windows, so e.g. time shift is allowed to exceed the maximum allowed time shift.
- [X] FDSNException in gatherer.gather_observed( Uknown Error (timeout))  
      *Added catch for FDSNException in gather_observed() and return st_obs=None*
- [X] Pyaflowa exit gracefully if no data is gathered for an entire event, currently finalize is throwing uncaught errors  
      *Pyaflowa counts successful processes, if none, finalize is skipped*
- [X] Fix: waveform plots not deleted if theyre not included in the composite pdf  
      *Moved purge outside loop and set it to delete all files inside the dir with the correct tag*
- [X] Pyaflowa waveform composites aren't made if a thrifty inversion is done because it skips s00  
      *removed the step count requirement in the if statement*
      
#### Questions
- [ ] Would multithreading calls to obspy.clients.fdsn.Client be useful? Since these calls are waiting on a webservice routine they may benefit from asynchronous thread calls, but would the speed-up be at all useful? Maybe if done en-masse for all stations and events at the same time with a single giant gather function. Perhaps a one-time gathering functionality is warranted here? Something like Gatherer.gather_all(event_list, station_list), to gather all EventXML's, all StationXML's, and all observed waveforms.
- [ ] Should we try to find a way to NOT repeatedly save StationXML files, because as of currently, StationXML files are gathered and stored for each event. Not a heavy storage demand, but not very elegant either. This is perhaps a good thing, though, because response information may change temporally, and we currently gather StationXML information based on event origin time, which means there is a change that these files are different. Need to discuss with someone.
- [ ] Is weighting adjoint sources by station proximity something that Pyatoa should do, how could it be implemented?
- [X] Fixed windowing might encounter some problems because the synthetic trace is changing, so the values of max_cc_, cc_shift and dlnA are not being re-evaluated. Can this be remedied? Can we add some functionality to Pyflex to reevaluate misfit values based on waveforms?
*Pyflex already had an option to recalculate window parameters for a given set of windows, this is now implemented in Pyatoa*
- [X] Should we be able to export ASDFDataSets to hardcoded directory structures. This would provide a form of 'backwards compatability' to old styles of tomography, but if the User already can access a Dataset, do they need this functionality?


#### General
- [ ] Standardize channel naming, perhaps the same as LASIF to push for consistency across tools.
- [X] Big change: Rename 'model' to 'iteration', and start counting iterations from 1 not 0, to match Seisflows.  
      Refer to 'iterations' as 'evaluations', in reference to function evaluations.  
      Use integers to refer to iterations and step counts, only format for print statements, removes a lot of difficult   
      Current implementation is wrong and misleading, as e.g. m00s03 refers to the line search for model 1.
- [X] Flag/cache system to tell the Manage/Gatherer/Pyaflowa that an FDSNNoDataException has been thrown, so it won't query FDSN in future iterations  
*Just set client to None in future iterations to stop gathering from FDSN, not as elegant but simpler*
- [ ] Change mapping to Cartopy or just drop mapping capabilities completely? Seems to make things more difficult for not much gain. Or have a completely indepenent module to 
     make maps all in one go given CMTSOLUTIONS and STATION file list. Or use ObsPy mapping functionality to replace Manager.srcrcvmap()? Current Basemap calling is a bit 
     rough.
- [ ] Make network calls optional in e.g. Inspector(), search for station name only and maybe if multiple stations have the same name, then require network code. Station should be enough to identify.
- [ ] Allow additional windowing criteria that avoids choosing windows for stations that are too close/too far from one another.
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
- [ ] Put location (LOC) wildcard and channel (CHA) wildcards in Config and make those accessible to the Manager when gathering data. That way these don't have to be hardcoded in the workflow scripts.
- [ ] Include a set() function to prevent incorrect parameter sets or overwriting read-only parameters
- [ ] ~~Include UTM projection into config and propogate into scripts~~
- [X] Change 'model_number' to model

#### Manager
- [ ] Manager shouldn't throw general ManagerError but actual explicit exceptions related to each part of the processing scheme?
- [ ] Calculate full waveform difference and save in ASDFDataSet, for use in variance reduction
- [ ] Remove window_amplitude_ratio() from Manager, Pyflex already has this with 'check_global_data_quality'
- [ ] Check if convolve_stf properly performs the time shift  
- [ ] Enforce zero padding front and back of waveforms for short source-receiver distance, as windowing ignores these because the stalta starts too early.
- [ ] Log statement stating fid when saving figures
- [ ] Manager.load() should allow loading misfit windows and adjoint sources as well, will need a warning statement to say waveforms are not processed
- [ ] Custom error messages, or specific error messages from Manager rather than returning 0 or 1
- [X] Functions return self to allow chaining
- [ ] ~~Include window weights into adjoint source calculation~~
     *Pyflex only uses User-input window weight functions, no instrinsic weighting is provided*
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
- [X] Save Pyflex/Pyadjoint Config parameters, not just the map name, incase map names change
- [ ] ~~Generate waveform plots, maps, from an ASDF dataset. As in remove the need to create a Manager just to make waveform plots, if a dataset already exists~~
*Not quite, but Manager.load() allows pretty quick access to data in an ASDFDataSet, and provides the machinery to make the plots. Only a few lines of code required*
- [X] Processing provenance saved into auxiliary_data?  
      *Saved processing stats from obspy stream into the Config object for each model/step*
- [X] Retain step count information for MisfitWindows and AdjointSources

#### Misc.
- [ ] read_stations should have built in functions to set the order of stations, alphabetical, by lat, by lon?
- [X] Move all hardcoded stuff into plugins, e.g. fault_coordinates, geonet_moment_tensor
- [X] Move tools and visualize out of utils dir into main dir

#### Visuals
- [ ] Manager plotter use a separate smaller subplot for window information, including sta/lta, rejected windows, etc.
- [ ] Kwargs in docstrings and raise TypeError for undefined kwargs
- [ ] Lower alpha of legend in waveform plotter
- [ ] Show the extent of the window search criteria from Pyflex
- [ ] Plot c_0 * stalta_waterlevel for water level rejection
- [ ] Allow plotting time dependent rejection criteria from Pyflex
- [X] Automatically combine map and waveform plot in matplotlib, rather than after generating .png's
     *Done with GridSpec and a MangerPlotter class, works pretty good!*
- [X] More info in manager plotter title, such as event depth, magnitude, srcrcv distance, baz
     *Not needed now that map and wav plots are combined*
- [X] Depth cross section of a Catalog object  
      *In pyaflowa artist *  
- [ ] ~~Plot output optim only show models, maybe iterations as smaller points, or plot both~~  
     *this is okay, Inspector can make more detailed plots if we want, just want a quick visualization**  
- [X] VTK plotter  
     *VTK plotter using mayavi has been started, but needs fine tuning, and implementation into the greater workflow*  
- [X] Remove rcParams and explicitely set all plot attributes in calls
- [x] Windows extend the height of the waveform plot

#### Pyaflowa (Deprecated):
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
- [X] Inspector should be able to append new data only, rather than having to aggregate from scratch/
- [ ] merge() or += or + function to combine two inspectors which renames iterations of the second inspector incase two different inversion legs are run
- [ ] Automatically create a list of windows corresponding to largest time shift, or misfit or dlnA
- [ ] misfits() should default to misfit per model, not step, needs to be done after the name change of model > iteration 
- [ ] save focal mechanism attributes from sources (if available)

#### Possible Inspector Features:
- [X] Merge() datasets if inversions are run in separate directories, allow one to be appended to the other with change in iteration label.
- [ ] Function to automatically create list of maximum 'key' (e.g. cc_shift_in_seconds') for a given iter/step/event/station etc.  

