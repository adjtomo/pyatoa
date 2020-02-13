## For Version 0.2.0

#### Bugs
- [ ] FDSNException in gatherer.gather_observed( Uknown Error (timeout))
- [ ] Pyflex value Error is being thrown (pyflex.window_selector() line 427 np.abs(noise).max() zero size array)
- [ ] Pyaflowa exit gracefully if no data is gathered for an entire event, currently finalize is throwing uncaught errors
- [ ] Fix: waveform plots not deleted if theyre not included in the composite pdf

#### General
- [ ] Proper docs for all functions. methods, classes 
- [ ] Finish writing tests with PyTest
- [ ] ipynb in gitignore
- [ ] remove large data files from test data
- [ ] make model number, step count formatters standard package wide functions

#### Manager
- [ ] write `__call__` attribute
- [ ] Initialize an empty Manager with an empty Config to remove the need to call Config separately
- [ ] Move window by amplitude in Manager.window() into its own function
- [ ] Check if convolve_stf properly performs the time shift
- [ ] Option to save processed streams in the dataset
- [ ] Change preprocess() function to take a Manager object, that way when the `overwrite` parameter is called, the User knows exactly
      how to write their new preprocess function
- [ ] Fixed windowing might encounter some problems because the synthetic trace is changing, so the values of max_cc_, cc_shift and dlnA       are not being re-evaluated. Can this be remedied?
- [ ] Weighting adjoint sources by station proximity?

#### ASDF
- [ ] ASDF export all data to directory structure, saving waveforms, inventories, aux data etc into individual 
- [ ] Generate waveform plots, maps, from an ASDF dataset. As in remove the need to create a Manager just to make 
      waveform plots, if a datset already exists
- [ ] Processing provenance saved into auxiliary_data?
- [ ] Replace station lat lon getters with ASDFDataSet.get_all_coordinates()

#### Misc.
- [ ] Move all hardcoded stuff into plugins, e.g. fault_coordinates, geonet_moment_tensor
- [ ] Move tools and visualize out of utils dir into main dir

#### Visuals
- [ ] Depth cross section of a Catalog object
- [ ] Plot output optim only show models, maybe iterations as smaller points, or plot both
- [ ] VTK plotter
- [ ] Remove rcParams and explicitely set all plot attributes in calls

#### Pyaflowa:
- [ ] fix windows within a single iteration
- [ ] source receiver json file with event-station information such as lat/lon, dist and BAz, only once for m00s00
- [ ] ensure end_pad >= PAR.NT * PAR.DT so that observations will be as long as synthetics
- [ ] sanity check parameters in initialize
- [ ] if given parameters (i.e. misspelled) are not used, notify the user somehow, maybe in seisflows
- [ ] auxiliary_data.Statistics: store time shift (min, max, time shift components), format misfit with 'e' notation not 'f'
- [ ] auxiliary_data store windows and adjoint sources by step lengths? Or just store s00 because that is the "initial" step

