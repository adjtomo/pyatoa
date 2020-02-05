## For Version 0.2.0

#### General
- [ ] Proper docs for all functions. methods, classes 
- [ ] Finish writing tests with PyTest
- [ ] ipynb in gitignore
- [ ] remove large data files from test data

#### Manager
- [ ] Initialize an empty Manager with an empty Config to remove the need to call Config separately
- [ ] Move window by amplitude in Manager.window() into its own function
- [ ] Check if convolve_stf properly performs the time shift
- [ ] Option to save processed streams in the dataset

#### ASDF
- [ ] ASDF export all data to directory structure, saving waveforms, inventories, aux data etc into individual 
- [ ] Generate waveform plots, maps, from an ASDF dataset. As in remove the need to create a Manager just to make 
      waveform plots, if a datset already exists

#### Misc.
- [ ] Move all hardcoded stuff into plugins, e.g. fault_coordinates, geonet_moment_tensor
- [ ] Move tools and visualize out of utils dir into main dir

#### Visuals
- [ ] Depth cross section of a Catalog object
- [ ] Plot output optim only show models, maybe iterations as smaller points, or plot both
