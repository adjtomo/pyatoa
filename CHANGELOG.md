# CHANGELOG

## v0.2.1

- Updated internal call structures to deal with Pyadjoint v0.2.1 API changes
- Changed internal test ASDFDataSet and created a script to generate new dataset
  because the old one had no way of being remade.
- New Docs + Example + Example data: Processing data with Pyatoa and MPI
- Remove GitHub Pip install links for PySEP, Pyflex and Pyadjoint

## v0.2.2

- Bugfix: Gatherer attempting to access a removed Config parameter
- Resolve PyPDF2 -> PyPDF dependency deprecation warning
- Bugfix: Manager.standardize() only resamples if required, otherwise small time shifting is introduced (Issue \#34)
