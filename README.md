# pyatoa
## Python's Adjoint Tomography Operational Assistance
### a misfit quantification package for the modern tomographer.

Pyatoa provides a level of abstraction over a few Python based tools that are useful in the adjoint tomography problem:

Obspy: for seismic data fetching, handling, processing and organization
PyASDF: an HDF5 wrapper for heirarchical data storage which removes the need for complex directory structures and naming schema
Pyflex: a Python port of Flexwin, an automatic time window selection algorithm
Pyadjoint: a package for calculating misfit and creating adjoint sources
Matplotlib: a 2D plotting library used to generate quick plots of waveforms, statistics, misfit etc. for quick data digestion
Basemap: A mapping library to conveniently show source receiver distributions, raypaths, etc. (deprecated, I hope to switch to Cartopy in the future)

Pyatoa provides some pre-fabricated classes and functions which take advantage of the underlying tools, allowing for a streamlined
misfit quantification workflow, that can be used bog-standard, or overwritten with some simple scripting.

The design philosophy of Pyatoa follows closely to Obspy, which provides the foundation for this package: that is, flexible classes and functions that can be used for rapid development, while also providing easy to use functionalities for quick data handling and visualization, but leave room for an easy transition to further, more detailed and advanced scientific research. For example, Pyatoa leaves all data in the form of Obspy streams, or easy to use Python classes (dictionaries, numpy arrays), which makes it dead simple to take the outputs of an inversion and push forward with more rigorous science.
