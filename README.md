# pyatoa
## Python's Adjoint Tomography Operational Assistance
### A misfit quantification package for the modern tomographer.

Pyatoa provides a level of abstraction over a few Python based tools that are useful in the adjoint tomography problem:

**Obspy:** for seismic data fetching, handling, processing and organization  
**PyASDF:** an HDF5 wrapper for heirarchical data storage which removes the need for complex directory structures and naming schema  
**Pyflex:** a Python port of Flexwin, an automatic time window selection algorithm  
**Pyadjoint:** a package for calculating misfit and creating adjoint sources  
**Matplotlib:** a 2D plotting library used to generate quick plots of waveforms, statistics, misfit etc. for quick data digestion  
**Basemap:** A mapping library to conveniently show source receiver distributions, raypaths, etc. (deprecated, I hope to switch to Cartopy in the future)  

Pyatoa provides pre-fabricated classes and functions which take advantage of these underlying tools, which allow for streamlined
misfit quantification workflows, that can be used bog-standard, or overwritten to the researchers needs with simple scripting.

The design philosophy of Pyatoa follows closely with Obspy; that is, flexible classes and functions that can be used for rapid development, quick data handling and simple but powerful visualization. This leaves room for transitions to more detailed scientific research. 

Detailed Documentation coming soon!
