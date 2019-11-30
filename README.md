# pyatoa  
### Python's Adjoint Tomography Operational Assistance
### A misfit quantification package for the modern tomographer

Pyatoa provides a level of abstraction over a few Python based tools that are useful in the adjoint tomography problem:

**Obspy:** for seismic data fetching, handling, processing and organization  
**Pyflex:** a Python port of Flexwin, an automatic time window selection algorithm  
**Pyadjoint:** a package for calculating misfit and creating adjoint sources  
**PyASDF:** heirarchical data storage that removes the need for complex directory structures and naming schema  
**Matplotlib:** plotting library for quick plots of waveforms, statistics, misfit etc.  
**Basemap:** A mapping library for source receiver distributions, raypaths, etc. (deprecated, Cartopy in the future)  

Pyatoa provides pre-fabricated classes and functions which take advantage of these underlying tools, it allows for streamlined
misfit quantification workflows that can be used bog-standard, or overwritten to a researchers needs with simple Python scripting.

The design philosophy of Pyatoa follows closely with Obspy; that is, flexible classes and functions that can be used for rapid development, quick data handling and simple but powerful visualization, while still leaving room for easy transition to more detailed scientific research. 

Pyatoa has been developed in conjunction with Seisflows, an automated workflow for seismic imaging. Included are plugin scripts that allow Pyatoa to be slotted in neatly to a Seisflows workflow template.

Detailed documentation coming soon!
