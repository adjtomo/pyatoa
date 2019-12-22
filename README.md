## Python's Adjoint Tomography Operations Assitant

---
Docs (in development): https://pyatoa.readthedocs.io/en/latest/

---

![Logo](pyatoa/docs/pyatoa_logo.png)
### A misfit quantification package for the modern tomographer

Pyatoa is a waveform-based misfit quantification package, which provides abstraction over a few Python based tools.

**[Obspy:](https://github.com/obspy/obspy/wiki)** for seismic data fetching, handling, processing and organization.    
**[Pyflex:](https://krischer.github.io/pyflex/)** a Python port of Flexwin, an automatic time window selection algorithm.  
**[Pyadjoint:](http://krischer.github.io/pyadjoint/)** a package for calculating misfit and creating adjoint sources.  
**[PyASDF:](https://seismicdata.github.io/pyasdf/)** heirarchical data storage for seismic data.  
**[Matplotlib:](https://matplotlib.org/)** 2D plotting library for visualization of waveforms, statistics, misfit etc.  
**[Basemap:](https://matplotlib.org/basemap/)** A mapping library for source receiver distributions, raypaths, etc. (deprecated, Cartopy in the future).  

Pyatoa can be used bog-standard, or customized with simple Python scripting, to gather, process, window, measure misfit, visualize, save and export waveform data. It was designed as a tool to be used in conjunction with [Seisflows](https://github.com/rmodrak/seisflows), an automated workflow for seismic inversions, and [Specfem3D Cartesian](https://geodynamics.org/cig/software/specfem3d/), a numerical solver for seismic wave propogation.

The design philosophy of Pyatoa follows closely with Obspy; that is, flexible classes and functions that can be used for rapid development, quick data handling and simple but powerful visualization, while still leaving room for easy transition to more detailed scientific research. 

