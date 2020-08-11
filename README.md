## Python's Adjoint Tomography Operations Assitant  
### Misfit assessment for the modern tomographer

---
Documentation can be found here (in development): https://pyatoa.rtfd.io

---
<p align="center">
  <img src="pyatoa/docs/pyatoa_logo_w_text.png" />
</p>

Pyatoa is a waveform-based misfit quantification package, with additional tools for measurement aggregation, and visualizations of inversion results. It is meant to facilitate the assessment of seismic inversions. Under the hood, Pyatoa is built on a few key Python packages:

**[ObsPy:](https://github.com/obspy/obspy/wiki)** for seismic data fetching, handling, processing and organization    
**[Pyflex:](https://krischer.github.io/pyflex/)** a Python port of Flexwin, for automatic time window selection  
**[Pyadjoint:](http://krischer.github.io/pyadjoint/)** evaluating misfit functions and creating adjoint sources  
**[PyASDF:](https://seismicdata.github.io/pyasdf/)** heirarchical, self-describing storage of seismic data  
**[Pandas:](https://pandas.pydata.org/)** large-scale aggregation and manipulation of measurement information

Pyatoa is meant to be used in scripting, interactive environments, or written into larger workflow tools. It is used to facilitate the gathering, processing, windowing, measurement, visualization, storage and assessment of waveform data in a tomographic inversion. Although applicable in a standalone maner, Pyatoa was designed as a tool to be used in conjunction with [Seisflows](https://github.com/rmodrak/seisflows), an automated workflow for seismic inversions, and [Specfem3D Cartesian](https://geodynamics.org/cig/software/specfem3d/), a numerical solver for seismic wave propogation.

The design philosophy of Pyatoa follows closely with Obspy; that is, flexible classes and functions that can be used for rapid development, repeatable measurements, quick data handling, and simple but powerful visualization, with enough flexibility for detailed scientific research. 

