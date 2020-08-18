## Python's Adjoint Tomography Operations Assitant  
### Misfit assessment for the modern tomographer

---
`Documentation` can be found on Read the Docs: https://pyatoa.rtfd.io (in development)

---
<p align="center">
  <img src="pyatoa/docs/pyatoa_logo_w_text.png" />
</p>

`Pyatoa` is a waveform-based misfit quantification package, with additional tools for measurement aggregation, and visualizations of inversion results. It is meant to facilitate the assessment of seismic inversions. Under the hood, Pyatoa is built on a few key Python packages:

**[ObsPy:](https://github.com/obspy/obspy/wiki)** for seismic data fetching, handling, processing and organization    
**[Pyflex:](https://krischer.github.io/pyflex/)** a Python port of Flexwin, for automatic time window selection  
**[Pyadjoint:](http://krischer.github.io/pyadjoint/)** evaluating misfit functions and creating adjoint sources  
**[PyASDF:](https://seismicdata.github.io/pyasdf/)** heirarchical, self-describing storage of seismic data  
**[Pandas:](https://pandas.pydata.org/)** large-scale aggregation and manipulation of measurement information

`Pyatoa` can be used in scripting, interactive Python environments, or written into larger workflow tools. Although applicable in a standalone maner, Pyatoa was designed as a tool to be used in conjunction with [SeisFlows](https://github.com/rmodrak/seisflows), an automated workflow for seismic inversions, and [SPECFEM3D Cartesian](https://geodynamics.org/cig/software/specfem3d/), a numerical solver for seismic wave propogation.

The design philosophy of `Pyatoa` follows closely with Obspy, i.e. custom-built objects that make tomography research flexible, rapid, and repeatable.

