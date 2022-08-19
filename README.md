## Python's Adjoint Tomography Operations Assitant  

[![Join the chat at https://gitter.im/pyatoa/community](https://badges.gitter.im/pyatoa/community.svg)](https://gitter.im/pyatoa/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Documentation Status](https://readthedocs.org/projects/pyatoa/badge/?version=latest)](https://pyatoa.readthedocs.io/en/latest/?badge=latest)

### Misfit assessment for the modern tomographer

---
`Documentation` can be found on Read the Docs: https://pyatoa.readthedocs.io/en/devel

---
<p align="center">
  <img src="docs/images/pyatoa_logo_w_text.png" />
</p>

**Pyatoa**\* is a waveform-based misfit quantification package, with additional tools for measurement aggregation, and visualizations of inversion results. It is meant to facilitate the assessment of seismic inversions. Under the hood, **Pyatoa** is built on, and provides a high-level API for, a few key Python packages:

**[ObsPy:](https://github.com/obspy/obspy/wiki)** for seismic data fetching, handling, processing and organization    
**[Pyflex:](https://krischer.github.io/pyflex/)** a Python port of Flexwin, for automatic time window selection  
**[Pyadjoint:](http://krischer.github.io/pyadjoint/)** evaluating misfit functions and creating adjoint sources  
**[PyASDF:](https://seismicdata.github.io/pyasdf/)** heirarchical, self-describing storage of seismic data  
**[Pandas:](https://pandas.pydata.org/)** large-scale aggregation and manipulation of measurement information

**Pyatoa** can be used in scripting, interactive Python environments, or written into larger workflow tools. Although applicable in a standalone maner, Pyatoa was designed as a tool to be used in conjunction with [SeisFlows3](https://github.com/adjtomo/seisflows3), an automated workflow for seismic inversions, and [SPECFEM3D Cartesian](https://geodynamics.org/cig/software/specfem3d/), a numerical solver for seismic wave propogation.

The design philosophy of **Pyatoa** is easy-to-use custom-built objects that make tomography research flexible, rapid, and repeatable.

<sub> \**pronounced Py-uh-toe-uh (ˈpaɪəˈtoʊə), inspired by the famed volcano Krakatoa* </sub>

