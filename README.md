#pyatoa
Adjoint Tomography Operational Assistance for Python
Pyatoa is a workflow management package for the adjoint tomography problem.
The idea behind the pacakge is that Pyatoa can take care of all steps in the adjoint problem from data retrieval to adjoint source creation.
Pyatoa builds an easy to use skeleton around a few key python packages: PyFlex for misfit window identification, PyAdjoint for misfit quantification and adjoint source creation and PyASDF for input/output.
This package is built around SPECFEM3D Cartesian, and expects, for example, synthetics and file naming schemes standard in the SPECFEM3D workflow.