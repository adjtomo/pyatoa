Overview
==============

What is Pyatoa?
~~~~~~~~~~~~~~~
At its core, Pyatoa is a waveform assessment package. Put very simply, it is 
designed to compare one set of wiggles to another set of similar wiggles.

The aim of Pyatoa is to simplify misfit quantification in adjoint tomography. 
It does this by providing a high-level API to reduce the amount of code required 
to accomplish repeatable tasks involved in waveform comparisons.

As the behavior of full waveform inversion is dependent on the input data, a 
high level of care must be taken in curating what is fed in. 
Pyatoa provides this attention to detail by carefully addressing the following
tasks: seismic data gathering, waveform preprocessing, 
time windowing, misfit quantification, and adjoint source creation. 

To ensure an inversion stays on the rails, Pyatoa includes a custom data
structure for data access and storage, as well as internal fault 
tolerance and sanity checks throughout. Detailed logging is meant to ensure that
Pyatoa is not a black box.

Outside the standard routines, Pyatoa also includes a measurement aggregation 
tool to simplify bulk measurement assessement, while a series of plotting 
routines facilitate the rapid visualization of standard inversion results.

Pyatoa is open-source and completely Python based.


What isn't Pyatoa?
~~~~~~~~~~~~~~~~~~

Pyatoa is not a standalone adjoint tomography tool, it does not have the 
capability to generate sythetic waveforms, submit jobs on HPC systems,  
interface with numerical models, etc. Rather, it is meant to augment the 
capabilities of external numerical solvers (e.g. SPECFEM3D) and workflow tools
(e.g. SeisFlows).

Pyatoa is not intellegent. Although it provides sanity checks and fault 
tolerance, it cannot do your science for you. Careful attention must be paid in 
properly choosing input data and processing parameters based on the problem at 
hand. Garbage in will always mean garbage out.

Pyatoa is not quiet, it logs almost every task that it performs. Although this
can be overwhelming, each log statement is important. If you are running Pyatoa
for the first time, be sure to read the logs to get an idea of what is going
on under the hood.


How do I use Pyatoa?
~~~~~~~~~~~~~~~~~~~~

Pyatoa was written following the design philosophy of ObsPy, that means it's 
meant to be used as a Python tool, NOT as a standalone command-line tool, or 
GUI based program.
Pyatoa can be invoked through scripting, or in interactive Python
environments such as the Python interpretor, IPython, Jupyter Notebooks, etc.


The following introductory notebooks should provide a quick overview of how 
Pyatoa and its underlying functionalities should be used. They are meant 
to provide a jumping off point to communicate how one might leverage Pyatoa 
during a tomographic inversion, however they are by no means all encompassing.
The later tutorials and API pages provide more detailed views of the package.

.. toctree::
   :maxdepth: 3
   
   quickstart
   prepwork
    
