Overview
==============

What is Pyatoa?
~~~~~~~~~~~~~~~
At its core, Pyatoa is a waveform assessment package. The short version: it's
designed to compare one set of wiggles to another set of (similar) wiggles. The
objects within this package are meant to facilitate, augment, or complement
this core functionality. These include:

- Metadata and waveform collection
- Waveform standardization and preprocesing
- Time windowing and adjoint source generation
- Hierarchical data storage of waveforms, metadata, and measurements
- Bulk measurement aggregation and analysis
- Waveform and measurement plotting
- Interface with broader workflow tools

The aim of Pyatoa is to simplify and automate misfit quantification in adjoint
tomography. It does so by providing a high-level API to reduce the amount of
code required to accomplish repeatable tasks involved in waveform comparisons.

As the behavior of full waveform inversion is dependent on the input data, a 
high level of care must be taken in curating what is fed in. To ensure an
inversion stays on the rails, Pyatoa includes a custom data
structure for data access and storage, as well as internal fault 
tolerance and sanity checks throughout. Detailed logging is meant to ensure that
Pyatoa is not a black box.

Outside the standard routines, Pyatoa also includes a measurement aggregation 
tool to simplify bulk measurement assessment, while a series of plotting
routines facilitate visualization of standard inversion results.

Pyatoa is open-source and completely Python based.



What isn't Pyatoa?
~~~~~~~~~~~~~~~~~~

Pyatoa is not a standalone adjoint tomography workflow tool, it does not have
the capability to generate sythetic waveforms, submit jobs on HPC systems,
interface with numerical models, etc. Rather, it was built to augment the
capabilities of external numerical solvers (e.g. SPECFEM3D) and workflow tools
(e.g. SeisFlows3).

Pyatoa is not smart. Although it provides quality checks along the way, it
cannot do your science for you. Careful attention must be paid in
properly choosing input data and processing parameters based on the problem at 
hand. Garbage in == garbage out.

Pyatoa is not quiet, it logs almost every task that it performs. Although this
can be overwhelming, each log statement hold some importance. If you are running
Pyatoa for the first time, be sure to read the logs to get an idea of what is
going on under the hood.


Why is Pyatoa (necessary)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some question that you might have regarding Pyatoa are: 

- Why is Pyatoa necessary? 
- If workflow tools like SeisFlows3 exist, what role does Pyatoa play?
- Is Pyatoa necessary to run seismic inversions? 

Well, SeisFlows3 was originally written (as SeisFlows) as an automated workflow 
tool for full waveform inversion, allowing for generalized interfacing with a 
number of numerical solvers and compute interfaces. 

However, within the original SeisFlows, we identified some key features that are
missing from the package but necessary for earthquake-based tomography, namely: 
data gathering, waveform preprocessing, windowing, 
flexible adjoint source creation, inversion assessment, and figure generation.

Now, these tasks can of course be performed manually with existing tools, e.g.,
one can preprocess with ObsPy, window with Flexwin, create adjoint sources 
with NumPy, assess an inversion with Pandas and generate figures with 
Matplotlib, Matlab, GMT etc., however the name of the game here is automation
and reproducibility. 

If everyone has their own individual codes to perform these tasks, then each 
researcher must effectively re-invent the wheel. Consequently, Pyatoa was 
developed to provide a high-level interface for users (and by extension 
SeisFlows3) to automate the above-named tasks as well as provide a platform for 
the tomography community to improve upon collectively.

If a user does not need the above capabilities, e.g., while running very simple
2D synthetic inversions where data-synthetic misfits can be computed along the
entire waveform, then Pyatoa is not necessary. If a user only needs to make 
data-synthetic comparisons on a set of similar waveforms, then only Pyatoa is 
needed. And if a user wants to automate an entire seismic inversion involving
real data, then we recommend using a combination of Pyatoa and SeisFlows3 (see
:doc:`Pyaflowa </pyaflowa>` documentation).


How do (I use) Pyatoa?
~~~~~~~~~~~~~~~~~~~~~~~

Pyatoa was written following the design philosophy of ObsPy, that means it's 
meant to be used as a Python tool, NOT as a standalone command-line tool, or 
GUI based program.
Pyatoa can be invoked through scripting, or in interactive Python
environments such as the Python interpreter, IPython, Jupyter Notebooks, etc.

The notebooks found in the introduction section should provide a quick overview
of how Pyatoa and its underlying functionalities should be used. They are meant
to provide a jumping off point to communicate how one might leverage Pyatoa 
during a tomographic inversion, however they are by no means all encompassing.
The later tutorials and API pages provide more detailed views of the package.

    
