Utilities
=========

``Pyatoa`` comes with utilities that facilitate interaction with a few
core packages.

--------------

ASDFDataSets
------------

By default PyASDF is a powerful package that natively contains almost
all the functions we need to store seismological data. However since
time windows and adjoint sources are non-standard seismological objects,
some additional utility functions are required to read and write them
to/from an ASDFDataSet.

--------------

SPECFEM I/O
-----------

``SPECFEM`` reads and writes its own native file systems. ``Pyatoa``
contains custom functions to read in SPECFEM file formats as ObsPy
objects, and conversely write ObsPy objects into SPECFEM-ready files.

Catalogs <=> CMTSOLUTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~

``SPECFEM`` requires earthquake sources to be input as CMTSOLUTION
files,


--------------

Seisflows
---------

--------------

Plotting models from .VTK files
-------------------------------

``SPECFEM`` outputs `Visualization ToolKit (.vtk) <https://vtk.org/>`__
files for visualizing meshes, models, gradients etc. Although standalone
programs can be used to visualize these models, Pyatoa uses the Python
package `Mayavi <https://docs.enthought.com/mayavi/mayavi/>`__ to
produce standard visualizations of Cartesian meshes. This facilitates
rapid assessment of models, but higher-quality model visualizations are
deferred to other packages like
`ParaView <https://www.paraview.org/>`__.

.. code:: ipython3

    # Allow mayavi visualizations to appear in notebook
    from mayavi import mlab
    mlab.init_notebook()
    
    from pyatoa.visuals.vtk_modeler import VTKModeler
    
    vm = VTKModeler()
    vm.load(fid="../tests/test_data/big_test_data/test_vs_model.vtk")
    vm.depth_slice(depth_km=5, show=True)


::


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    /tmp/ipykernel_9648/1554201936.py in <module>
          1 # Allow mayavi visualizations to appear in notebook
    ----> 2 from mayavi import mlab
          3 mlab.init_notebook()
          4 
          5 from pyatoa.visuals.vtk_modeler import VTKModeler


    ModuleNotFoundError: No module named 'mayavi'


--------------

Plotting waveform improvement using ASDFDataSets
------------------------------------------------

