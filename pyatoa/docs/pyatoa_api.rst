Pyatoa API
==========

Pyatoa’s main functionalities are contained in its ``core`` module.
Utility functions that are used throughout the package are contained in
the ``utils`` module. Visualization routines are found in the
``visuals`` module.

--------------

Core
----
The main classes of Pyatoa which provide core functionality.

.. toctree::
    :maxdepth: 1

    modules/core.config
    modules/core.manager
    modules/core.gatherer
    modules/core.inspector
    modules/core.pyaflowa
   

--------------

Utils
-----
Utility functions that support the core classes.

.. toctree::
    :maxdepth: 1

    modules/utils.calculate
    modules/utils.form
    modules/utils.images
    modules/utils.process
    modules/utils.read
    modules/utils.srcrcv
    modules/utils.window
    modules/utils.write


ASDF Utils
***********
Utilities designed specifically to facilitate manipulation of ASDFDataSets.

.. toctree::
    :maxdepth: 1

    modules/utils.asdf.add
    modules/utils.asdf.clean
    modules/utils.asdf.load
    modules/utils.asdf.write


--------------

Visuals
-------
Classes and associated functions for generating figures.

.. toctree::
    :maxdepth: 1

    modules/visuals.wave_maker
    modules/visuals.map_maker
    modules/visuals.insp_plot
    modules/visuals.comp_wave
    modules/visuals.improve_wave


--------------

Plugins
-------
A catch-all directory for functions, classes and scripts that don't fit above.

