=========
Tutorials
=========
Tutorial 1: core classes; general workflow
------------------------------------------

| Here we provide a brief introduction on how ``Pyatoa`` is intended to
  be used standalone.
| ``Pyatoa`` is normally very quiet to avoid unnecessary outputs,
  however it really babbles on about everything.
| We can use the ``logging`` module to see what’s going on under the
  hood.

.. code:: ipython3

    %pylab inline
    import obspy
    import pyatoa
    import logging
    
    logger = logging.getLogger("pyatoa")
    logger.setLevel(logging.DEBUG)

pyatoa.Config()
~~~~~~~~~~~~~~~

``Pyatoa`` is centralized around a ``Config`` object which controls all
the parameters deemed necessary in the misfit quantification workflow.
The ``Config`` object supports print statements to tell the User what
parameters are available, and how they are set.

.. code:: ipython3

    cfg = pyatoa.Config()
    print(cfg)

We can see here that we have some default parameters set, e.g. the
filter bands ``min_period`` and ``max_period``, as well as
configurations for the ``pyflex_map`` which specifies the parameters of
the ``pyflex_config``. We also can set the ``adj_src_type`` which
controls the parameters of the ``pyadjoint_config``. For now we will
ignore the rest of the parameters and move onto the ``Manager`` class.

pyatoa.Manager()
~~~~~~~~~~~~~~~~

The ``Manager`` class always requires a ``Config`` input parameter.
Printing the ``Manager`` shows us the ``Pyatoa`` data and workflow
status.

.. code:: ipython3

    mgmt = pyatoa.Manager(cfg)
    print(mgmt)

We can see above that we have no data collected and our workflow is
incomplete. Although ``Pyatoa`` comes with data gathering capabilities,
for this example we will just read in some test data. We will see that
once we set the data in the Manager class, the print statement updates
to show us how many Traces we have in our Streams, and the name of the
station in our inventory object.

.. code:: ipython3

    mgmt.inv = obspy.read_inventory("../tests/data/test_inv.xml")
    mgmt.st_obs = obspy.read("../tests/data/test_obs_data.ascii")
    mgmt.st_syn = obspy.read("../tests/data/test_syn_m00_data.ascii")
    print(f"OBSSERVED DATA\n\tsampling_rate:{mgmt.st_obs[0].stats.sampling_rate}, npts:{mgmt.st_obs[0].stats.npts}\n")
    print(f"SYNTHETIC DATA\n\tsampling_rate:{mgmt.st_syn[0].stats.sampling_rate}, npts:{mgmt.st_syn[0].stats.npts}\n")
    print(mgmt)

We now have the minimum data required to begin the workflow. However, we
cannot create misfit windows or measure misfit yet, because our traces
have different sampling rates, start and end-times, and spectral
content. If we try to run the window() and measure() functions, we will
be met with some check-stops telling us that we must first run other
functions.

.. code:: ipython3

    mgmt.window()
    mgmt.measure()

Data processing
~~~~~~~~~~~~~~~

The first step after gathering data is to standardize Streams.
standardize() will resample and trim the data so that sampling rate,
number of points, and start and end times are the same between the
Observed and Synthetic traces. Data by default should conform to
Synthetics, due to the requirements of a solver to receive inputs in the
same manner as its outputs, but running
``mgmt.standardize(standardize_to='obs')`` allows the User to override
this.

.. code:: ipython3

    mgmt.standardize()

.. code:: ipython3

    mgmt.preprocess()

.. code:: ipython3

    print(mgmt.st)

We can now see that after performing standardization and preprocessing,
3 of the workflow flags have been set ``True``, which means we can move
onto the the window() and measure() functions. Standardization and
preprocessing are done in place, which means the original data is not
retrievable. Accessing ``mgmt.st_obs`` returns the processed data.

Misfit Quantification
~~~~~~~~~~~~~~~~~~~~~

Now that the data has been processed, we can run our misfit
quantification:

.. code:: ipython3

    mgmt.window()
    mgmt.measure()

``Pyatoa`` has run ``Pyflex`` and recovered 1 window for the traces,
which is saved into the Manager as a dictionary. It has then run
``Pyadjoint`` to measure the misfit on the recovered windows from
``Pyflex``. We can take a look at the ``windows`` and ``adj_srcs`` to
see that they are saved as dictionaries following the
``Config.component_list``. Each entry of the ``windows`` dictionary is a
list of ``Window`` objects from ``Pyflex``. Each entry of the
``adj_srcs`` is an AdjointSource object from ``Pyadjoint``.

.. code:: ipython3

    print(f"Windows:\n{mgmt.windows}\n")
    print(f"Adjoint Sources:\n{mgmt.adj_srcs}\n")
    print(f"Adjoint Source Object:\n{mgmt.adj_srcs['E']}\n")
    print(f"Adjoint Source Data:\n{mgmt.adj_srcs['E'].adjoint_source}")

Plotting, mapping
~~~~~~~~~~~~~~~~~

Plot functionalities will plot the available streams in the Manager
class. If windows and adjoint sources are available for a given
component, they will also be plotted alongside the waveform data. Kwargs
can be passed to the plot function to change the default look of the
waveform plots. Information about the chosen misfit windows, such as the
time shift between Observed and Synthetic traces, will be annotated into
each window.

Mapping capabilities are also available. The mapping function will try
to include as much information as possible. In this case since we have
not included an Event in our Manager class, the mapper will only be able
to plot the station. Additionally, as we have not specified any map
corners, the Map will be plotted for the whole Earth.

.. code:: ipython3

    mgmt.plot()

Tutorial 2: custom configs for Pyflex and Pyadjoint
---------------------------------------------------

Config objects are written for both
`Pyflex <https://krischer.github.io/pyflex/_modules/pyflex/config.html>`__
and
`Pyadjoint <https://github.com/krischer/pyadjoint/blob/master/src/pyadjoint/config.py>`__
so that the User could fine tune their misfit measurement criteria. For
simplicitiy, these Config classes have been wrapped into the Pyatoa
Config class, either through map names or key word arguments. We can
take a look at the default settings by printing the Config class.

.. code:: ipython3

    from pprint import pprint
    cfg = pyatoa.Config()
    print(f"Pyflex Config Map: {cfg.pyflex_map}")
    pprint(vars(cfg.pyflex_config))
    
    print('\n')
    print(f"Pyadjoint Config: {cfg.adj_src_type}")
    pprint(vars(cfg.pyadjoint_config))

Pyflex Config
~~~~~~~~~~~~~

To override the default parameters of the Pyflex Config, two options are
given.

1) hardcode the source-code of Pyatoa to include a preset Config map in
   pyatoa/plugins/pyflex_config.set_pyflex_config(). There we have
   already set three pre-set maps, “example”, “alaska” and “hikurangi”.
   Users may define their own custom choices and map names.

2) pass keyword arguments through the Pyatoa Config. These are then
   passed to the Pyflex Config if matching attributes are found.

Below we show examples of both types of overwriting. For simplicity we
just look at the Pyflex Config parameter ``c_0`` to see how this works.

.. code:: ipython3

    cfg = pyatoa.Config()
    print(f'Default c_0: {cfg.pyflex_config.c_0}')
    cfg = pyatoa.Config(pyflex_map="incorrect_map_names_will_set_the_pyflex_config_to_default")
    print(f'Incorrect map name, no kwargs c_0: {cfg.pyflex_config.c_0}')
    cfg = pyatoa.Config(pyflex_map="hikurangi")
    print(f'Option 1, preset map c_0: {cfg.pyflex_config.c_0}')
    cfg = pyatoa.Config(pyflex_map=None, c_0=2.0)
    print(f'Option 2, kwargs c_0: {cfg.pyflex_config.c_0}')

Pyadjoint Config
~~~~~~~~~~~~~~~~

To overwrite the default Config parameters of pyadjoint, kwargs may be
passed through the Pyatoa ``Config`` class to overwrite the default
Pyadjoint Config parameters. Here we show examples for the
``phase_step`` parameter in the pyadjoint.Config() class.

Note: Adjoint source types must also be passed to the Config parameter.
Pyadjoint allows custom adjoint source types, so there are no checks for
correctly specified ``adj_src_type``. Attempting to run
Manager.measure() with an undefined ``adj_src_type`` will lead to errors
from Pyadjoint.

.. code:: ipython3

    cfg = pyatoa.Config()
    print(f'Default adj_src_type: {cfg.adj_src_type}')
    cfg = pyatoa.Config(adj_src_type="custom_adj_src_type")
    print(f'Custom adj_src_type: {cfg.adj_src_type}')
    print(f'Default phase_step: {cfg.pyadjoint_config.phase_step}')
    cfg = pyatoa.Config(phase_step=3)
    print(f'Kwargs phase_step: {cfg.pyadjoint_config.phase_step}')

Tutoral 3: dynamic data gathering
---------------------------------

Pyatoa allows for dynamic data gathering based on its mid-level
``Gatherer`` class which calls on the low-level ``Getter`` and
``Fetcher`` classes. These are all wrapped up in the ``Manager`` class
so that the User does not need to interact with the lower levels.

When data gathering, ``Manager`` always seeks the path of least
resistance; that is, the ``Manager`` will always search for data
internally (either via a given Pyasdf Dataset, or through directory
structures), before moving on to more costly external searches which use
FDSN to query webservices. We can see an example of this in the
following code.

Note: Data is fetched based on event origintime information.

Internal Fetching
'''''''''''''''''

To get the ``Manager`` to search internal pathways, paths must be given
in the proper format. 1) Observed data must be saved based on SEED
observatory convention, that is in a specific directory structure with
specific file naming, following the format:

| ``/path/to/waveforms/YEAR/NETWORK/STATION/CHANNEL/NN.SSS.LL.CCC.D.YYYY.DDD``
| for example
| ``/path/to/waveforms/2018/NZ/BFZ/HHZ.D/NZ.BFZ.10.HHZ.D.2018.049``

2) Synthetic data should be placed in a directory, with the convention
   ``path/to/synthetics/NN.SSS.CCC.sem?``

.. code:: ipython3

    import logging
    from pyatoa import Config, Manager
    
    logger = logging.getLogger("pyatoa")
    logger.setLevel(logging.DEBUG)
    
    cfg = Config(event_id="2018p130600",
                 cfgpaths={"waveforms":"../tests/data/test_directories/waveforms",
                           "synthetics":"../tests/data/test_directories/synthetics"}
                )
    mgmt = Manager(cfg)
    mgmt.gather(station_code="NZ.BFZ.??.HH?", choice=["st_obs", "st_syn"])
    print(mgmt)

External Getting
''''''''''''''''

We can see that we have collected data from the given directories. Since
data gathering happens using event origin times, the Event must be
gathered before the data. If no paths are given in the
``Config.cfgpaths`` parameter, and no pyasdf.ASDFDataSet is assigned to
the ``Manager``, then the Gatherer will query the FDSN webservice in
order to download data. This has been tested with New Zealand’s GeoNet
FDSN client, and IRIS’ FDSN client, however not with other webservices,
although no problems are expected.

Moment Tensors
''''''''''''''

The Gatherer will also try to fetch moment tensor information. So far
this only works for GeoNet and IRIS events. For GeoNet events, the
github repository containing John Ristau’s moment tensors will be read
in and appended. For IRIS events, GCMT will be queried for moment tensor
information.

.. code:: ipython3

    mgmt = Manager(Config(event_id="2018p130600", client="GEONET"))
    mgmt.gather('NZ.BFZ.??.HH?', choice=["inv", "st_obs"])
    print(mgmt)

.. code:: ipython3

    # e.g. we can look at the 2018-01-23 Mww7.9 Gulf Of Alaska event,
    # recorded at the Black Forest Observatory in Germany
    mgmt = Manager(pyatoa.Config(event_id="10607586", client="IRIS"))
    mgmt.gather("II.ERM.00.BHZ")
    print(mgmt)
    mgmt.st_obs.plot()

Tutorial 4: interacting with pyasdf
-----------------------------------

Pyasdf datasets are HDF5 datasets that allow for heirarchical storage of
seismic data. It provides a very clean and compact way to store all the
data that is collected during this workflow, including raw seismic
waveforms, event and response information, misfit window and adjoint
source information. To enable saving to a Pyasdf Dataset, you simply
need to set your target dataset as an input when calling the ``Manager``
class.

.. code:: ipython3

    import pyasdf
    import logging
    from pyatoa import Config, Manager
    
    logger = logging.getLogger("pyatoa")
    logger.setLevel(logging.DEBUG)
    
    ds = pyasdf.ASDFDataSet("test_dataset.h5")
    
    cfg = Config(event_id="2018p130600", client="GEONET", 
                 cfgpaths={"synthetics":"../tests/data/test_directories/synthetics",
                           "waveforms":"../tests/data/test_directories/waveforms"}
                 )
    cfg.write(write_to=ds)
    mgmt = Manager(cfg, ds=ds)
    mgmt.gather("NZ.BFZ.??.HH?")
    print(mgmt)

.. code:: ipython3

    print(f"DATASET:\n{ds}\n\n"
          f"EVENTS:\n{ds.events}\n\n"
          f"STATION:\n{ds.waveforms['NZ.BFZ']}\n\n"
          f"INVENTORY:\n{ds.waveforms['NZ.BFZ'].StationXML}\n\n"
          f"WAVEFORMS:\n{ds.waveforms['NZ.BFZ']['observed']}"
         )


.. code:: ipython3

    mgmt.standardize()
    mgmt.preprocess()
    mgmt.window()
    mgmt.measure()

.. code:: ipython3

    print(f"AUX DATA:\n{ds.auxiliary_data}\n\n"
          f"WINDOWS:\n{ds.auxiliary_data.MisfitWindows['None'].NZ_BFZ_E_0}\n\n"
          f"ADJOINT SOURCES:\n{ds.auxiliary_data.AdjointSources['None'].NZ_BFZ_HXE}"
         )


Tutorial 6: visualization
-------------------------

To be written

Tutorial 7: export to Specfem3D
-------------------------------

To be written

Tutorial 8: plugin to Seisflows
-------------------------------

To be written

Tutorial 9: misfit statistics
-----------------------------

To be written

