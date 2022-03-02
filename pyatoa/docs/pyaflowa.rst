Pyaflowa
========

:doc:`Pyaflowa </modules/core.pyaflowa` is Pyatoa's workflow management class. It is in charge of facilitating automation of the misfit quantification machinery defined by the Config, Manager and Gatherer classes. Pyaflowa provides a structure to process en-masse waveform data from multiple stations and multiple events using only a few function calls.

``Pyaflowa`` can be used standalone, and also provides the necessary
interface to work with SeisFlows3. When used within the SeisFlows3
preprocess module, Pyaflowa reduces the overhead required to include
Pyatoa functionality into a SeisFlows3 inversion.

Standalone (w/ SPECFEM3D)
-------------------------

Let’s say you want to run forward and adjoint simulations to generate
synthetic seismograms and sensitivity kernels using SPECFEM3D. To make
this happen, a user will need some specific pieces that correspond to
the problem at hand.

-  Forward simulations: DATA/CMTSOLUTION, DATA/STATIONS
-  Adjoint simulations: DATA/STATIONS_ADJOINT, /sem/*.adj

We can use ``Pyatoa`` to prepare the necessary files for SPECFEM3D
forward simulatinos, and then use ``Pyaflowa`` to: 1) compare synthetic
seismograms with real data for a number of stations and generate adjoint
sources which can be fed directly back into a SPECFEM3D adjoint
simulation.

Prep: Generate CMTSOLUTION
~~~~~~~~~~~~~~~~~~~~~~~~~~

First we must generate the source file used to represent the moment
tensor in our earthquake simulation. In this example, we’ll use my
favorite New Zealand `example event,
2018p130600 <https://www.geonet.org.nz/earthquake/2018p130600>`__, an
M5.2 that occurred in the central North Island, New Zealand while I was
doing my PhD in Wellington.

With ``ObsPy`` + ``Pyatoa`` we’ll be able to gather this event from FDSN
webservices and save the resulting moment tensor and event information
to the required `CMTSOLUTION
format <https://www.globalcmt.org/CMTsearch.html>`__ expected by
SPECFEM3D.

We will use the Pyatoa function :doc:`append_focal_mechanism </modules/core.gatherer/append_focal_mechanism>`

.. code:: ipython3

    import os
    from pyatoa import append_focal_mechanism
    from obspy.clients.fdsn import Client

.. code:: ipython3

    # We will just move to an empty directory for this example problem
    print(os.getcwd())
    docs_path = os.path.abspath("../tests/test_data/docs_data/pyaflowa_doc")
    if not os.path.exists(docs_path):
        os.mkdir(docs_path)
    os.chdir(docs_path)
    print(os.getcwd())


.. parsed-literal::

    /home/bchow/REPOSITORIES/pyatoa/pyatoa/docs
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc


.. code:: ipython3

    # Using ObsPy's FDSN Client, we can retrieve event information using an event id
    c = Client("GEONET")
    cat = c.get_events(eventid="2018p130600")
    
    # Using Pyatoa, we can append a focal mechanism from the GeoNet moment tensor catalog
    cat[0] = append_focal_mechanism(cat[0], client="GEONET")
    
    # ObsPy has built in support for writing CMTSOLUTION files expected by SPECFEM3D
    cat.write("CMTSOLUTION", format="CMTSOLUTION")


.. parsed-literal::

    /home/bchow/miniconda3/envs/docs/lib/python3.7/site-packages/obspy/io/cmtsolution/core.py:370: UserWarning: No body wave magnitude found. Will be replaced by the first magnitude in the event object.
      warnings.warn("No body wave magnitude found. Will be replaced by the "
    /home/bchow/miniconda3/envs/docs/lib/python3.7/site-packages/obspy/io/cmtsolution/core.py:376: UserWarning: No surface wave magnitude found. Will be replaced by the first magnitude in the event object.
      warnings.warn("No surface wave magnitude found. Will be replaced by "
    /home/bchow/miniconda3/envs/docs/lib/python3.7/site-packages/obspy/io/cmtsolution/core.py:397: UserWarning: Could not find a centroid origin. Will instead assume that the preferred or first origin is the centroid origin.
      warnings.warn("Could not find a centroid origin. Will instead "
    /home/bchow/miniconda3/envs/docs/lib/python3.7/site-packages/obspy/io/cmtsolution/core.py:410: UserWarning: Hypocentral origin will be identical to the centroid one.
      warnings.warn("Hypocentral origin will be identical to the "


.. code:: ipython3

    # Lets just have a look at the file that's been created, which is a CMTSOLUTION that is ready
    # to be used in SPECFEM3D
    !cat "CMTSOLUTION"


.. parsed-literal::

     PDE 2018 02 18 07 43 48.13  -39.9490  176.2995  20.6 5.2 5.2 NORTH ISLAND, NEW ZEALAND
    event name:           785552
    time shift:           0.0000
    half duration:        0.6989
    latitude:           -39.9490
    longitude:          176.2995
    depth:               20.5946
    Mrr:           -2.479380E+23
    Mtt:            1.314880E+23
    Mpp:            1.164500E+23
    Mrt:            5.032500E+22
    Mrp:            6.607700E+22
    Mtp:            9.359300E+22


Prep: Generate STATIONS file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SPECFEM3D also requires a STATIONS file which defines the locations of
receivers for simulation output. As in Step 1 we’ll generate a list of
stations using ``ObsPy`` and write them into the required STATIONS file
using ``Pyatoa`` and a corresponding obspy.Inventory object.

   | **NOTE:**
   | In the ObsPy function get_stations(), the reasoning behind the
     following arguments provided: \* **network = “NZ”** refers to the
     code for New Zealand’s permament seismic netnwork \* **station =
     “??Z”** means we only want 3 letter station codes that end in Z,
     which GeoNet usually reserves for broadband seismometers \*
     **channel = “HH?”** refers to a broadband (first ‘H’) seismometer
     (second ‘H’), for any available component (wildcard ‘?’), usually
     N/E/Z. This follows `SEED naming
     convention <https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/>`__.
     \* This **min** and **max latitude / longitude** defines a small
     region where we want to search for stations

.. code:: ipython3

    from pyatoa import write_stations

.. code:: ipython3

    inv = c.get_stations(network="NZ", station="??Z", channel="HH?",
                         minlatitude=-41, maxlatitude=-39,
                         minlongitude=173, maxlongitude=176)
    write_stations(inv, fid="STATIONS")


.. parsed-literal::

    /home/bchow/miniconda3/envs/docs/lib/python3.7/site-packages/obspy/io/stationxml/core.py:98: UserWarning: The StationXML file has version 1, ObsPy can read versions (1.0, 1.1). Proceed with caution.
      version, ", ".join(READABLE_VERSIONS)))


.. code:: ipython3

    # Let's have a look at the stations we picked up from FDSN
    print(inv)


.. parsed-literal::

    Inventory created at 2022-03-02T06:13:37.000000Z
    	Created by: Delta
    		    
    	Sending institution: GeoNet (WEL(GNS_Test))
    	Contains:
    		Networks (1):
    			NZ
    		Stations (4):
    			NZ.MRZ (Mangatainoka River)
    			NZ.TSZ (Takapari Road)
    			NZ.VRZ (Vera Road)
    			NZ.WAZ (Wanganui)
    		Channels (0):
    


.. code:: ipython3

    # And lets have a look at the STATIONS file that's been created
    !cat "STATIONS"


.. parsed-literal::

       MRZ    NZ    -40.6605    175.5785    0.0    0.0
       TSZ    NZ    -40.0586    175.9611    0.0    0.0
       VRZ    NZ    -39.1243    174.7585    0.0    0.0
       WAZ    NZ    -39.7546    174.9855    0.0    0.0


Forward Simulation: Generate synthetics using SPECFEM3D [external]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unfortunately this cannot be shown in a Jupyter notebook as generating
synthetics requires interfacing with the SPECFEM3D code, which usually
takes place on a cluster. In this example we assume this step has been
completed successfully, with resultant synthetic waveforms produced by
SPECFEM3D for the given event and stations defined above.

   | **NOTE:**
   | Output synthetic seismograms are expected to be formatted as
     two-column ASCII files, which I have pre-generated for this
     example. File names follow the expected output from SPECFEM3D.
     Adherance to this format is very important for running Pyaflowa.

..

   **NOTE:** By default synthetic waveform data is expected to be
   separated by event ID, e.g., PATH/TO/SYNTHETICS/{EVENT_ID}/*semd

.. code:: ipython3

    # Let's copy the premade synthetic data into our current working directory
    !ls ../../synthetics
    !mkdir -p synthetics
    !cp -r ../../synthetics/201?p?????? ./synthetics


.. parsed-literal::

    2012p242656  2018p130600  NZ.BFZ.BXE.semd  NZ.BFZ.BXN.semd  NZ.BFZ.BXZ.semd


.. code:: ipython3

    !head ./synthetics/2018p130600/NZ.MRZ.BX?.semd


.. parsed-literal::

    ==> ./synthetics/2018p130600/NZ.MRZ.BXE.semd <==
      -20.0000000       0.00000000    
      -19.9850006       0.00000000    
      -19.9699993       0.00000000    
      -19.9549999       0.00000000    
      -19.9400005       0.00000000    
      -19.9249992       0.00000000    
      -19.9099998       0.00000000    
      -19.8950005       0.00000000    
      -19.8799992       0.00000000    
      -19.8649998       0.00000000    
    
    ==> ./synthetics/2018p130600/NZ.MRZ.BXN.semd <==
      -20.0000000       0.00000000    
      -19.9850006       0.00000000    
      -19.9699993       0.00000000    
      -19.9549999       0.00000000    
      -19.9400005       0.00000000    
      -19.9249992       0.00000000    
      -19.9099998       0.00000000    
      -19.8950005       0.00000000    
      -19.8799992       0.00000000    
      -19.8649998       0.00000000    
    
    ==> ./synthetics/2018p130600/NZ.MRZ.BXZ.semd <==
      -20.0000000       0.00000000    
      -19.9850006       0.00000000    
      -19.9699993       0.00000000    
      -19.9549999       0.00000000    
      -19.9400005       0.00000000    
      -19.9249992       0.00000000    
      -19.9099998       0.00000000    
      -19.8950005       0.00000000    
      -19.8799992       0.00000000    
      -19.8649998       0.00000000    


Pyaflowa’s directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Pyaflowa`` abstracts away the enigmatic inner machinations of
``Pyatoa``. To do so it manages an internal directory structure to
search for inputs and store outputs.

When used standalone, ``Pyaflowa`` creates its own directory structure
within a given working directory. When used in conjunction with
SeisFlows3, ``Pyaflowa`` will work within the preset internal directory
structure of SeisFlows3 (see Pyaflowa + SeisFlows3).

Let’s start by initiating ``Pyaflowa``. As with any usage of Pyatoa, a
Config object is required to define internally used parameters which
will inturn be used to control gathering, waveform processing, and
misfit quantification.

.. code:: ipython3

    from pyatoa import Pyaflowa, Config

.. code:: ipython3

    cfg = Config(iteration=1, step_count=0, client="GEONET", min_period=10, max_period=30,
                 pyflex_preset="nznorth_10-30s")
    
    pf = Pyaflowa(structure="standalone", workdir="./", config=cfg)

.. code:: ipython3

    # We can take a look at Pyaflowa's DEFAULT internal directory structure with the path_structure attribute
    pf.path_structure




.. parsed-literal::

    cwd          : './'
    data         : './input/DATA'
    datasets     : './datasets'
    figures      : './figures'
    logs         : './logs'
    ds_file      : './datasets/{source_name}.h5'
    stations_file: './{source_name}/STATIONS'
    responses    : './input/responses'
    waveforms    : './input/waveforms'
    synthetics   : './input/synthetics/{source_name}'
    adjsrcs      : './adjsrcs/{source_name}'
    event_figures: './figures/{source_name}'



If you want to change the Pyaflowa directory structure from the default
values shown above, you can simply pass the keys of the
``path_structure`` attribute as keyword arguments in the initialization
of Pyaflowa. Let’s generate a non-standard path structure to point to
our existing data.

.. code:: ipython3

    # Make sure we're in the correct directory so that we don't start making dir. randomly
    assert os.path.basename(os.getcwd()) == "pyaflowa_doc"
    
    # Custom set directory structure
    base_path = os.getcwd()
    kwargs = {"workdir": base_path,
              "synthetics": os.path.join(base_path, "synthetics", "{source_name}"),
              "stations_file": os.path.join(base_path, "STATIONS"),
              "data": base_path,
              "datasets": base_path,
             }
    
    pf = Pyaflowa(structure="standalone", config=cfg, **kwargs)
    pf.path_structure




.. parsed-literal::

    cwd          : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    data         : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    datasets     : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    figures      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures'
    logs         : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/logs'
    ds_file      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/{source_name}.h5'
    stations_file: '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/STATIONS'
    responses    : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/responses'
    waveforms    : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/waveforms'
    synthetics   : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/{source_name}'
    adjsrcs      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/adjsrcs/{source_name}'
    event_figures: '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/{source_name}'



--------------

The IO (input/output) class
~~~~~~~~~~~~~~~~~~~~~~~~~~~

By running the Pyaflowa.setup() function, Pyaflowa will make the
required directory structure defined above. It will also return an
``IO`` object. This internally used object store information related to
paths, configurations and processing.

The user does **not** need to interact with the ``IO`` object, but we
can take a look at it for clarity. It contains the internal directory
structure used by ``Pyaflowa``, the ``Config`` object which will control
all of the Manager processing that will take place, and internal
attributes which keep track of how processing occurs.

.. code:: ipython3

    io = pf.setup(source_name="2018p130600")

.. code:: ipython3

    for key, val in io.items():
        print(key, val)


.. parsed-literal::

    paths cwd          : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    data         : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    datasets     : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    figures      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures'
    logs         : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/logs'
    ds_file      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/2018p130600.h5'
    stations_file: '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/STATIONS'
    responses    : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/responses'
    waveforms    : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/waveforms'
    synthetics   : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600'
    adjsrcs      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/adjsrcs/2018p130600'
    event_figures: '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600'
    
    logger <Logger pyatoa (DEBUG)>
    config CONFIG
        iteration:               i01
        step_count:              s00
        event_id:                2018p130600
    GATHER
        client:                  GEONET
        start_pad:               20
        end_pad:                 500
        save_to_ds:              True
    PROCESS
        min_period:              10.0
        max_period:              30.0
        filter_corners:          2.0
        unit_output:             DISP
        rotate_to_rtz:           False
        win_amp_ratio:           0.0
        synthetics_only:         False
    LABELS
        component_list:          ['E', 'N', 'Z']
        observed_tag:            observed
        synthetic_tag:           synthetic_i01s00
        paths:                   {'responses': '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/responses', 'waveforms': '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/waveforms', 'synthetics': '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600', 'events': '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'}
    EXTERNAL
        pyflex_preset:           nznorth_10-30s
        adj_src_type:            cc_traveltime_misfit
        pyflex_config:           <pyflex.config.Config object at 0x7f04e245a910>
        pyadjoint_config:        <pyadjoint.config.Config object at 0x7f04e245a490>
    
    misfit 0
    nwin 0
    stations 0
    processed 0
    exceptions 0
    plot_fids []


--------------

Running Pyaflowa (gather and process waveforms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Great, we’re all set up to run ``Pyaflowa``. Internally ``Pyaflowa``
knows the event, path structure and stations that we want to use for
misfit quantification. Now when we run it, ``Pyaflowa`` will instantiate
``Manager`` classes, attempt to gather data from disk or from web
services, preprocess data and synthetic waveforms according to the
``Config`` object, and generate misfit windows and adjoint sources.

The example code block below is a an example of what Pyaflowa is doing
under the hood: it simply abstracts commands that are used to run
processing for multiple stations. It also contains a few internal checks
to make sure unexpected errors don’t throw the processing step off the
rails.::

::

   from pyasdf import ASDFDataSet 
   from pyatoa import Manager

   with ASDFDataSet(io.paths.ds_file) as ds:
       mgmt = Manager(ds=ds, config=io.config)
       for code in ["NZ.BFZ.*.HH?"]:
           mgmt.gather(code=code)
           mgmt.flow()

.. code:: ipython3

    pf.process_event(source_name="2018p130600", loc="*", cha="HH?")


.. parsed-literal::

    [2022-03-02 14:08:56] - pyatoa - DEBUG: gathering event
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for event info
    [2022-03-02 14:08:56] - pyatoa - DEBUG: matching event found: 81EE9F
    [2022-03-02 14:08:56] - pyatoa - INFO: 
    ================================================================================
    
    NZ.MRZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:08:56] - pyatoa - DEBUG: gathering event
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for event info
    [2022-03-02 14:08:56] - pyatoa - DEBUG: matching event found: 81EE9F
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering data for NZ.MRZ.*.HH?
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:08:56] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:08:56] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:08:56] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:08:56] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:08:56] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:08:56] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.MRZ.BXE.semd
    [2022-03-02 14:08:56] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.MRZ.BXN.semd
    [2022-03-02 14:08:56] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.MRZ.BXZ.semd
    [2022-03-02 14:08:56] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:08:56] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:08:56] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:08:56] - pyatoa - DEBUG: shifting NZ.MRZ.10.HHE starttime by 0.001607s
    [2022-03-02 14:08:56] - pyatoa - DEBUG: shifting NZ.MRZ.10.HHN starttime by 0.001607s
    [2022-03-02 14:08:56] - pyatoa - DEBUG: shifting NZ.MRZ.10.HHZ starttime by 0.001607s
    [2022-03-02 14:08:56] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:08:56] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:08:56] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:56] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:08:56] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:08:56] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:56] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:08:56] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:56] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:08:56] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:56] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:08:56] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:08:56,898] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:56,899] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:56,899] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:56,900] - pyflex - INFO: Initial window selection yielded 2 possible windows.
    [2022-03-02 14:08:56,900] - pyflex - INFO: Rejection based on travel times retained 2 windows.
    [2022-03-02 14:08:56,901] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 264017440516.493011, Amplitude SNR: 1258959.382416
    [2022-03-02 14:08:56,901] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:08:56,901] - pyflex - INFO: Water level rejection retained 1 windows
    [2022-03-02 14:08:56,902] - pyflex - INFO: Single phase group rejection retained 1 windows
    [2022-03-02 14:08:56,902] - pyflex - INFO: Removing duplicates retains 1 windows.
    [2022-03-02 14:08:56,902] - pyflex - INFO: Rejection based on minimum window length retained 1 windows.
    [2022-03-02 14:08:56,902] - pyflex - INFO: SN amplitude ratio window rejection retained 1 windows
    [2022-03-02 14:08:56,906] - pyflex - INFO: Rejection based on data fit criteria retained 1 windows.
    [2022-03-02 14:08:56,906] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:56] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:08:57,036] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:57,036] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:57,037] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:57,038] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:57,038] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:57,039] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 77923017814.226379, Amplitude SNR: 910031.449339
    [2022-03-02 14:08:57,039] - pyflex - INFO: Rejection based on minimum window length retained 10 windows.
    [2022-03-02 14:08:57,039] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:57,039] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:57,040] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:57,040] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:57,040] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:57,043] - pyflex - DEBUG: Window rejected due to CC value: 0.635155
    [2022-03-02 14:08:57,043] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:57,044] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:57] - pyatoa - INFO: 1 window(s) selected for comp N
    [2022-03-02 14:08:57,172] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:57,172] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:57,173] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:57,174] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:57,174] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:57,174] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 13630877755876006.000000, Amplitude SNR: 439117916.987834
    [2022-03-02 14:08:57,175] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:08:57,175] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:08:57,175] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:08:57,175] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:08:57,176] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:08:57,176] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:08:57,182] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:57,183] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:57] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:08:57] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:08:57] - pyatoa - INFO: 3 window(s) total found
    [2022-03-02 14:08:57] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:08:57] - pyatoa - INFO: 0.366 misfit for comp E
    [2022-03-02 14:08:57] - pyatoa - INFO: 0.154 misfit for comp N
    [2022-03-02 14:08:57] - pyatoa - INFO: 0.095 misfit for comp Z
    [2022-03-02 14:08:57] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:08:57] - pyatoa - INFO: total misfit 0.614
    [2022-03-02 14:08:57] - pyatoa - INFO: 
    
    	OBS WAVS:  3
    	SYN WAVS:  3
    	WINDOWS:   3
    	MISFIT:    0.61
    
    [2022-03-02 14:08:57] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_MRZ.pdf
    [2022-03-02 14:08:57] - pyatoa - INFO: 
    ================================================================================
    
    FINALIZE
    
    ================================================================================
    [2022-03-02 14:08:58] - pyatoa - INFO: 
    ================================================================================
    
    NZ.TSZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering data for NZ.TSZ.*.HH?
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:08:58] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:08:58] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:08:58] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:08:58] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:08:58] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:08:58] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:08:58] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:08:58] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:08:58] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.TSZ.BXE.semd
    [2022-03-02 14:08:58] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.TSZ.BXN.semd
    [2022-03-02 14:08:58] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.TSZ.BXZ.semd
    [2022-03-02 14:08:58] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:08:58] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:08:58] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:08:58] - pyatoa - DEBUG: zero pad NZ.TSZ.10.HHE (0, 0) samples
    [2022-03-02 14:08:58] - pyatoa - DEBUG: new starttime NZ.TSZ.10.HHE: 2018-02-18T07:43:28.130000Z
    [2022-03-02 14:08:58] - pyatoa - DEBUG: zero pad NZ.TSZ.10.HHN (0, 0) samples
    [2022-03-02 14:08:58] - pyatoa - DEBUG: new starttime NZ.TSZ.10.HHN: 2018-02-18T07:43:28.130001Z
    [2022-03-02 14:08:58] - pyatoa - DEBUG: zero pad NZ.TSZ.10.HHZ (0, 0) samples
    [2022-03-02 14:08:58] - pyatoa - DEBUG: new starttime NZ.TSZ.10.HHZ: 2018-02-18T07:43:28.130001Z
    [2022-03-02 14:08:58] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:08:58] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:08:58] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:58] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:08:58] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:08:58] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:58] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:08:58] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:58] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:08:58] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:58] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:08:58] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:08:58,492] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:58,493] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:58,494] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:58] - pyatoa - WARNING: Cannot window, most likely because the source-receiver distance is too small w.r.t the minimum period
    [2022-03-02 14:08:58] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_TSZ.pdf
    [2022-03-02 14:08:59] - pyatoa - INFO: 
    ================================================================================
    
    NZ.VRZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering data for NZ.VRZ.*.HH?
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:08:59] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:08:59] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:08:59] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:08:59] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:08:59] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:08:59] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:08:59] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:08:59] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:08:59] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.VRZ.BXE.semd
    [2022-03-02 14:08:59] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.VRZ.BXN.semd
    [2022-03-02 14:08:59] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.VRZ.BXZ.semd
    [2022-03-02 14:08:59] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:08:59] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:08:59] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:08:59] - pyatoa - DEBUG: shifting NZ.VRZ.10.HHE starttime by 0.001578s
    [2022-03-02 14:08:59] - pyatoa - DEBUG: shifting NZ.VRZ.10.HHN starttime by 0.001578s
    [2022-03-02 14:08:59] - pyatoa - DEBUG: shifting NZ.VRZ.10.HHZ starttime by 0.001578s
    [2022-03-02 14:08:59] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:08:59] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:08:59] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:59] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:08:59] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:08:59] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:59] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:08:59] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:59] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:08:59] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:59] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:08:59] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:08:59,612] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:59,612] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:59,613] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:59,614] - pyflex - INFO: Initial window selection yielded 13 possible windows.
    [2022-03-02 14:08:59,614] - pyflex - INFO: Rejection based on travel times retained 13 windows.
    [2022-03-02 14:08:59,615] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 23466401785.970848, Amplitude SNR: 384703.121688
    [2022-03-02 14:08:59,615] - pyflex - INFO: Rejection based on minimum window length retained 13 windows.
    [2022-03-02 14:08:59,615] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:59,616] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:59,616] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:59,616] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:59,617] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:59,622] - pyflex - DEBUG: Window rejected due to CC value: 0.675487
    [2022-03-02 14:08:59,622] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:59,623] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:59] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:08:59,793] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:59,793] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:59,794] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:59,795] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:59,795] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:59,796] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 31526068792.362656, Amplitude SNR: 495295.417951
    [2022-03-02 14:08:59,796] - pyflex - INFO: Rejection based on minimum window length retained 10 windows.
    [2022-03-02 14:08:59,796] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:59,796] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:59,797] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:59,797] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:59,797] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:59,801] - pyflex - INFO: Rejection based on data fit criteria retained 3 windows.
    [2022-03-02 14:08:59,802] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:59] - pyatoa - INFO: 1 window(s) selected for comp N
    [2022-03-02 14:08:59,972] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:59,973] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:59,973] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:59,974] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:59,975] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:59,975] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 26562645048.905380, Amplitude SNR: 378503.886560
    [2022-03-02 14:08:59,975] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:08:59,975] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:08:59,976] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:08:59,976] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:08:59,976] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:08:59,977] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:08:59,987] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:59,987] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:59] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:08:59] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:08:59] - pyatoa - INFO: 3 window(s) total found
    [2022-03-02 14:08:59] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:09:00] - pyatoa - INFO: 0.198 misfit for comp E
    [2022-03-02 14:09:00] - pyatoa - INFO: 0.011 misfit for comp N
    [2022-03-02 14:09:00] - pyatoa - INFO: 0.065 misfit for comp Z
    [2022-03-02 14:09:00] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:09:00] - pyatoa - INFO: total misfit 0.275
    [2022-03-02 14:09:00] - pyatoa - INFO: 
    
    	OBS WAVS:  3
    	SYN WAVS:  3
    	WINDOWS:   3
    	MISFIT:    0.27
    
    [2022-03-02 14:09:00] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_VRZ.pdf
    [2022-03-02 14:09:00] - pyatoa - INFO: 
    ================================================================================
    
    FINALIZE
    
    ================================================================================
    [2022-03-02 14:09:00] - pyatoa - INFO: 
    ================================================================================
    
    NZ.WAZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering data for NZ.WAZ.*.HH?
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:09:00] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:09:00] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:09:00] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:09:00] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:09:00] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:09:00] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:09:00] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:09:00] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:09:00] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.WAZ.BXE.semd
    [2022-03-02 14:09:00] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.WAZ.BXN.semd
    [2022-03-02 14:09:00] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.WAZ.BXZ.semd
    [2022-03-02 14:09:00] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:09:00] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:09:00] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:09:01] - pyatoa - DEBUG: shifting NZ.WAZ.10.HHE starttime by 0.001607s
    [2022-03-02 14:09:01] - pyatoa - DEBUG: shifting NZ.WAZ.10.HHN starttime by 0.001607s
    [2022-03-02 14:09:01] - pyatoa - DEBUG: shifting NZ.WAZ.10.HHZ starttime by 0.001607s
    [2022-03-02 14:09:01] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:09:01] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:09:01] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:01] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:09:01] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:09:01] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:01] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:09:01] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:01] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:09:01] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:01] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:09:01] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:09:01,233] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:01,233] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:01,234] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:01,235] - pyflex - INFO: Initial window selection yielded 22 possible windows.
    [2022-03-02 14:09:01,235] - pyflex - INFO: Rejection based on travel times retained 22 windows.
    [2022-03-02 14:09:01,236] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 28382657308.924114, Amplitude SNR: 494242.727553
    [2022-03-02 14:09:01,236] - pyflex - INFO: Rejection based on minimum window length retained 20 windows.
    [2022-03-02 14:09:01,236] - pyflex - INFO: Water level rejection retained 8 windows
    [2022-03-02 14:09:01,237] - pyflex - INFO: Single phase group rejection retained 8 windows
    [2022-03-02 14:09:01,237] - pyflex - INFO: Removing duplicates retains 4 windows.
    [2022-03-02 14:09:01,237] - pyflex - INFO: Rejection based on minimum window length retained 4 windows.
    [2022-03-02 14:09:01,238] - pyflex - INFO: SN amplitude ratio window rejection retained 4 windows
    [2022-03-02 14:09:01,246] - pyflex - INFO: Rejection based on data fit criteria retained 4 windows.
    [2022-03-02 14:09:01,247] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:01] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:09:01,382] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:01,382] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:01,383] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:01,384] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:01,384] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:01,385] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 60358322015.940735, Amplitude SNR: 576780.923206
    [2022-03-02 14:09:01,385] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:01,385] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:01,385] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:01,386] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:01,386] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:01,386] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:01,389] - pyflex - DEBUG: Window rejected due to amplitude fit: 3.871042
    [2022-03-02 14:09:01,390] - pyflex - DEBUG: Window rejected due to amplitude fit: 3.793857
    [2022-03-02 14:09:01,390] - pyflex - INFO: Rejection based on data fit criteria retained 0 windows.
    [2022-03-02 14:09:01,390] - pyflex - INFO: Weighted interval schedule optimization retained 0 windows.
    [2022-03-02 14:09:01] - pyatoa - INFO: 0 window(s) selected for comp N
    [2022-03-02 14:09:01,567] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:01,568] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:01,568] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:01,569] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:01,570] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:01,570] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 4763134899.314076, Amplitude SNR: 166126.234705
    [2022-03-02 14:09:01,570] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:01,571] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:01,571] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:01,571] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:01,571] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:01,572] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:01,578] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:09:01,578] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:01] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:09:01] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:09:01] - pyatoa - INFO: 2 window(s) total found
    [2022-03-02 14:09:01] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:09:01] - pyatoa - INFO: 1.037 misfit for comp E
    [2022-03-02 14:09:01] - pyatoa - INFO: 0.293 misfit for comp Z
    [2022-03-02 14:09:01] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:09:01] - pyatoa - INFO: total misfit 1.329
    [2022-03-02 14:09:01] - pyatoa - INFO: 
    
    	OBS WAVS:  3
    	SYN WAVS:  3
    	WINDOWS:   2
    	MISFIT:    1.33
    
    [2022-03-02 14:09:01] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_WAZ.pdf
    [2022-03-02 14:09:02] - pyatoa - INFO: 
    ================================================================================
    
    FINALIZE
    
    ================================================================================
    [2022-03-02 14:09:02] - pyatoa - INFO: creating single .pdf file of all output figures
    [2022-03-02 14:09:02] - pyatoa - INFO: generating STATIONS_ADJOINT file for SPECFEM
    [2022-03-02 14:09:02] - pyatoa - INFO: 
    ================================================================================
    
    SUMMARY
    
    ================================================================================
    SOURCE NAME: 2018p130600
    STATIONS: 3 / 4
    WINDOWS: 8
    RAW MISFIT: 2.22
    UNEXPECTED ERRORS: 0




.. parsed-literal::

    0.138628125



Inspect Pyaflowa outputs
~~~~~~~~~~~~~~~~~~~~~~~~

Iwe have a look at the work directory, we can see the outputs of the
Pyaflowa workflow, which will be: \* An ASDFDataSet with waveforms,
metadata, misfit windows and adjoint sources \* Waveform figures for all
the stations processed \* Adjoint source ASCII files (.adj) required for
a SPECFEM3D adjoint simulation \* STATIONS_ADJOINT file required for a
SPECFEM3D adjoint simulation \* The output log which shows the

.. code:: ipython3

    # Here is the working directory with all the inputs and outputs
    !ls


.. parsed-literal::

    2012p242656.h5	adjsrcs      figures  logs	STATIONS_ADJOINT
    2018p130600.h5	CMTSOLUTION  input    STATIONS	synthetics


The ASDFDataSet contains all the data and metadata collected and created during the workflow. This can be viewed using the functionalities of PyASDF, which is detailed further in the :doc:`Data Storage </storage>` documentation page.

.. code:: ipython3

    # Each event will output adjoint source files that can be fed directly into SPECFEM3D
    !ls adjsrcs/2018p130600


.. parsed-literal::

    NZ.BFZ.BXE.adj	NZ.MRZ.BXE.adj	NZ.VRZ.BXE.adj	NZ.WAZ.BXE.adj
    NZ.BFZ.BXN.adj	NZ.MRZ.BXN.adj	NZ.VRZ.BXN.adj	NZ.WAZ.BXN.adj
    NZ.BFZ.BXZ.adj	NZ.MRZ.BXZ.adj	NZ.VRZ.BXZ.adj	NZ.WAZ.BXZ.adj


.. code:: ipython3

    # Adjoint source files are created as two-column ASCII files, in the same manner as the synthetics 
    # generated by SPECFEM3D
    !head adjsrcs/2018p130600/NZ.MRZ.*.adj


.. parsed-literal::

    ==> adjsrcs/2018p130600/NZ.MRZ.BXE.adj <==
    -2.000000000000000000e+01 0.000000000000000000e+00
    -1.998499999999999943e+01 0.000000000000000000e+00
    -1.996999999999999886e+01 0.000000000000000000e+00
    -1.995499999999999829e+01 0.000000000000000000e+00
    -1.994000000000000128e+01 0.000000000000000000e+00
    -1.992500000000000071e+01 0.000000000000000000e+00
    -1.991000000000000014e+01 0.000000000000000000e+00
    -1.989499999999999957e+01 0.000000000000000000e+00
    -1.987999999999999901e+01 0.000000000000000000e+00
    -1.986499999999999844e+01 0.000000000000000000e+00
    
    ==> adjsrcs/2018p130600/NZ.MRZ.BXN.adj <==
    -2.000000000000000000e+01 0.000000000000000000e+00
    -1.998499999999999943e+01 0.000000000000000000e+00
    -1.996999999999999886e+01 0.000000000000000000e+00
    -1.995499999999999829e+01 0.000000000000000000e+00
    -1.994000000000000128e+01 0.000000000000000000e+00
    -1.992500000000000071e+01 0.000000000000000000e+00
    -1.991000000000000014e+01 0.000000000000000000e+00
    -1.989499999999999957e+01 0.000000000000000000e+00
    -1.987999999999999901e+01 0.000000000000000000e+00
    -1.986499999999999844e+01 0.000000000000000000e+00
    
    ==> adjsrcs/2018p130600/NZ.MRZ.BXZ.adj <==
    -2.000000000000000000e+01 0.000000000000000000e+00
    -1.998499999999999943e+01 0.000000000000000000e+00
    -1.996999999999999886e+01 0.000000000000000000e+00
    -1.995499999999999829e+01 0.000000000000000000e+00
    -1.994000000000000128e+01 0.000000000000000000e+00
    -1.992500000000000071e+01 0.000000000000000000e+00
    -1.991000000000000014e+01 0.000000000000000000e+00
    -1.989499999999999957e+01 0.000000000000000000e+00
    -1.987999999999999901e+01 0.000000000000000000e+00
    -1.986499999999999844e+01 0.000000000000000000e+00


.. code:: ipython3

    # A composite PDF of all waveform figures for each source-receiver pair will be generated
    # for the user to quickly evaluate data-synthetic misfit graphically
    !ls figures/2018p130600


.. parsed-literal::

    i01s00_2018p130600.pdf


.. code:: ipython3

    # This doesn't work, image needs to be stored within the dir. that I opened jupyter notebook with
    # from IPython.display import IFrame
    # IFrame("figures/2018p130600/i01s00_2018p130600.pdf", width=600, height=300)

.. code:: ipython3

    # Text log files help the user keep track of all processing steps, and misfit information
    !ls logs


.. parsed-literal::

    i01s00_2012p242656.log	i01s00_2018p130600.log


.. code:: ipython3

    !cat logs/i01s00_2018p130600.log


.. parsed-literal::

    [2022-03-02 14:08:56] - pyatoa - INFO: 
    ================================================================================
    
    NZ.MRZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:08:56] - pyatoa - DEBUG: gathering event
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for event info
    [2022-03-02 14:08:56] - pyatoa - DEBUG: matching event found: 81EE9F
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering data for NZ.MRZ.*.HH?
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:08:56] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:08:56] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:08:56] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:08:56] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:08:56] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:08:56] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:08:56] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:08:56] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.MRZ.BXE.semd
    [2022-03-02 14:08:56] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.MRZ.BXN.semd
    [2022-03-02 14:08:56] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.MRZ.BXZ.semd
    [2022-03-02 14:08:56] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:08:56] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:08:56] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:08:56] - pyatoa - DEBUG: shifting NZ.MRZ.10.HHE starttime by 0.001607s
    [2022-03-02 14:08:56] - pyatoa - DEBUG: shifting NZ.MRZ.10.HHN starttime by 0.001607s
    [2022-03-02 14:08:56] - pyatoa - DEBUG: shifting NZ.MRZ.10.HHZ starttime by 0.001607s
    [2022-03-02 14:08:56] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:08:56] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:08:56] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:56] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:08:56] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:08:56] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:56] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:08:56] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:56] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:08:56] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:56] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:08:56] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:08:56] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:56] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:56] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:56] - pyflex - INFO: Initial window selection yielded 2 possible windows.
    [2022-03-02 14:08:56] - pyflex - INFO: Rejection based on travel times retained 2 windows.
    [2022-03-02 14:08:56] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 264017440516.493011, Amplitude SNR: 1258959.382416
    [2022-03-02 14:08:56] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:08:56] - pyflex - INFO: Water level rejection retained 1 windows
    [2022-03-02 14:08:56] - pyflex - INFO: Single phase group rejection retained 1 windows
    [2022-03-02 14:08:56] - pyflex - INFO: Removing duplicates retains 1 windows.
    [2022-03-02 14:08:56] - pyflex - INFO: Rejection based on minimum window length retained 1 windows.
    [2022-03-02 14:08:56] - pyflex - INFO: SN amplitude ratio window rejection retained 1 windows
    [2022-03-02 14:08:56] - pyflex - INFO: Rejection based on data fit criteria retained 1 windows.
    [2022-03-02 14:08:56] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:56] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:08:57] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:57] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:57] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:57] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 77923017814.226379, Amplitude SNR: 910031.449339
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on minimum window length retained 10 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:57] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:57] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:57] - pyflex - DEBUG: Window rejected due to CC value: 0.635155
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:57] - pyatoa - INFO: 1 window(s) selected for comp N
    [2022-03-02 14:08:57] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:57] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:57] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:57] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 13630877755876006.000000, Amplitude SNR: 439117916.987834
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:08:57] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:08:57] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:08:57] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:57] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:57] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:08:57] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:08:57] - pyatoa - INFO: 3 window(s) total found
    [2022-03-02 14:08:57] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:08:57] - pyatoa - INFO: 0.366 misfit for comp E
    [2022-03-02 14:08:57] - pyatoa - INFO: 0.154 misfit for comp N
    [2022-03-02 14:08:57] - pyatoa - INFO: 0.095 misfit for comp Z
    [2022-03-02 14:08:57] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:08:57] - pyatoa - INFO: total misfit 0.614
    [2022-03-02 14:08:57] - pyatoa - INFO: 
    
    	OBS WAVS:  3
    	SYN WAVS:  3
    	WINDOWS:   3
    	MISFIT:    0.61
    
    [2022-03-02 14:08:57] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_MRZ.pdf
    [2022-03-02 14:08:57] - pyatoa - INFO: 
    ================================================================================
    
    FINALIZE
    
    ================================================================================
    [2022-03-02 14:08:58] - pyatoa - INFO: 
    ================================================================================
    
    NZ.TSZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering data for NZ.TSZ.*.HH?
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:08:58] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:08:58] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:08:58] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:08:58] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:08:58] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:08:58] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:08:58] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:08:58] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:08:58] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:08:58] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.TSZ.BXE.semd
    [2022-03-02 14:08:58] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.TSZ.BXN.semd
    [2022-03-02 14:08:58] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.TSZ.BXZ.semd
    [2022-03-02 14:08:58] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:08:58] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:08:58] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:08:58] - pyatoa - DEBUG: zero pad NZ.TSZ.10.HHE (0, 0) samples
    [2022-03-02 14:08:58] - pyatoa - DEBUG: new starttime NZ.TSZ.10.HHE: 2018-02-18T07:43:28.130000Z
    [2022-03-02 14:08:58] - pyatoa - DEBUG: zero pad NZ.TSZ.10.HHN (0, 0) samples
    [2022-03-02 14:08:58] - pyatoa - DEBUG: new starttime NZ.TSZ.10.HHN: 2018-02-18T07:43:28.130001Z
    [2022-03-02 14:08:58] - pyatoa - DEBUG: zero pad NZ.TSZ.10.HHZ (0, 0) samples
    [2022-03-02 14:08:58] - pyatoa - DEBUG: new starttime NZ.TSZ.10.HHZ: 2018-02-18T07:43:28.130001Z
    [2022-03-02 14:08:58] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:08:58] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:08:58] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:58] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:08:58] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:08:58] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:58] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:08:58] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:58] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:08:58] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:58] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:08:58] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:08:58] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:58] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:58] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:58] - pyatoa - WARNING: Cannot window, most likely because the source-receiver distance is too small w.r.t the minimum period
    [2022-03-02 14:08:58] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_TSZ.pdf
    [2022-03-02 14:08:59] - pyatoa - INFO: 
    ================================================================================
    
    NZ.VRZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering data for NZ.VRZ.*.HH?
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:08:59] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:08:59] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:08:59] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:08:59] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:08:59] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:08:59] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:08:59] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:08:59] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:08:59] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:08:59] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.VRZ.BXE.semd
    [2022-03-02 14:08:59] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.VRZ.BXN.semd
    [2022-03-02 14:08:59] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.VRZ.BXZ.semd
    [2022-03-02 14:08:59] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:08:59] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:08:59] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:08:59] - pyatoa - DEBUG: shifting NZ.VRZ.10.HHE starttime by 0.001578s
    [2022-03-02 14:08:59] - pyatoa - DEBUG: shifting NZ.VRZ.10.HHN starttime by 0.001578s
    [2022-03-02 14:08:59] - pyatoa - DEBUG: shifting NZ.VRZ.10.HHZ starttime by 0.001578s
    [2022-03-02 14:08:59] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:08:59] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:08:59] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:59] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:08:59] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:08:59] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:59] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:08:59] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:59] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:08:59] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:59] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:08:59] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:08:59] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:59] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:59] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:59] - pyflex - INFO: Initial window selection yielded 13 possible windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on travel times retained 13 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 23466401785.970848, Amplitude SNR: 384703.121688
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on minimum window length retained 13 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:59] - pyflex - DEBUG: Window rejected due to CC value: 0.675487
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:59] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:08:59] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:59] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:59] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:59] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 31526068792.362656, Amplitude SNR: 495295.417951
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on minimum window length retained 10 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on data fit criteria retained 3 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:59] - pyatoa - INFO: 1 window(s) selected for comp N
    [2022-03-02 14:08:59] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:59] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:59] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:59] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 26562645048.905380, Amplitude SNR: 378503.886560
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:08:59] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:59] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:59] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:08:59] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:08:59] - pyatoa - INFO: 3 window(s) total found
    [2022-03-02 14:08:59] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:09:00] - pyatoa - INFO: 0.198 misfit for comp E
    [2022-03-02 14:09:00] - pyatoa - INFO: 0.011 misfit for comp N
    [2022-03-02 14:09:00] - pyatoa - INFO: 0.065 misfit for comp Z
    [2022-03-02 14:09:00] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:09:00] - pyatoa - INFO: total misfit 0.275
    [2022-03-02 14:09:00] - pyatoa - INFO: 
    
    	OBS WAVS:  3
    	SYN WAVS:  3
    	WINDOWS:   3
    	MISFIT:    0.27
    
    [2022-03-02 14:09:00] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_VRZ.pdf
    [2022-03-02 14:09:00] - pyatoa - INFO: 
    ================================================================================
    
    FINALIZE
    
    ================================================================================
    [2022-03-02 14:09:00] - pyatoa - INFO: 
    ================================================================================
    
    NZ.WAZ.*.HH?
    
    ================================================================================
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering data for NZ.WAZ.*.HH?
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering observed waveforms
    [2022-03-02 14:09:00] - pyatoa - INFO: searching ASDFDataSet for observations
    [2022-03-02 14:09:00] - pyatoa - INFO: matching observed waveforms found
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering StationXML
    [2022-03-02 14:09:00] - pyatoa - INFO: searching ASDFDataSet for station info
    [2022-03-02 14:09:00] - pyatoa - INFO: matching StationXML found
    [2022-03-02 14:09:00] - pyatoa - INFO: saved to ASDFDataSet
    [2022-03-02 14:09:00] - pyatoa - INFO: gathering synthetic waveforms
    [2022-03-02 14:09:00] - pyatoa - INFO: searching ASDFDataSet for synthetics
    [2022-03-02 14:09:00] - pyatoa - INFO: searching local filesystem for synthetics
    [2022-03-02 14:09:00] - pyatoa - DEBUG: searching for synthetics: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/{net}.{sta}.*{cmp}.sem{dva}
    [2022-03-02 14:09:00] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.WAZ.BXE.semd
    [2022-03-02 14:09:00] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.WAZ.BXN.semd
    [2022-03-02 14:09:00] - pyatoa - INFO: retrieved synthetics locally:
    /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/2018p130600/NZ.WAZ.BXZ.semd
    [2022-03-02 14:09:00] - pyatoa - INFO: matching synthetic waveforms found
    [2022-03-02 14:09:00] - pyatoa - INFO: saved to ASDFDataSet with tag 'synthetic_i01s00'
    [2022-03-02 14:09:00] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:09:01] - pyatoa - DEBUG: shifting NZ.WAZ.10.HHE starttime by 0.001607s
    [2022-03-02 14:09:01] - pyatoa - DEBUG: shifting NZ.WAZ.10.HHN starttime by 0.001607s
    [2022-03-02 14:09:01] - pyatoa - DEBUG: shifting NZ.WAZ.10.HHZ starttime by 0.001607s
    [2022-03-02 14:09:01] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:09:01] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:09:01] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:01] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:09:01] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:09:01] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:01] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:09:01] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:01] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:09:01] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:01] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:09:01] - pyatoa - INFO: running Pyflex w/ map: nznorth_10-30s
    [2022-03-02 14:09:01] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:01] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:01] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:01] - pyflex - INFO: Initial window selection yielded 22 possible windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on travel times retained 22 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 28382657308.924114, Amplitude SNR: 494242.727553
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on minimum window length retained 20 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Water level rejection retained 8 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Single phase group rejection retained 8 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Removing duplicates retains 4 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on minimum window length retained 4 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: SN amplitude ratio window rejection retained 4 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on data fit criteria retained 4 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:01] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:09:01] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:01] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:01] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:01] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 60358322015.940735, Amplitude SNR: 576780.923206
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:01] - pyflex - DEBUG: Window rejected due to amplitude fit: 3.871042
    [2022-03-02 14:09:01] - pyflex - DEBUG: Window rejected due to amplitude fit: 3.793857
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on data fit criteria retained 0 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Weighted interval schedule optimization retained 0 windows.
    [2022-03-02 14:09:01] - pyatoa - INFO: 0 window(s) selected for comp N
    [2022-03-02 14:09:01] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:01] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:01] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:01] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 4763134899.314076, Amplitude SNR: 166126.234705
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:01] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:09:01] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:01] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:09:01] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:09:01] - pyatoa - INFO: 2 window(s) total found
    [2022-03-02 14:09:01] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:09:01] - pyatoa - INFO: 1.037 misfit for comp E
    [2022-03-02 14:09:01] - pyatoa - INFO: 0.293 misfit for comp Z
    [2022-03-02 14:09:01] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:09:01] - pyatoa - INFO: total misfit 1.329
    [2022-03-02 14:09:01] - pyatoa - INFO: 
    
    	OBS WAVS:  3
    	SYN WAVS:  3
    	WINDOWS:   2
    	MISFIT:    1.33
    
    [2022-03-02 14:09:01] - pyatoa - INFO: saving figure to: /home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/2018p130600/i01_s00_NZ_WAZ.pdf
    [2022-03-02 14:09:02] - pyatoa - INFO: 
    ================================================================================
    
    FINALIZE
    
    ================================================================================
    [2022-03-02 14:09:02] - pyatoa - INFO: creating single .pdf file of all output figures
    [2022-03-02 14:09:02] - pyatoa - INFO: generating STATIONS_ADJOINT file for SPECFEM
    [2022-03-02 14:09:02] - pyatoa - INFO: 
    ================================================================================
    
    SUMMARY
    
    ================================================================================
    SOURCE NAME: 2018p130600
    STATIONS: 3 / 4
    WINDOWS: 8
    RAW MISFIT: 2.22
    UNEXPECTED ERRORS: 0
    ===================
    
    SUMMARY
    
    ================================================================================
    SOURCE NAME: 2018p130600
    STATIONS: 3 / 4
    WINDOWS: 8
    RAW MISFIT: 2.22
    UNEXPECTED ERRORS: 0


Multi-processing with Pyaflowa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pyaflowa allows for multi-processing using Python’s concurrent.futures.
This means that multiple events can be processed in parallel,
potentially allowing for large speed up when running waveform
processing. We’ll include another event
`2012p242656 <https://www.geonet.org.nz/earthquake/2012p242656>`__, in
our work directory and show how simple it is to use this functionality.

.. code:: ipython3

    !ls synthetics/2012p242656


.. parsed-literal::

    NZ.MRZ.BXE.semd  NZ.TSZ.BXE.semd  NZ.VRZ.BXE.semd  NZ.WAZ.BXE.semd
    NZ.MRZ.BXN.semd  NZ.TSZ.BXN.semd  NZ.VRZ.BXN.semd  NZ.WAZ.BXN.semd
    NZ.MRZ.BXZ.semd  NZ.TSZ.BXZ.semd  NZ.VRZ.BXZ.semd  NZ.WAZ.BXZ.semd


.. code:: ipython3

    # Here we initiate almost the same key word arguments, however we add a formatting statement in the synthetics
    # kwarg so that Pyaflowa knows which directory to search when looking for synthetics.
    # Here the STATION FILE is the same
    pf.path_structure




.. parsed-literal::

    cwd          : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    data         : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    datasets     : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc'
    figures      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures'
    logs         : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/logs'
    ds_file      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/{source_name}.h5'
    stations_file: '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/STATIONS'
    responses    : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/responses'
    waveforms    : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/input/waveforms'
    synthetics   : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/synthetics/{source_name}'
    adjsrcs      : '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/adjsrcs/{source_name}'
    event_figures: '/home/bchow/REPOSITORIES/pyatoa/pyatoa/tests/test_data/docs_data/pyaflowa_doc/figures/{source_name}'



.. code:: ipython3

    # Now, all we need to do to invoke multiprocessing is input a list of source names to the run function
    pf.multi_event_process(source_names=["2018p130600", "2012p242656"], cha="HH?", loc="*")


.. parsed-literal::

    Beginning parallel processing of 2 events...


.. parsed-literal::

    [2022-03-02 14:09:03] - pyatoa - DEBUG: gathering event
    [2022-03-02 14:09:03] - pyatoa - INFO: searching ASDFDataSet for event info
    [2022-03-02 14:09:03] - pyatoa - DEBUG: gathering event
    [2022-03-02 14:09:03] - pyatoa - INFO: searching ASDFDataSet for event info
    [2022-03-02 14:09:03] - pyatoa - DEBUG: matching event found: 81EE9F
    [2022-03-02 14:09:03] - pyatoa - DEBUG: matching event found: 81EE9F
    [2022-03-02 14:09:03,733] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:03,734] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:03,735] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:03,736] - pyflex - INFO: Initial window selection yielded 2 possible windows.
    [2022-03-02 14:09:03,736] - pyflex - INFO: Rejection based on travel times retained 2 windows.
    [2022-03-02 14:09:03,737] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 264017440516.493011, Amplitude SNR: 1258959.382416
    [2022-03-02 14:09:03,737] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:03,738] - pyflex - INFO: Water level rejection retained 1 windows
    [2022-03-02 14:09:03,738] - pyflex - INFO: Single phase group rejection retained 1 windows
    [2022-03-02 14:09:03,738] - pyflex - INFO: Removing duplicates retains 1 windows.
    [2022-03-02 14:09:03,739] - pyflex - INFO: Rejection based on minimum window length retained 1 windows.
    [2022-03-02 14:09:03,739] - pyflex - INFO: SN amplitude ratio window rejection retained 1 windows
    [2022-03-02 14:09:03,743] - pyflex - INFO: Rejection based on data fit criteria retained 1 windows.
    [2022-03-02 14:09:03,744] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:03,875] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:03,876] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:03,877] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:03,878] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:03,878] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:03,879] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 77923017814.226379, Amplitude SNR: 910031.449339
    [2022-03-02 14:09:03,879] - pyflex - INFO: Rejection based on minimum window length retained 10 windows.
    [2022-03-02 14:09:03,879] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:09:03,880] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:09:03,880] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:09:03,880] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:09:03,881] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:09:03,884] - pyflex - DEBUG: Window rejected due to CC value: 0.635155
    [2022-03-02 14:09:03,885] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:09:03,885] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:04,013] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:04,014] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:04,015] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:04,016] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:04,016] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:04,017] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 13630877755876006.000000, Amplitude SNR: 439117916.987834
    [2022-03-02 14:09:04,017] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:04,018] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:04,018] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:04,018] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:04,019] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:04,019] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:04,026] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:09:04,026] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:05,248] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:05,248] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:05,249] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:06,410] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:06,411] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:06,412] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:06,413] - pyflex - INFO: Initial window selection yielded 13 possible windows.
    [2022-03-02 14:09:06,413] - pyflex - INFO: Rejection based on travel times retained 13 windows.
    [2022-03-02 14:09:06,414] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 23466401785.970848, Amplitude SNR: 384703.121688
    [2022-03-02 14:09:06,414] - pyflex - INFO: Rejection based on minimum window length retained 13 windows.
    [2022-03-02 14:09:06,414] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:09:06,415] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:09:06,415] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:09:06,416] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:09:06,416] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:09:06,422] - pyflex - DEBUG: Window rejected due to CC value: 0.675487
    [2022-03-02 14:09:06,423] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:09:06,423] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:06,595] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:06,596] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:06,597] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:06,598] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:06,598] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:06,599] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 31526068792.362656, Amplitude SNR: 495295.417951
    [2022-03-02 14:09:06,599] - pyflex - INFO: Rejection based on minimum window length retained 10 windows.
    [2022-03-02 14:09:06,599] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:09:06,600] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:09:06,600] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:09:06,600] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:09:06,601] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:09:06,605] - pyflex - INFO: Rejection based on data fit criteria retained 3 windows.
    [2022-03-02 14:09:06,606] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:06,775] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:06,776] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:06,777] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:06,778] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:06,779] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:06,779] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 26562645048.905380, Amplitude SNR: 378503.886560
    [2022-03-02 14:09:06,779] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:06,780] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:06,780] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:06,781] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:06,781] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:06,781] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:06,791] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:09:06,791] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:08,021] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:08,022] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:08,023] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:08,024] - pyflex - INFO: Initial window selection yielded 22 possible windows.
    [2022-03-02 14:09:08,025] - pyflex - INFO: Rejection based on travel times retained 22 windows.
    [2022-03-02 14:09:08,025] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 28382657308.924114, Amplitude SNR: 494242.727553
    [2022-03-02 14:09:08,026] - pyflex - INFO: Rejection based on minimum window length retained 20 windows.
    [2022-03-02 14:09:08,026] - pyflex - INFO: Water level rejection retained 8 windows
    [2022-03-02 14:09:08,026] - pyflex - INFO: Single phase group rejection retained 8 windows
    [2022-03-02 14:09:08,027] - pyflex - INFO: Removing duplicates retains 4 windows.
    [2022-03-02 14:09:08,027] - pyflex - INFO: Rejection based on minimum window length retained 4 windows.
    [2022-03-02 14:09:08,028] - pyflex - INFO: SN amplitude ratio window rejection retained 4 windows
    [2022-03-02 14:09:08,037] - pyflex - INFO: Rejection based on data fit criteria retained 4 windows.
    [2022-03-02 14:09:08,037] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:09:08,173] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:08,173] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:08,174] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:08,175] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:08,176] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:08,176] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 60358322015.940735, Amplitude SNR: 576780.923206
    [2022-03-02 14:09:08,177] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:08,177] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:08,177] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:08,178] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:08,178] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:08,178] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:08,182] - pyflex - DEBUG: Window rejected due to amplitude fit: 3.871042
    [2022-03-02 14:09:08,183] - pyflex - DEBUG: Window rejected due to amplitude fit: 3.793857
    [2022-03-02 14:09:08,183] - pyflex - INFO: Rejection based on data fit criteria retained 0 windows.
    [2022-03-02 14:09:08,183] - pyflex - INFO: Weighted interval schedule optimization retained 0 windows.
    [2022-03-02 14:09:08,320] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:09:08,321] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:09:08,322] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:09:08,323] - pyflex - INFO: Initial window selection yielded 10 possible windows.
    [2022-03-02 14:09:08,324] - pyflex - INFO: Rejection based on travel times retained 10 windows.
    [2022-03-02 14:09:08,324] - pyflex - INFO: Global SNR checks passed. Integrated SNR: 4763134899.314076, Amplitude SNR: 166126.234705
    [2022-03-02 14:09:08,325] - pyflex - INFO: Rejection based on minimum window length retained 9 windows.
    [2022-03-02 14:09:08,325] - pyflex - INFO: Water level rejection retained 3 windows
    [2022-03-02 14:09:08,325] - pyflex - INFO: Single phase group rejection retained 3 windows
    [2022-03-02 14:09:08,326] - pyflex - INFO: Removing duplicates retains 2 windows.
    [2022-03-02 14:09:08,326] - pyflex - INFO: Rejection based on minimum window length retained 2 windows.
    [2022-03-02 14:09:08,326] - pyflex - INFO: SN amplitude ratio window rejection retained 2 windows
    [2022-03-02 14:09:08,333] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:09:08,334] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.


::


    ---------------------------------------------------------------------------

    _RemoteTraceback                          Traceback (most recent call last)

    _RemoteTraceback: 
    """
    Traceback (most recent call last):
      File "/home/bchow/miniconda3/envs/docs/lib/python3.7/concurrent/futures/process.py", line 239, in _process_worker
        r = call_item.fn(*call_item.args, **call_item.kwargs)
      File "/home/bchow/miniconda3/envs/docs/lib/python3.7/concurrent/futures/process.py", line 198, in _process_chunk
        return [fn(*args) for args in chunk]
      File "/home/bchow/miniconda3/envs/docs/lib/python3.7/concurrent/futures/process.py", line 198, in <listcomp>
        return [fn(*args) for args in chunk]
      File "/home/bchow/REPOSITORIES/pyatoa/pyatoa/core/pyaflowa.py", line 420, in _process_event_multiprocess_true
        return self.process_event(*args, **kwargs, multiprocess=True)
      File "/home/bchow/REPOSITORIES/pyatoa/pyatoa/core/pyaflowa.py", line 407, in process_event
        io=io, **kwargs)
      File "/home/bchow/REPOSITORIES/pyatoa/pyatoa/core/pyaflowa.py", line 625, in process_station
        mgmt.plot(corners=self.map_corners, show=False, save=save)
      File "/home/bchow/REPOSITORIES/pyatoa/pyatoa/core/manager.py", line 1078, in plot
        mp.plot(corners=corners, show=show, save=save, **kwargs)
      File "/home/bchow/REPOSITORIES/pyatoa/pyatoa/visuals/mgmt_plot.py", line 149, in plot
        **kwargs)
      File "/home/bchow/REPOSITORIES/pyatoa/pyatoa/visuals/wave_maker.py", line 533, in plot
        syn = self.st_syn.select(component=comp)[0]
      File "/home/bchow/miniconda3/envs/docs/lib/python3.7/site-packages/obspy/core/stream.py", line 649, in __getitem__
        return self.traces.__getitem__(index)
    IndexError: list index out of range
    """

    
    The above exception was the direct cause of the following exception:


    IndexError                                Traceback (most recent call last)

    /tmp/ipykernel_78179/2476287388.py in <module>
          1 # Now, all we need to do to invoke multiprocessing is input a list of source names to the run function
    ----> 2 pf.multi_event_process(source_names=["2018p130600", "2012p242656"], cha="HH?", loc="*")
    

    ~/REPOSITORIES/pyatoa/pyatoa/core/pyaflowa.py in multi_event_process(self, source_names, max_workers, **kwargs)
        441                     source_names,
        442                     executor.map(self._process_event_multiprocess_true,
    --> 443                                  source_names)
        444                     ):
        445                 misfits[os.path.basename(source_name)] = misfit


    ~/miniconda3/envs/docs/lib/python3.7/concurrent/futures/process.py in _chain_from_iterable_of_lists(iterable)
        481     careful not to keep references to yielded objects.
        482     """
    --> 483     for element in iterable:
        484         element.reverse()
        485         while element:


    ~/miniconda3/envs/docs/lib/python3.7/concurrent/futures/_base.py in result_iterator()
        596                     # Careful not to keep a reference to the popped future
        597                     if timeout is None:
    --> 598                         yield fs.pop().result()
        599                     else:
        600                         yield fs.pop().result(end_time - time.monotonic())


    ~/miniconda3/envs/docs/lib/python3.7/concurrent/futures/_base.py in result(self, timeout)
        426                 raise CancelledError()
        427             elif self._state == FINISHED:
    --> 428                 return self.__get_result()
        429 
        430             self._condition.wait(timeout)


    ~/miniconda3/envs/docs/lib/python3.7/concurrent/futures/_base.py in __get_result(self)
        382     def __get_result(self):
        383         if self._exception:
    --> 384             raise self._exception
        385         else:
        386             return self._result


    IndexError: list index out of range


Unfortunately I don’t think it’s possible to run parallel tasks in a
Jupyter notebook and I’m not prepared to put in the effort to figure out
how to do this, so I hope you can take my word for it that this simply
executes two parallel processing steps using the concurrent.futures
multiprocessing machinery. neat! \__\_

Automating earthquake-based FWT with Pyaflowa and SeisFlows3
------------------------------------------------------------

Earthquake-based full waveform tomography (FWT) is a complicated
procedure involving performing data-synthetic comparisons for large
numbers of source-receiver pairs. This is ideally done in parallel as
this processing step is (mostly) identical and independent for each
source-receiver pair.

There are many available workflow tools for automating FWT, but here we
show the combination of ``Pyatoa`` and
`SeisFlows3 <https://github.com/bch0w/seisflows3>`__, a Python-based
workflow tool for automating full waveform tomography, adjoint
tomography, or full waveform inversion (FWI) on a variety of compute
systems.

``Pyaflowa`` has already been implemented as a `SeisFlows3 preprocess
module <https://github.com/bch0w/seisflows3/blob/master/seisflows3/preprocess/pyatoa.py>`__.
You can see Pyaflowa being called in the function
`seisflows3.preprocess.prepare_eval_grad() <https://github.com/bch0w/seisflows3/blob/master/seisflows3/preprocess/pyatoa.py#L162>`__.
You can see that all of the complexity of Pyatoa is contained within
itself, such that SeisFlows3 only needs to make essentially two calls to
unleash all of the power of Pyatoa.

In order to interact with ``Pyaflowa`` during a SeisFlows3 inversion,
the user edits the `SeisFlows3 parameters.yaml
file <https://github.com/bch0w/pyatoa/blob/master/pyatoa/tests/test_data/test_seisflows_parameters.yaml>`__.
Since we cannot run SeisFlows3 within this jupyter notebook (yet), we
will simply invoke some of the SeisFlows specific calls in Pyaflowa to
investigate what is going on inside the black box.

Let’s start by implement Pyaflowa in the same fashion that the
`SeisFlows3 preprocess module
does <https://github.com/bch0w/seisflows3/blob/master/seisflows3/preprocess/pyatoa.py#L180>`__.
The following cell defines two functions that are slightly modified from
SeisFlows3 to allow us to interact with some SeisFlows3 objects without
having to invoke the entire package

.. code:: ipython3

    import re, yaml
    
    class Dict(object):
        """
        At runtime, SeisFlows defines its parameters as Dict objects, 
        Re-defined dictionary-like object for holding parameters or paths
        Allows for easier access of dictionary items, does not allow resets of
        attributes once defined, only allows updates through new dictionaries.
        """
        def __iter__(self):
            return iter(sorted(self.__dict__.keys()))
    
        def __getattr__(self, key):
            return self.__dict__[key]
    
        def __getitem__(self, key):
            return self.__dict__[key]
    
        def __setattr__(self, key, val):
            if key in self.__dict__:
                raise TypeError("Once defined, parameters cannot be changed.")
            self.__dict__[key] = val
    
        def __delattr__(self, key):
            if key in self.__dict__:
                raise TypeError("Once defined, parameters cannot be deleted.")
            raise KeyError
    
        def update(self, newdict):
            super(Dict, self).__setattr__('__dict__', newdict)
    
        def __init__(self, newdict):
            self.update(newdict)
        
        def __str__(self):
            """
            Pretty print dictionaries and first level nested dictionaries
            """
            str_ = "{"
            for key, item in vars(self).items():
                if isinstance(item, str):
                    str_ += f"{key}: '{item}'\n"
                elif isinstance(item, dict):
                    str_ += key + ": {\n"
                    for key2, item2 in vars(self)[key].items():
                        if isinstance(item2, str):
                            str_ += f"\t{key2}: '{item2}'\n"
                        else:
                            str_ += f"\t{key2}: {item2}\n"
                    str_ += "}\n"
                else:
                    str_ += f"{key}: {item}\n"
            str_ += "}"
            return str_
        
        
    def loadyaml(filename):
        """
        Define how the PyYaml yaml loading function behaves. 
        Replaces None and inf strings with NoneType and numpy.inf respectively
        
        :type filename: str
        :param filename: .yaml file to load in
        """
        # work around PyYAML bugs
        yaml.SafeLoader.add_implicit_resolver(
            u'tag:yaml.org,2002:float',
            re.compile(u'''^(?:
             [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
            |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
            |\\.[0-9_]+(?:[eE][-+][0-9]+)?
            |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
            |[-+]?\\.(?:inf|Inf|INF)
            |\\.(?:nan|NaN|NAN))$''', re.X),
            list(u'-+0123456789.'))
    
        with open(filename, 'r') as f:
            mydict = yaml.safe_load(f)
    
        if mydict is None:
            mydict = dict()
    
        # Replace 'None' and 'inf' values to match expectations
        for key, val in mydict.items():
            if val == "None":
                mydict[key] = None
            if val == "inf":
                mydict[key] = np.inf
    
        return mydict

.. code:: ipython3

    # Here we instantiate the SeisFlows3 PAR and PATH variables, which are
    # dictionaries that contain all of the necessary user-defined information
    # to run an inversion
    PAR = loadyaml("../tests/test_data/test_seisflows_parameters.yaml")
    PATH = PAR["PATHS"]
    PAR.pop("PATHS")


::


    ---------------------------------------------------------------------------

    FileNotFoundError                         Traceback (most recent call last)

    /tmp/ipykernel_78179/1015256936.py in <module>
          2 # dictionaries that contain all of the necessary user-defined information
          3 # to run an inversion
    ----> 4 PAR = loadyaml("../tests/test_data/test_seisflows_parameters.yaml")
          5 PATH = PAR["PATHS"]
          6 PAR.pop("PATHS")


    /tmp/ipykernel_78179/2408363588.py in loadyaml(filename)
         75         list(u'-+0123456789.'))
         76 
    ---> 77     with open(filename, 'r') as f:
         78         mydict = yaml.safe_load(f)
         79 


    FileNotFoundError: [Errno 2] No such file or directory: '../tests/test_data/test_seisflows_parameters.yaml'


.. code:: ipython3

    # Print out the PAR and PATH variables to see what they define
    for p in [PAR, PATH]:
        for key, val in p.items():
            print(f"{key:<15}: {val}")
        print("\n")


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    /tmp/ipykernel_78179/2443140216.py in <module>
          1 # Print out the PAR and PATH variables to see what they define
    ----> 2 for p in [PAR, PATH]:
          3     for key, val in p.items():
          4         print(f"{key:<15}: {val}")
          5     print("\n")


    NameError: name 'PAR' is not defined


.. code:: ipython3

    pyaflowa = Pyaflowa(structure="seisflows", sfpaths=PATH, 
                        sfpar=PAR, iteration=1, step_count=0)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    /tmp/ipykernel_78179/811192365.py in <module>
    ----> 1 pyaflowa = Pyaflowa(structure="seisflows", sfpaths=PATH, 
          2                     sfpar=PAR, iteration=1, step_count=0)


    NameError: name 'PATH' is not defined


