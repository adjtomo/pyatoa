Tutorial
========
A Pyatoa crash course.

Jumping right in, we'll need to import ``Pyatoa`` and relevant packages

::
    
    import obspy
    import pyatoa
    import pyasdf


Config
______
The adjoint tomography process requires some problem-dependent user-defined choices. These are encapsulated in a ``config`` object. 
To set the configuration, we call the ``Config`` class, with a few parameters.

::

    config = pyatoa.Config(event_id="2018p130600",
                           model_number="m00",
                           rotate_to_rtz=False,
                           unit_output="DISP"
                           )
    print(config)

::

    CONFIG
        Model Number:          m00
        Event ID:              2018p130600
        Minimum Filter Period: 10.0
        Maximum Filter Period: 30.0
        Filter Corners:        4.0
        Rotate to RTZ:         False
        Unit Output:           DISP
        Pyflex Config:         [0.08, 15.0, 1.0, 0.8, 0.7, 4.0, 0.0, 1.0, 2.0, 3.0, 10.0]
        Adjoint Source Type:   multitaper_misfit
        Paths to waveforms:    None
        Paths to synthetics:   None
        Paths to responses:    None

We can see Pyatoa set some default parameters.

ASDF Dataset
____________
``PyASDF`` is an HDF5 wrapper that Pyatoa uses to take care of input/output. It is used to avoid complicated directory structures. 
Let's create a PyASDF dataset in our working directory, and write our config so that we know what parameters we've set.

::

    ds = pyasdf.ASDFDataSet("./example_dataset.h5")
    config.write_to_asdf(ds=ds)
    ds

::
    
    ASDF file [format version: 1.0.2]: 'TEST.h5' (142.4 KB)
        Contains 1 event(s)
        Contains waveform data from 0 station(s).


Logging
_______
Logging allows Pyatoa to tell you what it's doing under the hood, great for prototyping and debugging. 

::
    
    import logging
    logger = logging.getLogger("pyatoa")
    logger.setLevel(logging.DEBUG)

Manager
_______
The main player in ``Pyatoa`` is the ``Manager``, the  workflow manager which controls all other supplementary classes. 
Upon calling the ``Manager`` we tell it where the Config and ASDF dataset are.

::

    mgmt = pyatoa.Manager(config=config, ds=ds)

::

    [2019-01-23 16:06:53,536] - pyatoa - INFO: initiating gatherer
    [2019-01-23 16:06:53,536] - pyatoa - INFO: gathering event information
    [2019-01-23 16:06:54,392] - pyatoa - INFO: geonet moment tensor found for event: 2018p130600
    [2019-01-23 16:06:54,393] - pyatoa - INFO: appending GeoNet moment tensor information to event
    [2019-01-23 16:06:54,393] - pyatoa - INFO: event got from external
    
If given an event ID, the Manager will automatically fetch moment tensor information and store it for future use.

::
     
    print(mgmt)

::

    CRATE
        Event:                     smi:nz.org.geonet/2018p130600
        Inventory:                 False
        Observed Stream(s):        False
        Synthetic Stream(s):       False
    MANAGER
        Obs Data Preprocessed:       False
        Syn Data Preprocessed:       False
        Synthetic Data Shifted:      False
        Pyflex runned:               False
        Pyadjoint runned:            False

We can see that the Manager holds all the data we need for the tomography problem in a ``Crate``. 
The Manager itself has a workflow checklist so that the User can keep track of what steps have been completed.
The Manager expects only one event and one station at a time.

Data Gathering
______________
    
Great, now we can gather some data. We give the Manager a station code and let it search for data.

::

    # Station codes must be in the form NN.SSSS.LL.CCC 
    # (N=network, S=station, L=location, C=channel). 
    # Wildcards okay
    mgmt.gather_data(station_code="NZ.PUZ.*.HH?")
    
Because we did not specify any data paths for Pyatoa to search,
it will automatically use ObsPy to search for data via FDSN webservices 


::

    [2019-01-24 12:43:28,201] - pyatoa - INFO: GATHERING NZ.PUZ.*.HH? for 2018p130600
    [2019-01-24 12:43:28,201] - pyatoa - INFO: gathering station information
    [2019-01-24 12:43:28,201] - pyatoa - INFO: internal station information not found, searching ext.
    [2019-01-24 12:43:28,492] - pyatoa - INFO: gathering observation waveforms
    [2019-01-24 12:43:28,492] - pyatoa - INFO: internal observation data unavailable, searching ext.
    [2019-01-24 12:43:31,620] - pyatoa - INFO: stream got external NZ.PUZ.*.HH?
    [2019-01-24 12:43:31,620] - pyatoa - INFO: gathering synthetic waveforms
    No synthetic waveforms for NZ.PUZ.*.HH? found for given event

::
    
    print(mgmt)
    CRATE
        Event:                     smi:nz.org.geonet/2018p130600
        Inventory:                 NZ.PUZ
        Observed Stream(s):        3
        Synthetic Stream(s):       False
    MANAGER
        Obs Data Preprocessed:       False
        Syn Data Preprocessed:       False
        Synthetic Data Shifted:      False
        Pyflex runned:               False
        Pyadjoint runned:            False

The Manager has successfully retrieved our 3-component observation waveforms (trimmed to the event time),
as well as the response information, via FDSN webservices.
Synthetics are user provided from their solver of choice, and their location specified in the Config.
Since we didn't specify, Pyatoa returns no synthetic streams.


The waveforms are accessible as ObsPy stream objects, which gives the User the flexibility to manipulate them outside the workflow.

::

    print(mgmt.st_obs)
    3 Trace(s) in Stream:
    NZ.PUZ.10.HHE | 2018-02-18T07:43:28.128389Z - 2018-02-18T07:52:08.128389Z | 100.0 Hz, 52001 samples
    NZ.PUZ.10.HHN | 2018-02-18T07:43:28.128389Z - 2018-02-18T07:52:08.128389Z | 100.0 Hz, 52001 samples
    NZ.PUZ.10.HHZ | 2018-02-18T07:43:28.128389Z - 2018-02-18T07:52:08.128389Z | 100.0 Hz, 52001 samples


For the sake of this example, we will give the ``Config`` a local path to synthetic data
so we can continue along the workflow.

::

    config.paths["synthetics"] = ["/Users/chowbr/Documents/subduction/seismic"]
    print(mgmt.config)

::

    CONFIG
        Model Number:          m00
        Event ID:              2018p130600
        Minimum Filter Period: 10.0
        Maximum Filter Period: 30.0
        Filter Corners:        4.0
        Rotate to RTZ:         False
        Unit Output:           DISP
        Pyflex Config:         [0.08, 15.0, 1.0, 0.8, 0.7, 4.0, 0.0, 1.0, 2.0, 3.0, 10.0]
        Adjoint Source Type:   multitaper_misfit
        Paths to waveforms:    []
        Paths to synthetics:   ['/Users/chowbr/Documents/subduction/seismic']
        Paths to responses:    []
