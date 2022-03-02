Data Storage
============

Pyatoa stores data using `PyASDF
ASDFDataSets <https://seismicdata.github.io/pyasdf/asdf_data_set.html>`__,
which are seismological data structures built upon the HDF5 file format.

Datasets are hierarchical (tree-like), portable, compressible, and
self-describing or containing both data and metadata. They are built
around ObsPy objects, removing any need for conversions in the
transition from data storage to data processing.

An ``ASDFDataSet`` can be passed directly to the ``Manager`` class. By
default, gathered data and processed results will automatically be
stored inside the dataset following a pre-defined naming convention.
Naming schemes are set using parameters in the ``Config`` object.

Below we show how data is saved throughout a workflow, and how it can be
accessed using PyASDF and Pyatoa.

For a detailed tutorial on the ``ASDFDataSet``, see:
https://seismicdata.github.io/pyasdf/tutorial.html

.. code:: ipython3

    import os
    import obspy
    from pyatoa import Config, Manager, logger
    from pyasdf import ASDFDataSet
    
    logger.setLevel("DEBUG")
    
    # Load in the test data
    inv = obspy.read_inventory("../tests/test_data/test_dataless_NZ_BFZ.xml")
    cat = obspy.read_events("../tests/test_data/test_catalog_2018p130600.xml")
    event = cat[0]
    st_obs = obspy.read("../tests/test_data/test_obs_data_NZ_BFZ_2018p130600.ascii")
    st_syn = obspy.read("../tests/test_data/test_syn_data_NZ_BFZ_2018p130600.ascii")

--------------

Initializing
------------

| First we must open a new ``ASDFDataSet`` file. We will fill it with
  data from the ``Manager``.
| ``ASDFDataSet``\ s can also be used as a context manager, using the
  ``with`` argument. This ensures the file is closed after use.

.. code:: ipython3

    # Fresh dataset: make sure we aren't trying to work with existing data
    ds_fid = "../tests/test_data/docs_data/test_ASDFDataSet.h5"
    if os.path.exists(ds_fid):
        os.remove(ds_fid)
    
    ds = ASDFDataSet(ds_fid)
    print(ds)


.. parsed-literal::

    ASDF file [format version: 1.0.3]: '../tests/test_data/docs_data/test_ASDFDataSet.h5' (96.0 bytes)
    	Contains 0 event(s)
    	Contains waveform data from 0 station(s).


| We can pass the ``ASDFDataSet`` ds directly to the initialization of
  the ``Manager`` class.
| The string representation of the ``Manager`` class shows us that the
  ``ASDFDataSet`` has been attached, by showing the name of the dataset.

   **NOTE:** In Pyatoa, by convention, each event gets its own
   ``ASDFDataSet``; each ``ASDFDataSet`` should be named using a unique
   event identifier. This ensures that files are kept a reasonable size
   and avoids the need for more complicated internal naming schemes.

.. code:: ipython3

    mgmt = Manager(ds=ds, config=Config(), inv=inv, event=event, st_obs=st_obs, st_syn=st_syn)
    print(mgmt)


.. parsed-literal::

    [2022-03-02 14:09:16] - pyatoa - DEBUG: Component list set to E/N/Z


.. parsed-literal::

    Manager Data
        dataset   [ds]:        test_ASDFDataSet.h5
        quakeml   [event]:     smi:nz.org.geonet/2018p130600
        station   [inv]:       NZ.BFZ
        observed  [st_obs]:    3
        synthetic [st_syn]:    3
    Stats & Status
        half_dur:              0.6989458964552759
        time_offset_sec:       None
        standardized:          False
        obs_processed:         False
        syn_processed:         False
        nwin   [windows]:      None
        misfit [adjsrcs]:      None
    


--------------

Manually writing data
---------------------

We can save the current Manager data using the ``Manager.write()``
function. The Pyatoa ``Config`` object can also be written to the
``ASDFDataSet`` using the ``Config.write()`` function.

Once written, we see the ``ASDFDataSet`` has been populated with event
and station metadata, waveform data, and Config information.

.. code:: ipython3

    mgmt.write()
    mgmt.config.write(write_to=ds)

.. code:: ipython3

    ds




.. parsed-literal::

    ASDF file [format version: 1.0.3]: '../tests/test_data/docs_data/test_ASDFDataSet.h5' (495.4 KB)
    	Contains 1 event(s)
    	Contains waveform data from 1 station(s).
    	Contains 1 type(s) of auxiliary data: Configs



.. code:: ipython3

    ds.events




.. parsed-literal::

    1 Event(s) in Catalog:
    2018-02-18T07:43:48.127644Z | -39.949, +176.300 | 5.156706293 M  | manual



.. code:: ipython3

    ds.waveforms.list()




.. parsed-literal::

    ['NZ.BFZ']



.. code:: ipython3

    ds.auxiliary_data.Configs




.. parsed-literal::

    1 auxiliary data item(s) of type 'Configs' available:
    	default



--------------

Automatically written data
--------------------------

| During a Pyatoa workflow, individual functions will automatically
  write their outputs into the given ``ASDFDataSet``.
| Here the log statements show the ``Manager.window()`` and
  ``Manager.measure()`` functions saving their outputs into the data
  set.

.. code:: ipython3

    mgmt.standardize().preprocess();


.. parsed-literal::

    [2022-03-02 14:09:16] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:09:16] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHE (0, 0) samples
    [2022-03-02 14:09:16] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHE: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:16] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHN (0, 0) samples
    [2022-03-02 14:09:16] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHN: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:16] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHZ (0, 0) samples
    [2022-03-02 14:09:16] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHZ: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:16] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:09:16] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:09:16] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:16] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:09:16] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:09:16] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:16] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:09:16] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:16] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:09:16] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:16] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)


.. code:: ipython3

    mgmt.window();


.. parsed-literal::

    [2022-03-02 14:09:16] - pyatoa - INFO: running Pyflex w/ map: default
    [2022-03-02 14:09:16] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:09:17] - pyatoa - INFO: 1 window(s) selected for comp N
    [2022-03-02 14:09:17] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:09:17] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:09:17] - pyatoa - INFO: 3 window(s) total found


.. code:: ipython3

    mgmt.measure();


.. parsed-literal::

    [2022-03-02 14:09:17] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:09:17] - pyatoa - INFO: 0.365 misfit for comp E
    [2022-03-02 14:09:17] - pyatoa - INFO: 1.620 misfit for comp N
    [2022-03-02 14:09:17] - pyatoa - INFO: 0.004 misfit for comp Z
    [2022-03-02 14:09:17] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:09:17] - pyatoa - INFO: total misfit 1.989


--------------

Accessing saved data using PyASDF
---------------------------------

| All saved data can be accessed using ``ASDFDataSet`` attributes.
| For a more thorough explanation of accessing data with an
  ``ASDFDataSet``, see: https://seismicdata.github.io/pyasdf/index.html

**Event metadata** is stored as an ObsPy ``Catalog`` object in the
``ASDFDataSet.events`` attribute.

.. code:: ipython3

    ds.events[0]




.. parsed-literal::

    Event:	2018-02-18T07:43:48.127644Z | -39.949, +176.300 | 5.156706293 M  | manual
    
    	                  resource_id: ResourceIdentifier(id="smi:nz.org.geonet/2018p130600")
    	                   event_type: 'earthquake'
    	                creation_info: CreationInfo(agency_id='WEL(GNS_Primary)', author='scevent@kseqp01.geonet.org.nz', creation_time=UTCDateTime(2018, 2, 18, 7, 44, 9, 156454))
    	          preferred_origin_id: ResourceIdentifier(id="smi:nz.org.geonet/Origin#20180226021110.13419.62761")
    	       preferred_magnitude_id: ResourceIdentifier(id="smi:nz.org.geonet/Origin#20180226021110.13419.62761#netMag.M")
    	 preferred_focal_mechanism_id: ResourceIdentifier(id="smi:local/ad83e11b-cc91-4de7-9cd0-5c51f99e1062")
    	                         ---------
    	           event_descriptions: 1 Elements
    	             focal_mechanisms: 1 Elements
    	                      origins: 1 Elements
    	                   magnitudes: 3 Elements



--------------

| **Waveforms** are stored as ObsPy ``Stream`` objects, and **station
  metadata** is stored as ObsPy ``Inventory`` objects.
| They are stored together in the ``ASDFDataSet.waveforms`` attribute.

.. code:: ipython3

    ds.waveforms.NZ_BFZ.StationXML




.. parsed-literal::

    Inventory created at 2020-02-02T22:21:59.000000Z
    	Created by: Delta
    		    None
    	Sending institution: GeoNet (WEL(GNS_Test))
    	Contains:
    		Networks (1):
    			NZ
    		Stations (1):
    			NZ.BFZ (Birch Farm)
    		Channels (3):
    			NZ.BFZ.10.HHZ, NZ.BFZ.10.HHN, NZ.BFZ.10.HHE



.. code:: ipython3

    ds.waveforms.NZ_BFZ.observed + ds.waveforms.NZ_BFZ.synthetic




.. parsed-literal::

    6 Trace(s) in Stream:
    NZ.BFZ.10.HHE | 2018-02-18T07:43:28.128394Z - 2018-02-18T07:49:38.128394Z | 100.0 Hz, 37001 samples
    NZ.BFZ.10.HHN | 2018-02-18T07:43:28.128394Z - 2018-02-18T07:49:38.128394Z | 100.0 Hz, 37001 samples
    NZ.BFZ.10.HHZ | 2018-02-18T07:43:28.128394Z - 2018-02-18T07:49:38.128394Z | 100.0 Hz, 37001 samples
    NZ.BFZ..BXE   | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples
    NZ.BFZ..BXN   | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples
    NZ.BFZ..BXZ   | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples



--------------

**Misfit windows**, **Adjoint Sources**, and **Configuration
parameters** are stored in the ``ADSFDataSet.auxiliary_data`` attribute.

.. code:: ipython3

    ds.auxiliary_data




.. parsed-literal::

    Data set contains the following auxiliary data types:
    	AdjointSources (1 item(s))
    	Configs (1 item(s))
    	MisfitWindows (1 item(s))



If no ``iteration`` or ``step_count`` attributes are provided to the
``Config`` object, auxiliary data will be stored using the ``default``
tag.

.. code:: ipython3

    ds.auxiliary_data.MisfitWindows




.. parsed-literal::

    1 auxiliary data sub group(s) of type 'MisfitWindows' available:
    	default



.. code:: ipython3

    ds.auxiliary_data.MisfitWindows['default']




.. parsed-literal::

    3 auxiliary data item(s) of type 'MisfitWindows/default' available:
    	NZ_BFZ_E_0
    	NZ_BFZ_N_0
    	NZ_BFZ_Z_0



.. code:: ipython3

    ds.auxiliary_data.MisfitWindows.default.NZ_BFZ_E_0




.. parsed-literal::

    Auxiliary Data of Type 'MisfitWindows'
    	Path: 'default/NZ_BFZ_E_0'
    	Data shape: '(2,)', dtype: 'int64'
    	Parameters:
    		absolute_endtime: 2018-02-18T07:44:45.197644Z
    		absolute_starttime: 2018-02-18T07:43:42.437644Z
    		cc_shift_in_samples: 36
    		cc_shift_in_seconds: 1.08
    		center_index: 1523
    		channel_id: NZ.BFZ.10.HHE
    		dlnA: -0.70965282411
    		dt: 0.03
    		left_index: 477
    		max_cc_value: 0.871536711295
    		min_period: 10.0
    		phase_arrival_P: 15.2626263355
    		phase_arrival_Pn: 15.1318939626
    		phase_arrival_S: 25.7016469855
    		phase_arrival_Sn: 25.6750945772
    		phase_arrival_p: 14.0460406583
    		phase_arrival_s: 23.6216670031
    		phase_arrival_sP: 18.7800304978
    		relative_endtime: 77.07
    		relative_starttime: 14.31
    		right_index: 2569
    		time_of_first_sample: 2018-02-18T07:43:28.127644Z
    		window_weight: 5.46976440008



.. code:: ipython3

    ds.auxiliary_data.AdjointSources




.. parsed-literal::

    1 auxiliary data sub group(s) of type 'AdjointSources' available:
    	default



.. code:: ipython3

    ds.auxiliary_data.AdjointSources.default




.. parsed-literal::

    3 auxiliary data item(s) of type 'AdjointSources/default' available:
    	NZ_BFZ_BXE
    	NZ_BFZ_BXN
    	NZ_BFZ_BXZ



.. code:: ipython3

    ds.auxiliary_data.AdjointSources.default.NZ_BFZ_BXE




.. parsed-literal::

    Auxiliary Data of Type 'AdjointSources'
    	Path: 'default/NZ_BFZ_BXE'
    	Data shape: '(10000, 2)', dtype: 'float64'
    	Parameters:
    		adj_src_type: cc_traveltime_misfit
    		component: BXE
    		dt: 0.03
    		location: 10
    		max_period: 30.0
    		min_period: 10.0
    		misfit: 0.36539741683
    		network: NZ
    		starttime: 2018-02-18T07:43:28.127644Z
    		station: BFZ



--------------

Re-loading data using the Manager
---------------------------------

Data previously saved into an ``ASDFDataSet`` can be loaded back into a
``Manager`` class using the ``Manager.load()`` function. The ``load()``
function will search for matching metadata, waveforms and configuration
parameters, based on the ``path`` argument provided.

.. code:: ipython3

    mgmt = Manager(ds=ds)
    mgmt.load(code="NZ.BFZ", path="default")


.. parsed-literal::

    [2022-03-02 14:09:17] - pyatoa - INFO: no config provided, initiating default
    [2022-03-02 14:09:17] - pyatoa - DEBUG: Component list set to E/N/Z
    [2022-03-02 14:09:17] - pyatoa - INFO: loading config from dataset default




.. parsed-literal::

    Manager Data
        dataset   [ds]:        test_ASDFDataSet.h5
        quakeml   [event]:     smi:nz.org.geonet/2018p130600
        station   [inv]:       NZ.BFZ
        observed  [st_obs]:    3
        synthetic [st_syn]:    3
    Stats & Status
        half_dur:              0.6989458964552759
        time_offset_sec:       None
        standardized:          False
        obs_processed:         False
        syn_processed:         False
        nwin   [windows]:      None
        misfit [adjsrcs]:      None



.. code:: ipython3

    !pwd


.. parsed-literal::

    /home/bchow/REPOSITORIES/pyatoa/pyatoa/docs


Misfit windows and adjoint sources are not explicitely re-loaded.
Windows can be loaded using optional arguments in the
``Manager.window()`` function.

--------------

Saving data during an inversion
-------------------------------

For each function evaluation, a new set of synthetic waveforms, misfit
windows, adjoint sources and (potentially) configuration parameters, are
defined. Therefore, unique tags are required to save and load this
information in a reliable manner.

Pyatoa tags using the ``Config.iteration`` and ``Config.step_count``
attributes to define unique tags during an inversion.

.. code:: ipython3

    # Set the config iteration and step_count parameters
    cfg = Config(iteration=1, step_count=0)
    
    # Remove the previously created dataset
    os.remove(ds_fid)
    ds = ASDFDataSet(ds_fid)
    
    cfg.write(write_to=ds)
    mgmt = Manager(ds=ds, config=cfg, inv=inv, event=event, st_obs=st_obs, st_syn=st_syn)
    mgmt.write()
    mgmt.flow()


.. parsed-literal::

    [2022-03-02 14:09:17] - pyatoa - DEBUG: Component list set to E/N/Z
    [2022-03-02 14:09:17] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:09:17] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHE (0, 0) samples
    [2022-03-02 14:09:17] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHE: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:17] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHN (0, 0) samples
    [2022-03-02 14:09:17] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHN: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:17] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHZ (0, 0) samples
    [2022-03-02 14:09:17] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHZ: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:17] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:09:17] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:09:17] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:17] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:09:17] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:09:17] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:17] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:09:17] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:17] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:09:17] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:17] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:09:17] - pyatoa - INFO: running Pyflex w/ map: default
    [2022-03-02 14:09:17] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:09:17] - pyatoa - INFO: 1 window(s) selected for comp N
    [2022-03-02 14:09:17] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:09:17] - pyatoa - DEBUG: saving misfit windows to ASDFDataSet
    [2022-03-02 14:09:17] - pyatoa - INFO: 3 window(s) total found
    [2022-03-02 14:09:17] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:09:17] - pyatoa - INFO: 0.365 misfit for comp E
    [2022-03-02 14:09:17] - pyatoa - INFO: 1.620 misfit for comp N
    [2022-03-02 14:09:17] - pyatoa - INFO: 0.004 misfit for comp Z
    [2022-03-02 14:09:17] - pyatoa - DEBUG: saving adjoint sources to ASDFDataSet
    [2022-03-02 14:09:17] - pyatoa - INFO: total misfit 1.989


The ``ASDFDataSet`` is now populated with appropriately tagged data,
denoting which function evaluation it belongs to.

.. code:: ipython3

    ds.waveforms.NZ_BFZ




.. parsed-literal::

    Contents of the data set for station NZ.BFZ:
        - Has a StationXML file
        - 2 Waveform Tag(s):
            observed
            synthetic_i01s00



.. code:: ipython3

    ds




.. parsed-literal::

    ASDF file [format version: 1.0.3]: '../tests/test_data/docs_data/test_ASDFDataSet.h5' (576.9 KB)
    	Contains 1 event(s)
    	Contains waveform data from 1 station(s).
    	Contains 3 type(s) of auxiliary data: AdjointSources, Configs, MisfitWindows



.. code:: ipython3

    ds.waveforms.NZ_BFZ.synthetic_i01s00




.. parsed-literal::

    3 Trace(s) in Stream:
    NZ.BFZ..BXE | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples
    NZ.BFZ..BXN | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples
    NZ.BFZ..BXZ | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples



Auxiliary data will be tagged in a similar fashion, making it simple to
re-access specific function evaluations.

.. code:: ipython3

    ds.auxiliary_data.MisfitWindows




.. parsed-literal::

    1 auxiliary data sub group(s) of type 'MisfitWindows' available:
    	i01



.. code:: ipython3

    ds.auxiliary_data.MisfitWindows.i01




.. parsed-literal::

    1 auxiliary data sub group(s) of type 'MisfitWindows/i01' available:
    	s00



.. code:: ipython3

    ds.auxiliary_data.MisfitWindows.i01.s00




.. parsed-literal::

    3 auxiliary data item(s) of type 'MisfitWindows/i01/s00' available:
    	NZ_BFZ_E_0
    	NZ_BFZ_N_0
    	NZ_BFZ_Z_0



Using the ``Manager.load()`` function, we can specify the unique
``path`` to determine which function evaluation we want to retrieve data
from.

.. code:: ipython3

    mgmt = Manager(ds=ds)
    mgmt.load("NZ.BFZ", path="i01/s00", synthetic_tag="synthetic_i01s00")
    mgmt.standardize().preprocess()


.. parsed-literal::

    [2022-03-02 14:09:17] - pyatoa - INFO: no config provided, initiating default
    [2022-03-02 14:09:17] - pyatoa - DEBUG: Component list set to E/N/Z
    [2022-03-02 14:09:17] - pyatoa - INFO: loading config from dataset i01/s00
    [2022-03-02 14:09:18] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:09:18] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHE (0, 0) samples
    [2022-03-02 14:09:18] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHE: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:18] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHN (0, 0) samples
    [2022-03-02 14:09:18] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHN: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:18] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHZ (0, 0) samples
    [2022-03-02 14:09:18] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHZ: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:09:18] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:09:18] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:09:18] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:18] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:09:18] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:09:18] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:18] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:09:18] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:09:18] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:09:18] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:09:18] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)




.. parsed-literal::

    Manager Data
        dataset   [ds]:        test_ASDFDataSet.h5
        quakeml   [event]:     smi:nz.org.geonet/2018p130600
        station   [inv]:       NZ.BFZ
        observed  [st_obs]:    3
        synthetic [st_syn]:    3
    Stats & Status
        half_dur:              0.6989458964552759
        time_offset_sec:       -20.0
        standardized:          True
        obs_processed:         True
        syn_processed:         True
        nwin   [windows]:      None
        misfit [adjsrcs]:      None



| We can now load in previously retrieved windows from the dataset,
  using the ``Manager.window()`` function.
| Windows misfit criteria will be re-evaluated using the current set of
  data. We can turn off automatic window saving using the optional
  ``save`` argument.

.. code:: ipython3

    mgmt.window(fix_windows=True, iteration=1, step_count=0, save=False)


.. parsed-literal::

    [2022-03-02 14:09:18] - pyatoa - INFO: retrieving windows from dataset
    [2022-03-02 14:09:18] - pyatoa - DEBUG: searching for windows in i01s00
    [2022-03-02 14:09:18] - pyatoa - DEBUG: 3 window(s) found in dataset for NZ.BFZ
    [2022-03-02 14:09:18] - pyatoa - DEBUG: recalculating window criteria
    [2022-03-02 14:09:18] - pyatoa - DEBUG: E0_old - cc:0.87 / dt:36.0 / dlnA:-0.71
    [2022-03-02 14:09:18] - pyatoa - DEBUG: E0_new - cc:0.87 / dt:36.0 / dlnA:-0.71
    [2022-03-02 14:09:18] - pyatoa - DEBUG: N0_old - cc:0.99 / dt:63.0 / dlnA:-0.83
    [2022-03-02 14:09:18] - pyatoa - DEBUG: N0_new - cc:0.99 / dt:63.0 / dlnA:-0.83
    [2022-03-02 14:09:18] - pyatoa - DEBUG: Z0_old - cc:0.99 / dt:0.0 / dlnA:-0.90
    [2022-03-02 14:09:18] - pyatoa - DEBUG: Z0_new - cc:0.99 / dt:0.0 / dlnA:-0.90
    [2022-03-02 14:09:18] - pyatoa - INFO: 3 window(s) total found




.. parsed-literal::

    Manager Data
        dataset   [ds]:        test_ASDFDataSet.h5
        quakeml   [event]:     smi:nz.org.geonet/2018p130600
        station   [inv]:       NZ.BFZ
        observed  [st_obs]:    3
        synthetic [st_syn]:    3
    Stats & Status
        half_dur:              0.6989458964552759
        time_offset_sec:       -20.0
        standardized:          True
        obs_processed:         True
        syn_processed:         True
        nwin   [windows]:      3
        misfit [adjsrcs]:      None



*easy peasy mate*
