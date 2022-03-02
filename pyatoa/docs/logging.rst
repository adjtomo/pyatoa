Logging
=======

Pyatoa comes with a detailed ``logger``, which has varying levels of
output information. For simple statements denoting the status of the
workflow, use the ``INFO`` setting. For detailed logging which details
each step of a Pyatoa workflow, including output values, use the
``DEBUG`` setting. By default, the logger is set to ``WARNING`` only,
which only outputs information when something unexpected occurs.

To instantiate a logger, just import and set. Pyflex and Pyadjoint
loggers can be set in the same manner.

.. code:: ipython3

    import obspy
    from pyatoa import logger, Manager, Config
    from pyflex import logger as pflogger
    from pyadjoint import logger as palogger
    logger.setLevel("DEBUG")
    pflogger.setLevel("DEBUG")
    palogger.setLevel("DEBUG")
    
    # Read in test data
    inv = obspy.read_inventory("../tests/test_data/test_dataless_NZ_BFZ.xml")
    cat = obspy.read_events("../tests/test_data/test_catalog_2018p130600.xml")
    event = cat[0]
    st_obs = obspy.read("../tests/test_data/test_obs_data_NZ_BFZ_2018p130600.ascii")
    st_syn = obspy.read("../tests/test_data/test_syn_data_NZ_BFZ_2018p130600.ascii")
    
    mgmt = Manager(config=Config(), inv=inv, event=event, st_obs=st_obs, st_syn=st_syn)
    mgmt.flow()


.. parsed-literal::

    [2022-03-02 14:08:43] - pyatoa - DEBUG: Component list set to E/N/Z
    [2022-03-02 14:08:43] - pyatoa - INFO: standardizing streams
    [2022-03-02 14:08:43] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHE (0, 0) samples
    [2022-03-02 14:08:43] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHE: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:08:43] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHN (0, 0) samples
    [2022-03-02 14:08:43] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHN: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:08:43] - pyatoa - DEBUG: zero pad NZ.BFZ.10.HHZ (0, 0) samples
    [2022-03-02 14:08:43] - pyatoa - DEBUG: new starttime NZ.BFZ.10.HHZ: 2018-02-18T07:43:28.127644Z
    [2022-03-02 14:08:43] - pyatoa - DEBUG: time offset is -20.0s
    [2022-03-02 14:08:43] - pyatoa - INFO: preprocessing observation data
    [2022-03-02 14:08:43] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:43] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-02 14:08:43] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-02 14:08:43] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:43] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-02 14:08:43] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-02 14:08:43] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-02 14:08:43] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-02 14:08:43] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-02 14:08:43] - pyatoa - INFO: running Pyflex w/ map: default
    [2022-03-02 14:08:43,662] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:43,663] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:43,663] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:43,664] - pyflex - INFO: Initial window selection yielded 4 possible windows.
    [2022-03-02 14:08:43,664] - pyflex - INFO: Rejection based on travel times retained 4 windows.
    [2022-03-02 14:08:43,664] - pyflex - INFO: Rejection based on minimum window length retained 4 windows.
    [2022-03-02 14:08:43,665] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:43,665] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:43,665] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:43,665] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:43,666] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:43,666] - pyflex - INFO: Rejection based on data fit criteria retained 3 windows.
    [2022-03-02 14:08:43,667] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:43] - pyatoa - INFO: 1 window(s) selected for comp E
    [2022-03-02 14:08:43,762] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:43,762] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:43,763] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:43,764] - pyflex - INFO: Initial window selection yielded 4 possible windows.
    [2022-03-02 14:08:43,764] - pyflex - INFO: Rejection based on travel times retained 4 windows.
    [2022-03-02 14:08:43,764] - pyflex - INFO: Rejection based on minimum window length retained 4 windows.
    [2022-03-02 14:08:43,764] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:43,764] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:43,765] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:43,765] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:43,765] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:43,766] - pyflex - DEBUG: Window rejected due to amplitude fit: -1.394198
    [2022-03-02 14:08:43,766] - pyflex - INFO: Rejection based on data fit criteria retained 2 windows.
    [2022-03-02 14:08:43,766] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:43] - pyatoa - INFO: 1 window(s) selected for comp N
    [2022-03-02 14:08:43,862] - pyflex - INFO: Calculated travel times.
    [2022-03-02 14:08:43,863] - pyflex - INFO: Calculating envelope of synthetics.
    [2022-03-02 14:08:43,863] - pyflex - INFO: Calculating STA/LTA.
    [2022-03-02 14:08:43,864] - pyflex - INFO: Initial window selection yielded 4 possible windows.
    [2022-03-02 14:08:43,864] - pyflex - INFO: Rejection based on travel times retained 4 windows.
    [2022-03-02 14:08:43,864] - pyflex - INFO: Rejection based on minimum window length retained 4 windows.
    [2022-03-02 14:08:43,865] - pyflex - INFO: Water level rejection retained 4 windows
    [2022-03-02 14:08:43,865] - pyflex - INFO: Single phase group rejection retained 4 windows
    [2022-03-02 14:08:43,865] - pyflex - INFO: Removing duplicates retains 3 windows.
    [2022-03-02 14:08:43,865] - pyflex - INFO: Rejection based on minimum window length retained 3 windows.
    [2022-03-02 14:08:43,865] - pyflex - INFO: SN amplitude ratio window rejection retained 3 windows
    [2022-03-02 14:08:43,866] - pyflex - DEBUG: Window rejected due to amplitude fit: -1.356086
    [2022-03-02 14:08:43,866] - pyflex - DEBUG: Window rejected due to amplitude fit: -1.367611
    [2022-03-02 14:08:43,867] - pyflex - INFO: Rejection based on data fit criteria retained 1 windows.
    [2022-03-02 14:08:43,867] - pyflex - INFO: Weighted interval schedule optimization retained 1 windows.
    [2022-03-02 14:08:43] - pyatoa - INFO: 1 window(s) selected for comp Z
    [2022-03-02 14:08:43] - pyatoa - WARNING: Manager has no ASDFDataSet, cannot save windows
    [2022-03-02 14:08:43] - pyatoa - INFO: 3 window(s) total found
    [2022-03-02 14:08:43] - pyatoa - DEBUG: running Pyadjoint w/ type: cc_traveltime_misfit
    [2022-03-02 14:08:43] - pyatoa - INFO: 0.365 misfit for comp E
    [2022-03-02 14:08:43] - pyatoa - INFO: 1.620 misfit for comp N
    [2022-03-02 14:08:43] - pyatoa - INFO: 0.004 misfit for comp Z
    [2022-03-02 14:08:43] - pyatoa - WARNING: Manager has no ASDFDataSet, cannot save adjoint sources
    [2022-03-02 14:08:43] - pyatoa - INFO: total misfit 1.989


