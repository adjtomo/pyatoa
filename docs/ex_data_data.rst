Data-Data Misfit
================

The following code snippet downloads data for a given earthquake from two
nearby stations and generates misfit windows and adjoint sources for each.

The example event used is an `M6.0 in northern Alaska
<https://ds.iris.edu/ds/nodes/dmc/tools/event/10934221>`__. We will use two
stations in interior Alaska, `IM.IL31 and TA.K27K
<https://ds.iris.edu/wilber3/find_stations/10934221>`__.

.. note::

    In this example we show data-data misfit due to the difficulty of generating
    synthetic data for a small example problem.

.. code:: python

    from obspy import UTCDateTime
    from obspy.clients.fdsn import Client
    from pyatoa.utils.srcrcv import merge_inventories

    # First we download event metadata from IRIS
    c = Client("IRIS")
    event = c.get_events(eventid="10934221")[0]
    origintime = event.preferred_origin().time

    # Then we grab station metadata for IM.IL31 and TA.K27K
    inv = merge_inventories(
            inv_a=c.get_stations(network="IM", station="IL31", location="*",
                                 channel="BH?", starttime=origintime,
                                 endtime=origintime + 600, level="response")
    inv

    # Then we grab waveforms and station metadata
    st_1 = c.get_waveforms(network="IM", station="IL31", location="*",
                           channel="BH?", starttime=origintime,
                           endtime=origintime + 600)
    st_2 = c.get_waveforms(network="TA", station="K27K", location="*",
                           channel="BH?", starttime=origintime,
                           endtime=origintime + 600)



