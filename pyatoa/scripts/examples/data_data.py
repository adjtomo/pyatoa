"""
Pyatoa Example: Data-Data misfit quantification

Downloads data for an earthquake in Alaska and quantifies misfit, generating
an adjoint source and misfit windows. See docs page for explanation:

    https://pyatoa.readthedocs.io/data-data.html
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from pyatoa.utils.srcrcv import merge_inventories

# First we download event metadata from IRIS
c = Client("IRIS")
event = c.get_events(eventid="10934221")[0]
origintime = event.preferred_origin().time

# Then we grab station metadata for IM.IL31 and TA.K27K. Use a Pyatoa utility
# to combine the two inventories
inv = merge_inventories(
        inv_a=c.get_stations(network="IM", station="IL31", location="*",
                             channel="BH?", starttime=origintime,
                             endtime=origintime + 600, level="response"),
        inv_b=c.get_stations(network="TA", station="K27K", location="*",
                             channel="BH?", starttime=origintime,
                             endtime=origintime + 600, level="response")
)

# Then we grab waveforms and station metadata
st_1 = c.get_waveforms(network="IM", station="IL31", location="*",
                       channel="BH?", starttime=origintime,
                       endtime=origintime + 600)
st_2 = c.get_waveforms(network="TA", station="K27K", location="*",
                       channel="BH?", starttime=origintime,
                       endtime=origintime + 600)

