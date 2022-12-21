"""
Pyatoa Example: Data-Data misfit quantification

Downloads data for an earthquake in Alaska and quantifies misfit, generating
an adjoint source and misfit windows. See docs page for explanation:

https://pyatoa.readthedocs.io/data-data.html
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from pyatoa import Config, Manager
from pyatoa.utils.srcrcv import merge_inventories

# First we download event metadata from IRIS
c = Client("IRIS")
event = c.get_events(eventid="10934221")[0]
origintime = event.preferred_origin().time

# Then we grab station metadata for IM.IL31 and TA.K27K. Use a Pyatoa utility
# to combine the two inventories
inv = merge_inventories(
        inv_a=c.get_stations(network="TA", station="POKR", location="",
                             channel="BH?", starttime=origintime,
                             endtime=origintime + 300, level="response"),
        inv_b=c.get_stations(network="TA", station="J25K", location="",
                             channel="BH?", starttime=origintime,
                             endtime=origintime + 300, level="response")
)

# Then we grab waveforms and station metadata
st_1 = c.get_waveforms(network="TA", station="POKR", location="",
                       channel="BH?", starttime=origintime,
                       endtime=origintime + 300)
st_2 = c.get_waveforms(network="TA", station="J25K", location="",
                       channel="BH?", starttime=origintime,
                       endtime=origintime + 300)

# Now we initiate Pyatoa Config object which controls processing
pyflex_cfg = {"stalta_waterlevel": 0.05, 
              "tshift_acceptance_level": 25,
              "cc_acceptance_level": 0.6
              }
cfg = Config(min_period=10, max_period=30,
             unit_output="DISP", rotate_to_rtz=True,
             pyflex_preset="default", adj_src_type="cc_traveltime_misfit",
             st_obs_type="obs", st_syn_type="obs", **pyflex_cfg
             )

# Provide Manager with Config and data, let it do the rest
mgmt = Manager(config=cfg, st_obs=st_1, st_syn=st_2, 
               inv=inv, event=event
               )
mgmt.write()
mgmt.flow()
mgmt.plot()

# Output adjoint sources
mgmt.write_adjsrcs(path="./", write_blanks=True)


