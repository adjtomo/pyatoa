"""
It's useful to generate a fully loaded Manager object for testing purposes. 
Load data from the test directory and run the Manager workflow to achieve this.
"""
import ipdb
from IPython import embed
from pyatoa import Manager, Config
from obspy import read, read_events, read_inventory


# Generate ObsPy classes of data and dataless
st_obs = read("../tests/test_data/test_obs_data_NZ_BFZ_2018p130600.ascii")
st_syn = read("../tests/test_data/test_syn_data_NZ_BFZ_2018p130600.ascii")
event = read_events("../tests/test_data/test_catalog_2018p130600.xml")[0]
inv = read_inventory("../tests/test_data/test_dataless_NZ_BFZ.xml")

# Create the Manager using read in data
cfg = Config(event_id="2018p130600")
mgmt = Manager(config=cfg, event=event, st_obs=st_obs, st_syn=st_syn, inv=inv)
mgmt.flow()

ipdb.set_trace()
embed(colors="neutral")
