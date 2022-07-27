"""
Example Script 1: Data processing with the Pyatoa Manager

The simplest working example for testing out the misfit quantification 
capabilities of Pyatoa. In order, this script:

1. Reads in waveforms (observed and synthetic), metadata, earthquake metadata. 
2. Preprocesses, standardizes, windows and measures the two waveforms
3. Produces a waveform misfit figure, and a source-receiver map

.. rubric:
    python 1_manager_flow.py
"""
import os
from obspy import read, read_events, read_inventory
from pyatoa import Manager, Config, logger

logger.setLevel("DEBUG")

test_data = "../../tests/test_data/"

# Generate ObsPy classes of data and dataless
st_obs = read(os.path.join(test_data, "test_obs_data_NZ_BFZ_2018p130600.ascii"))
st_syn = read(os.path.join(test_data, "test_syn_data_NZ_BFZ_2018p130600.ascii"))
event = read_events(os.path.join(test_data, "test_catalog_2018p130600.xml"))[0]
inv = read_inventory(os.path.join(test_data, "test_dataless_NZ_BFZ.xml"))

# Create the Manager using read in data
cfg = Config(event_id="2018p130600")
mgmt = Manager(config=cfg, event=event, st_obs=st_obs, st_syn=st_syn, inv=inv)
mgmt.flow()
mgmt.plot(save="out1_output_figure.png")

# Pretty print parts of the manager for user to view
border = f"\n{'='*80}\n"
def pprint(header, value):
    print(f"{border}{header}{border}")
    print(value)

pprint("OBSERVED", mgmt.st_obs)
pprint("SYNTHETIC", mgmt.st_syn)
pprint("INVENTORY", mgmt.inv)
pprint("EVENT", mgmt.event)
pprint("MANAGER", mgmt)
pprint("STATS", mgmt.stats)
