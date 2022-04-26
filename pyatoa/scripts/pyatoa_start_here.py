"""
The simplest working example for testing out the misfit quantification 
capabilities of Pyatoa. In order, this script:

1. Reads in waveforms (observed and synthetic), metadata, earthquake metadata. 
2. Preprocesses, standardizes, windows and measures the two waveforms
3. Produces a waveform misfit figure, and a source-receiver map

.. rubric:
    python pyatoa_start_here.py
"""
from pyatoa import Manager, Config, logger
from obspy import read, read_events, read_inventory

logger.setLevel("DEBUG")

# Generate ObsPy classes of data and dataless
st_obs = read("../tests/test_data/test_obs_data_NZ_BFZ_2018p130600.ascii")
st_syn = read("../tests/test_data/test_syn_data_NZ_BFZ_2018p130600.ascii")
event = read_events("../tests/test_data/test_catalog_2018p130600.xml")[0]
inv = read_inventory("../tests/test_data/test_dataless_NZ_BFZ.xml")

# Create the Manager using read in data
cfg = Config(event_id="2018p130600")
mgmt = Manager(config=cfg, event=event, st_obs=st_obs, st_syn=st_syn, inv=inv)
mgmt.flow()
mgmt.plot(save="output_figure_pyatoa_start_here.png")

