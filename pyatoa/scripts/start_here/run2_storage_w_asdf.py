"""
Example Script 2: Data/Metadata storage w/ ASDF (Adaptable Seismic Data Format)

Same functionality as Example Script 1, but results are stored in ASDFDataSets
which are HDF5 files designed for seismic data

.. rubric:
    python 2_storage_w_asdf.py
"""
import os
from obspy import read, read_events, read_inventory
from pyasdf import ASDFDataSet
from pyatoa import Manager, Config, logger

logger.setLevel("DEBUG")

test_data = "../../tests/test_data/"

# Create an ASDFDataSet to be filled by the Manager
ds = ASDFDataSet("out2_example_asdf_dataset.h5", mode="w")

# Generate ObsPy classes of data and dataless
st_obs = read(os.path.join(test_data, "test_obs_data_NZ_BFZ_2018p130600.ascii"))
st_syn = read(os.path.join(test_data, "test_syn_data_NZ_BFZ_2018p130600.ascii"))
event = read_events(os.path.join(test_data, "test_catalog_2018p130600.xml"))[0]
inv = read_inventory(os.path.join(test_data, "test_dataless_NZ_BFZ.xml"))

# Create the Manager using read in data
cfg = Config(event_id="2018p130600")
mgmt = Manager(ds=ds, config=cfg, event=event, st_obs=st_obs, st_syn=st_syn, 
               inv=inv)
mgmt.write()
mgmt.flow()

# Pretty print parts of the manager for user to view
border = f"\n{'='*80}\n"
def pprint(header, value):
    print(f"{border}{header}{border}")
    print(value)

pprint("MANAGER > mgmt", mgmt)
pprint("DATASET > ds", ds)
pprint("WAVEFORMS > ds.waveforms.NZ_BFZ", ds.waveforms.NZ_BFZ)
pprint("OBSERVED > ds.waveforms.NZ_BFZ.observed", ds.waveforms.NZ_BFZ.observed)
pprint("SYNTHETIC > ds.waveforms.NZ_BFZ.synthetic", 
       ds.waveforms.NZ_BFZ.synthetic)
pprint("RESULTS > ds.auxiliary_data", ds.auxiliary_data)
pprint("WINDOWS > ds.auxiliary_data.MisfitWindows.default", 
       ds.auxiliary_data.MisfitWindows.default)
pprint("EXAMPLE WINDOW > ds.auxiliary_data.MisfitWindows.default.NZ_BFZ_E_0",
       ds.auxiliary_data.MisfitWindows.default.NZ_BFZ_E_0)
pprint("ADJSRCS > ds.auxiliary_data.AdjointSources.default", 
       ds.auxiliary_data.AdjointSources.default)
pprint("EXAMPLE ADJSRC > ds.auxiliary_data.AdjointSources.default.NZ_BFZ_BXE", 
       ds.auxiliary_data.AdjointSources.default.NZ_BFZ_BXE)

