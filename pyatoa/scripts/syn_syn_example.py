"""
Pyatoa Example: Synthetic-Synthetic misfit quantification

Generate adjoint sources for pairs of synthetic data generated from two 
different models using SPECFEM 2D/3D/3D_GLOBE.

.. note::

    It is expected that the User has already generated their synthetic data
    using SPECFEM and that the same STATIONS list was used 

https://pyatoa.readthedocs.io/syn-syn.html
"""
import os
from obspy import UTCDateTime
from pysep.utils.io import read_sem
from pyatoa import Config, Manager


# Determine where the synthetics are stored
path_syn = "./SGF/AAS000000"
path_obs = "./EGF/AAS000000"
path_out = "./SEM"

if not os.path.exists(path_out):
    os.mkdirs(path_out)

# Now we initiate Pyatoa Config object which controls processing
pyflex_cfg = {"stalta_waterlevel": 0.05, 
              "tshift_acceptance_level": 25,
              "cc_acceptance_level": 0.6
              }
cfg = Config(min_period=10, max_period=200, component_list=["Y"],
             pyflex_preset="default", adj_src_type="multitaper",
             st_obs_type="syn", st_syn_type="syn", **pyflex_cfg
             )
# Used for determining the time offset T0 of the synthetics
dummy_origintime = UTCDateTime("2000-01-01T00:00:00")

# Loop through all synthetics and process
for fid_syn, fid_obs in zip(sorted(os.listdir(path_syn)), 
                            sorted(os.listdir(path_obs))):

    st_syn = read_sem(os.path.join(path_syn, fid_syn), 
                      origintime=dummy_origintime) 
    st_obs = read_sem(os.path.join(path_obs, fid_obs), 
                      origintime=dummy_origintime) 

    # Provide Manager with Config and data, let it do the rest
    mgmt = Manager(config=cfg, st_obs=st_obs, st_syn=st_syn)
    mgmt.flow()
    mgmt.plot(show=False, save=os.path.join(path_out, "figures",
                                            f"{st_obs[0].get_id()}.png"))

    # Output adjoint sources
    mgmt.stats.time_offset_sec = st_obs[0].stats.starttime - dummy_origintime
    mgmt.write_adjsrcs(path=path_out, write_blanks=True)


