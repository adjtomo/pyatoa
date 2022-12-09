from obspy import read
from pyatoa import Config, Manager, logger

st = read()
mgmt = Manager(st_obs=st.select(component="Z"),
               st_syn=st.select(component="N"))
mgmt.standardize()
mgmt.plot(show=True)
