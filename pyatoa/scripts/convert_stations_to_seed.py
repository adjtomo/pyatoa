"""
Short utility script to convert a SPECFEM STATIONS file to a collection of 
per-channel StationXML files required by the Pyatoa automatic data gatherer.

See also pyatoa.utils.form.convert_stations_to_seed()
"""
import os
from pyatoa.utils.read import read_stations
from pyatoa.utils.write import write_inv_seed

stations_file = "./STATIONS"
path_to_save_seed = "./seed"

# Create the output directory
if not os.path.exists(path_to_save_seed):
    os.mkdir(path_to_save_seed)

# Read SPECFEM Stations file
inv = read_stations(stations_file)

# Write inventory as a collection of StationXML files
write_inv_seed(inv=inv, path=path_to_save_seed)
