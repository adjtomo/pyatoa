"""
Utility function to delete auxiliary data from all datasets in the work dir.
Allows for user defined parameters to control how much data to delete
"""
import sys
from glob import glob
from pyasdf import ASDFDataSet
from pyatoa.utils.asdf.clean import clean_dataset

h5_files = glob("*.h5")
if not len(h5_files):
    sys.exit("No HDF5 files found in directory")

# Check available iterations and step_count, assuming all files are the same
with ASDFDataSet(h5_files[0]) as ds:
    for iter_ in ds.auxiliary_data.Configs.list():
        steps_ = ds.auxiliary_data.Configs[iter_].list()
        print(f"{iter_}: {steps_}")

# Choose which parts to delete based on user-input
iteration = input("iteration: ")
if not iteration:
    iteration = None
step_count = input("step_count: ")
if not step_count:
    step_count = None

# Make sure that the user is okay with the decision
check = input(f"Deleting data from iteration ({iteration}) and "
              f"step_count ({step_count}) for {len(h5_files)} files. "
              f"Ok? y/[n]: ")
if check != "y":
    sys.exit(0)

for i, fid in enumerate(h5_files):
    print(f"{i:0>2}/{len(h5_files):0>2}: {fid}")
    with ASDFDataSet(fid) as ds:
        clean_dataset(ds, iteration=iteration, step_count=step_count)

print("finished. thank you.")
