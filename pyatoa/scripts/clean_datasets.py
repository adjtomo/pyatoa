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
print("="*79)
with ASDFDataSet(h5_files[0]) as ds:
    for iter_ in ds.auxiliary_data.Configs.list():
        steps_ = ds.auxiliary_data.Configs[iter_].list()
        print(f"{iter_}: {steps_}")
print("="*79)

# Choose which parts to delete based on user-input
iterations = input("iterations (csv): ")
if not iterations:
    iterations = None
else:
    iterations = iterations.split(",")

step_counts = input("step_counts (csv): ")
if not step_counts:
    step_counts = None
else:
    step_counts = step_counts.split(",")

# Make sure that the user is okay with the decision
check = input(f"Deleting data from iteration ({iterations}) and "
              f"step_count ({step_counts}) for {len(h5_files)} files. "
              f"Ok? y/[n]: ")
if check != "y":
    sys.exit(0)

for i, fid in enumerate(h5_files):
    print(f"{i+1:0>2}/{len(h5_files):0>2}: {fid}")
    with ASDFDataSet(fid) as ds:
        # Simply remove all synthetic related data
        if iterations is None and step_counts is None:
            clean_ds(ds, iteration=None, step_count=None)

        # Remove specific iterations and/or step counts
        elif iterations is not None:
            for iteration in iterations:
                if step_counts is not None:
                    for step_count in step_counts:
                        clean_dataset(ds, iteration=iteration, 
                                      step_count=step_count)
                else:
                    clean_dataset(ds, iteration=iteration, step_count=None)

print("finished. thank you.")
