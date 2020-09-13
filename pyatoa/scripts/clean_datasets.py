"""
Utility function to delete auxiliary data from all datasets in the work dir.
"""
from glob import glob
from pyasdf import ASDFDataSet as asdf
from pyatoa.utils.asdf.clean import clean_dataset

for fid in glob("*.h5"):
    with asdf(fid) as ds:
        clean_dataset(ds)
        print(fid)
