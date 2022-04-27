"""
Example Script 3: Bulk assessment with the Pyatoa Inspector

The Inspector combs through ASDFDataSets (created by the Pyatoa Manager) and
aggregates metadata and information on misfit windows and adjoint sources into
Pandas DataFrames. These can be used to explore inversion results, and make
pre-designed figures.

.. rubric:
    python run3_inspector_assessment.py
"""
import os
from pyatoa import Inspector, logger

logger.setLevel("DEBUG")

example_dataset = "./out2_example_asdf_dataset.h5"

assert(os.path.exists(example_dataset)), (
        "Please run 'run2_storage_w_asdf.py' to generate the necessary "
        "example data for this quick start example"
        )

insp = Inspector("out3_example_inspector")
insp.discover()

def pprint(header, value):
    print(f"{border}{header}{border}")
    print(value)

import ipdb;ipdb.set_trace()

