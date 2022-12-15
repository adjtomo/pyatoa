"""
It's useful to generate a fully loaded Manager object for testing purposes. 
Load data from the test directory and run the Manager workflow to achieve this.
"""
import os
from pyatoa import Config, Inspector
from obspy import read, read_events, read_inventory


# Assuming directory structure within the package is static
_root_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
_test_data_dir = os.path.join(_root_dir, "..", "tests", "test_data")


def load_example_data():
    """
    Returns example test data that can be used to preload a Manager
    """
    cfg = Config(event_id="2018p130600")
    st_obs = read(os.path.join(_test_data_dir, 
                               "test_obs_data_NZ_BFZ_2018p130600.ascii"))
    st_syn = read(os.path.join(_test_data_dir, 
                               "test_syn_data_NZ_BFZ_2018p130600.ascii"))
    event = read_events(os.path.join(_test_data_dir, 
                        "test_catalog_2018p130600.xml"))[0]
    inv = read_inventory(os.path.join(_test_data_dir, 
                                      "test_dataless_NZ_BFZ.xml"))
    
    return cfg, st_obs, st_syn, event, inv


def load_example_inspector():
    """
    Returns example window and misfit information that can be used to preload
    an Inspector
    """
    insp = Inspector(tag="test_inspector")
    insp.read(path=_test_data_dir)

    return insp