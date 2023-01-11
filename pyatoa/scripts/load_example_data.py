"""
It's useful to generate a fully loaded Manager object for testing purposes. 
Load data from the test directory and run the Manager workflow to achieve this.
"""
import os

# Assuming directory structure within the package is static
_root_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
_test_data_dir = os.path.join(_root_dir, "..", "tests", "test_data")


def load_example_data():
    """
    Returns example test data that can be used to preload a Manager
    """
    from pyatoa import Config
    from obspy import read, read_events, read_inventory

    cfg = Config(event_id="2018p130600", min_period=8., max_period=20.,
                 iteration=1, step_count=0)
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
    """Returns an example Inspector loaded with some waveform measurements"""
    # Avoids circular import
    from pyatoa import Inspector

    insp = Inspector(tag="test_inspector")
    insp.read(path=_test_data_dir)

    return insp


def load_example_asdfdataset():
    """Returns an example ASDFDataSet with a few waveforms etc"""
    from pyasdf import ASDFDataSet

    # Read only format so we don't mess up the example data
    ds = ASDFDataSet(os.path.join(_test_data_dir, "test_ASDFDataSet.h5"),
                     mode="r")

    return ds


def generate_example_asdfdataset():
    """Create the test_ASDFDataSet file, which sometimes needs to be re-made if 
    the package or dependencies change"""
    from pyasdf import ASDFDataSet
    from obspy import Catalog
    from pyatoa import Manager

    cfg, st_obs, st_syn, event, inv = load_example_data()

    with ASDFDataSet("test_ASDFDataSet.h5") as ds:
        cfg.write(write_to=ds)
        ds.events = Catalog(events=[event])
        mgmt = Manager(config=cfg, ds=ds, st_obs=st_obs, st_syn=st_syn, 
                       event=event, inv=inv)
        mgmt.write(ds=ds)
        mgmt.flow()

