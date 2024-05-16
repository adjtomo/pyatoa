"""
Test the PyASDF ASDFDataSet auxiliary utilities

!!! TO DO: utils.asdf.write tests
"""
import os
import pytest
from obspy import read, read_events, read_inventory
from pyasdf import ASDFDataSet
from pyatoa import Config, Manager
from pyatoa.utils.asdf import (add, clean, load)


@pytest.fixture
def empty_dataset(tmpdir):
    """Re-used test data pointing to STATIONS file"""
    fid = os.path.join(tmpdir, "empty_dataset.h5")
    # Make sure the dataset is actually empty
    if os.path.exists(fid):
        os.remove(fid)
    return ASDFDataSet(fid)


@pytest.fixture
def dataset():
    """Filled ASDFDataSet"""
    return ASDFDataSet("./test_data/test_ASDFDataSet.h5")


@pytest.fixture
def st_obs():
    """
    Raw observed waveforms from station NZ.BFZ.HH? for New Zealand event
    2018p130600 (GeoNet event id)
    """
    return read("./test_data/test_obs_data_NZ_BFZ_2018p130600.ascii")


@pytest.fixture
def st_syn():
    """
    Synthetic data stream generated using Specfem3D and focal mechanism for
    2018p130600. Minimum resolved period roughly 10s.
    """
    return read("./test_data/test_syn_data_NZ_BFZ_2018p130600.ascii")


@pytest.fixture
def cat():
    """
    ObsPy Event Catalog for New Zealand based event with
    GeoNet Event ID: 2018p130600
    """
    return read_events("./test_data/test_catalog_2018p130600.xml")


@pytest.fixture
def event(cat):
    """
    Event from Catalog
    """
    return cat[0]


@pytest.fixture
def inv():
    """
    StationXML information for station NZ.BFZ.HH?
    """
    return read_inventory("./test_data/test_dataless_NZ_BFZ.xml")


@pytest.fixture
def config():
    """
    Default Pyatoa Config object
    """
    return Config(event_id="2018p130600")


@pytest.fixture
def mgmt_pre(config, event, st_obs, st_syn, inv):
    """
    A manager filled with data but pre-workflow
    """
    return Manager(config=config, event=event, st_obs=st_obs, st_syn=st_syn,
                   inv=inv)


@pytest.fixture
def mgmt_post(mgmt_pre):
    """
    A manager that has completed the full workflow
    """
    mgmt_pre.standardize()
    mgmt_pre.preprocess(remove_response=True, output="DISP")
    mgmt_pre.window()
    mgmt_pre.measure()

    return mgmt_pre


def test_add_misfit_windows(empty_dataset, mgmt_post):
    """
    Test adding misfit windows to an ASDFDataSet
    :return:
    """
    path = mgmt_post.config.aux_path
    add.add_misfit_windows(windows=mgmt_post.windows, ds=empty_dataset,
                           path=path)

    windows = empty_dataset.auxiliary_data.MisfitWindows[path]
    assert(len(windows.list()) == 3)

    for comp, window_list in mgmt_post.windows.items():
        assert(len(window_list) == 1)
        window_check = window_list[0]
        # Dynamically construct station tags to access the individual parameters
        net, sta, loc, cha = window_check.channel_id.split(".")
        tag = f"{net}_{sta}_{comp}_0"
        par = windows[tag].parameters

        # Check a few parameters
        assert(par["left_index"] == window_check.left)
        assert(par["right_index"] == window_check.right)
        assert(par["absolute_starttime"] ==
               str(window_check.absolute_starttime)
               )

        # Retrieve and check phase arrivals
        # !!! This doesn't work because PyFlex is returning multiple entries
        # !!! for the same phase. Need to look into this.
        # for arrival in window_check.phase_arrivals:
        #     assert(par[f"phase_arrival_{arrival['name']}"] == arrival["time"])


def test_add_adjoint_sources(empty_dataset, mgmt_post):
    """
    Test adding adjoint sources to an ASDFDataSet
    :return:
    """
    path = mgmt_post.config.aux_path
    add.add_adjoint_sources(adjsrcs=mgmt_post.adjsrcs, ds=empty_dataset,
                            path=path,
                            time_offset=mgmt_post.stats.time_offset_sec
                            )

    adjsrcs = empty_dataset.auxiliary_data.AdjointSources[path]
    assert(len(adjsrcs.list()) == 3)

    # Just check misfit value is correct
    for comp, adjsrc in mgmt_post.adjsrcs.items():
        tag = f"NZ_BFZ_BX{comp}"
        misfit_check = adjsrcs[tag].parameters["misfit"]
        assert(adjsrc.misfit == misfit_check)


def test_clean_dataset(empty_dataset, mgmt_post):
    """
    Test dataset clean functions. Need to perform tasks on a dataset we create
    here, otherwise we may permanently affect test data if we use a pre-built
    ASDFDataSet
    """
    assert(not hasattr(empty_dataset.auxiliary_data, "MisfitWindows"))

    mgmt_post.write_to_dataset(ds=empty_dataset)
    assert(hasattr(empty_dataset.auxiliary_data, "MisfitWindows"))
    clean.clean_dataset(empty_dataset, fix_windows=False)

    assert(not hasattr(empty_dataset.auxiliary_data, "MisfitWindows"))


def test_clean_dataset_fix_windows(empty_dataset, mgmt_post):
    """
    Test cleaning a dataset but retaining windows
    :return:
    """
    assert(not hasattr(empty_dataset.auxiliary_data, "MisfitWindows"))

    mgmt_post.write_to_dataset(ds=empty_dataset)
    assert(hasattr(empty_dataset.auxiliary_data, "MisfitWindows"))
    clean.clean_dataset(empty_dataset, fix_windows=True)

    assert(hasattr(empty_dataset.auxiliary_data, "MisfitWindows"))
    assert(not hasattr(empty_dataset.auxiliary_data, "AdjointSources"))


def test_load_windows(dataset):
    """
    Test the function that returns windows in the Pyflex output format from
    an ASDFDataSet
    """
    windows = load.load_windows(ds=dataset, net="NZ", sta="BFZ",
                                iteration="i01", step_count="s00")
    ds_windows = dataset.auxiliary_data.MisfitWindows.i01.s00
    check_val = ds_windows.NZ_BFZ_N_0.parameters["max_cc_value"]
    assert(windows["N"][0].max_cc_value == check_val)


def test_load_previous_windows(empty_dataset, mgmt_post):
    """
    Test the function that returns windows in the Pyflex output format from
    an ASDFDataSet
    """
    # Add two sets of windows to the dataset
    for iter_ in ["i01", "i02"]:
        mgmt_post.config.iteration = iter_
        mgmt_post.config.step_count = "s00"
        mgmt_post.write_to_dataset(ds=empty_dataset)

    windows = empty_dataset.auxiliary_data.MisfitWindows
    check = load.previous_windows(windows, iteration="i02", step_count="s00")

    # Access a private attribute to check that we have the correct iter/step
    path = check._AuxiliaryDataAccessor__auxiliary_data_type
    _, iter_, step = path.split("/")
    assert(iter_ == "i01")
    assert(step == "s00")


def test_load_adjsrcs(dataset):
    """
    Test the function that returns windows in the Pyflex output format from
    an ASDFDataSet

    :param dataset:
    :return:
    """
    adjsrcs = load.load_adjsrcs(ds=dataset, net="NZ", sta="BFZ",
                                iteration="i01", step_count="s00")
    ds_adjsrcs = dataset.auxiliary_data.AdjointSources.i01.s00
    check_val = ds_adjsrcs.NZ_BFZ_BXN.parameters["misfit"]
    assert(adjsrcs["N"].misfit == check_val)



