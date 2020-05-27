"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import os
import json
import logging
import numpy as np
from IPython import embed
from pyasdf import ASDFDataSet
from pyatoa import Config, Manager, logger
from pyatoa.core.manager import ManagerError
from pyatoa.core.gatherer import GathererNoDataException
from obspy import read, read_events, read_inventory


logger.setLevel("DEBUG")


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
    return Config(event_id="2018p130600", client="GEONET")


@pytest.fixture
def window():
    """
    Pre-gathered window for this specific source-receiver configuration
    """
    return json.load(open(
        "./test_data/test_window_NZ_BFZ_N_0_2018p130600.json"))["windows"][0]

@pytest.fixture
def adjoint_source():
    """
    Pre-gathered adjoint source for specific configuration to be compared with 
    calculated adjoint source
    """
    return np.loadtxt(
        "./test_data/test_adjoint_source_NZ_BFZ_N_2018p130600.adj", unpack=True
        )

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
    mgmt_pre.preprocess()
    mgmt_pre.window()
    mgmt_pre.measure()

    return mgmt_pre


def test_property_st_shows_correct_number_streams(st_obs, st_syn):
    """
    Ensure that the stream property correctly evaluates the internal streams
    """
    mgmt = Manager(st_obs=st_obs)
    assert(len(mgmt.st) == 3)
    mgmt.st_syn = st_syn
    assert(len(mgmt.st) == 6)


def test_internal_flag_checks(mgmt_pre):
    """
    Ensure assertions in Manager _check catch incorrectly specified parameters
    """
    # embed(colors="neutral")
    assert(mgmt_pre.stats.event_id == "smi:nz.org.geonet/2018p130600")
    assert(mgmt_pre.stats.len_obs == mgmt_pre.stats.len_syn == 3)
    assert(mgmt_pre.stats.standardized == False)
    assert(mgmt_pre.stats.inv_name == "NZ.BFZ")

def test_gather_pass_exception():
    """
    Gathering functionality is tested in test_gather, but test the gather() 
    function by raising exceptions and seeing how Pyatoa passes then
    """
    pass


def test_read_write_from_asdfdataset(tmpdir, mgmt_pre, config):
    """
    Write a Manager into an ASDFDataSet and then read it back
    """
    with ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5")) as ds:
        mgmt_pre.ds = ds
        mgmt_pre.write(write_to="ds")

        # Load data back from dataset
        mgmt_loaded = Manager(ds=ds, config=config)
        mgmt_loaded.load("NZ.BFZ")

        # Manager has no equivalent represenation so just check some ids
        assert(mgmt_pre.stats.event_id == mgmt_loaded.stats.event_id)
        assert(mgmt_pre.stats.len_obs == mgmt_loaded.stats.len_obs)
        assert(mgmt_pre.stats.len_syn == mgmt_loaded.stats.len_syn)
        assert(mgmt_pre.stats.inv_name == mgmt_loaded.stats.inv_name)


def test_standardize_to_synthetics(mgmt_pre):
    """
    Ensure that standardizing streams provides the correct functionality.
    """
    # Values to check standardization against
    syn_starttime = mgmt_pre.st_syn[0].stats.starttime
    syn_samp_rate = mgmt_pre.st_syn[0].stats.sampling_rate
    syn_npts = mgmt_pre.st_syn[0].stats.npts

    # Ensure these don't match the observed values
    with pytest.raises(AssertionError):
        for tr in mgmt_pre.st_obs:
            assert(tr.stats.starttime == syn_starttime)
            assert(tr.stats.sampling_rate == syn_samp_rate)
            assert(tr.stats.npts == syn_npts)

    mgmt_pre.standardize()
    # Ensure stat is set correctly
    assert(mgmt_pre.stats.standardized == True)

    # Ensure that these values now match the observed values
    for tr in mgmt_pre.st_obs:
        assert(tr.stats.starttime == syn_starttime)
        assert(tr.stats.sampling_rate == syn_samp_rate)
        assert(tr.stats.npts == syn_npts)


def test_standardize_raises_manager_error(mgmt_pre):
    """
    Manager will raise manager error if no traces present
    """
    # Make sure Manager won't standardize if waveforms don't match
    mgmt_pre.st_syn = None
    mgmt_pre.stats.len_syn = 0
    with pytest.raises(ManagerError):
        mgmt_pre.standardize()

def test_preprocess_dont_rotate(mgmt_pre):
    """
    Standard preprocessing, dont rotate components, just filter, check if 
    filtering worked
    """
    mgmt_pre.standardize().preprocess()
    assert(mgmt_pre.stats.obs_filtered)
    assert(mgmt_pre.stats.syn_filtered)


def test_preprocess_rotate_to_rtz(mgmt_pre):
    """
    Standard preprocessing but rotate components based on the backazimuth
    """
    mgmt_pre.config.rotate_to_rtz = True
    mgmt_pre.standardize().preprocess()
    # Make sure component names were changed
    for tr in mgmt_pre.st:
        assert(tr.id[-1] in ["R", "T", "Z"])
    # Make sure the backazimuth calculated by ObsPy matches whats expected for 
    # this source-receiver configuration, only to two significant digits
    assert(float(f"{mgmt_pre.baz:.2f}") == 3.21)


# def test_convolve_source_time_function(mgmt_pre):
#     """
#     Ensure that convolution with the source time function works by checking the 
#     data is changed when running stf_convolve

#     TO DO:
#         de-convolve the STF signal and test against original data?
#     """
#     st_obs_data = {tr.id: tr.data for tr in mgmt_pre.st_obs}
#     st_syn_data = {tr.id: tr.data for tr in mgmt_pre.st_syn}

#     mgmt_pre._convolve_source_time_function()

#     # Ensure the synthetic data has been convolved
#     with pytest.raises(AssertionError):
#         for tr in mgmt_pre.st_syn:
#             assert((st_syn_data[tr.id] != tr.data).all())
#     # Ensure observed data remains the same
#     for tr in mgmt_pre.st_obs:
#         assert((st_obs_data[tr.id] == tr.data).all())


# def test_convolve_source_time_function_synthetic_synthetic_case(mgmt_pre):
#     """
#     Ensure that the synthetic-synthetic case convolves both streams
#     """
#     st_data = {tr.id: tr.data for tr in mgmt_pre.st}

#     mgmt_pre.config.synthetics_only = True
#     mgmt_pre._convolve_source_time_function()

#     # Ensure the synthetic data has been convolved
#     with pytest.raises(AssertionError):
#         for tr in mgmt_pre.st:
#             assert((tr.data == st_data[tr.id]).all())


def test_window(mgmt_pre, window):
    """
    Ensure windows functionality works as advertised
    """
    assert(mgmt_pre.config.pyflex_preset == "default")

    # Check that error is raised if insufficient workflow progress
    with pytest.raises(ManagerError):
        mgmt_pre.window()
    mgmt_pre.standardize().preprocess().window()

    # Ensure the correct number of windows are chosen
    for comp, nwin in {"N": 1, "E": 1}.items():
        assert(len(mgmt_pre.windows[comp]) == nwin)
    with pytest.raises(KeyError):
        mgmt_pre.windows["Z"]

    # Assert that the gathered window matches the test data window that was
    # saved from Pyflex. Names will be different
    win_check = mgmt_pre.windows["N"][0]

    assert(win_check.left == window["left_index"])
    assert(win_check.right == window["right_index"])
    assert(win_check.max_cc_value == window["max_cc_value"])
    assert(win_check.cc_shift == window["cc_shift_in_samples"])
    assert(win_check.dlnA == window["dlnA"])

def test_save_windows_and_use_fixed_windows_from_dataset(tmpdir, mgmt_post):
    """
    Test fixed windowing by adding a set of windows to a temp dataset and then
    retrieving them and comparing. Also test save_windows() at the same time
    """
    with ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5")) as ds:
        mgmt_post.ds = ds
        # Explicitely set the model and step count
        mgmt_post.config.model = 0
        mgmt_post.config.step = 0
        mgmt_post.config.save_to_ds = True
        saved_windows = mgmt_post.windows
        mgmt_post.save_windows()  # saved to path 'm00/s00'

        # Delete windows, iterate step, retrieve fixed windows
        mgmt_post.windows = None
        mgmt_post.config.step += 1
        mgmt_post.window(fix_windows=True)

        # Just check some parameter for each window to make sure all goods
        for comp in mgmt_post.windows:
            for w, window in enumerate(mgmt_post.windows[comp]):
                for attr in ["left", "right", "cc_shift"]:
                    assert(getattr(window, attr) == 
                                getattr(saved_windows[comp][w], attr))

def test_save_adjsrcs(tmpdir, mgmt_post, adjoint_source):
    """
    Checks that adjoint sources can be written to dataset and will match the 
    formatting required by Specfem3D
    """
    t_check, data_check = adjoint_source
    with ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5")) as ds:
        mgmt_post.ds = ds
        mgmt_post.save_adjsrcs()

        assert(hasattr(ds.auxiliary_data.AdjointSources.default, "NZ_BFZ_BXN"))
        adj_src = ds.auxiliary_data.AdjointSources.default.NZ_BFZ_BXN.data
        t, data = adj_src.value.T

        # Ensure the data and time arrays are as expected
        assert((t == t_check).all())
        assert((data == data_check).all())


def test_get_path_for_aux_data(mgmt_pre):
    """
    Ensure that path naming works as advertise
    """
    assert(mgmt_pre._get_path_for_aux_data() == "default")        
    mgmt_pre.config.model = "m00"
    assert(mgmt_pre._get_path_for_aux_data() == "m00")  
    mgmt_pre.config.step = "s00"
    assert(mgmt_pre._get_path_for_aux_data() == "m00/s00")  


def test_format_windows():
    """
    """
    pass











