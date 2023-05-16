"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import os
import json
import numpy as np
from pyasdf import ASDFDataSet
from pyatoa import Config, Manager, logger
from pyatoa.core.manager import ManagerError
from obspy import read, read_events, read_inventory


# Turn off the logger for tests
logger.propogate = False
logger.setLevel("CRITICAL")


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


# @pytest.fixture
# def window():
#     """
#     Pre-gathered window for this specific source-receiver configuration
#     This doesn't actually match the data... strange. Not used.
#     """
#     return json.load(open(
#         "./test_data/test_window_NZ_BFZ_N_0_2018p130600.json"))["windows"][0]
#
#
# @pytest.fixture
# def adjoint_source():
#     """
#     Pre-gathered adjoint source for specific configuration to be compared with
#     calculated adjoint source
#     """
#     return np.loadtxt(
#         "./test_data/test_adjoint_source_NZ_BFZ_N_2018p130600.adj", unpack=True
#         )


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


def test_read_write_from_asdfdataset(tmpdir, mgmt_pre, config):
    """
    Write a Manager into an ASDFDataSet and then read it back
    """
    with ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5")) as ds:
        mgmt_pre.ds = ds
        mgmt_pre.write()

        # Load data back from dataset
        mgmt_loaded = Manager(ds=ds, config=config)
        mgmt_loaded.load("NZ.BFZ", path="default")

        # Manager has no equivalence represenation so just check some ids
        assert(mgmt_pre.stats.event_id == mgmt_loaded.stats.event_id)
        assert(mgmt_pre.stats.len_obs == mgmt_loaded.stats.len_obs)
        assert(mgmt_pre.stats.len_syn == mgmt_loaded.stats.len_syn)
        assert(mgmt_pre.stats.inv_name == mgmt_loaded.stats.inv_name)


def test_standardize_to_synthetics(mgmt_pre):
    """
    Ensure that standardizing streams performs three main tasks, trimming
    origin times, matching sampling rates, and matching number of points.
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

    # Ensure that these values now match the observed values
    assert(mgmt_pre.stats.standardized == True)
    assert(mgmt_pre.stats.time_offset_sec == -20.)

    for tr in mgmt_pre.st_obs:
        assert(tr.stats.starttime == syn_starttime)
        assert(tr.stats.sampling_rate == syn_samp_rate)
        assert(tr.stats.npts == syn_npts)


def test_standardize_raises_manager_error(mgmt_pre):
    """
    Asser that Manager will raise an error if user tries to standardize with
    no traces present
    """
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
    assert mgmt_pre.stats.obs_processed
    assert mgmt_pre.stats.syn_processed


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


def test_preprocess_overwrite(mgmt_pre):
    """
    Apply an overwriting preprocessing function to ensure functionality works
    """
    # First ensure that only functions can be passed
    with pytest.raises(AssertionError):
        mgmt_pre.preprocess(overwrite="not a function")

    def preproc_fx(st, choice, value=1, **kwargs):
        """Zero out the data for an easy check on results"""
        for tr in st:
            tr.data *= value
        return st

    mgmt_pre.preprocess(overwrite=preproc_fx, value=0)

    for tr in mgmt_pre.st:
        assert(not tr.data.any())


def test_select_window(mgmt_pre):
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


def test_save_and_retrieve_windows(tmpdir, mgmt_post):
    """
    Test retrieve_windows() and save_windows() by saving windows into a
    scratch dataset and retrieving them back. Window criteria will be
    recalculated but since the waveforms are the same, the values will be the
    same as before.
    """
    with ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5")) as ds:
        mgmt_post.ds = ds
        # Explicitely set the model and step count
        mgmt_post.config.iteration = 0
        mgmt_post.config.step_count = 0
        mgmt_post.config.save_to_ds = True
        saved_windows = mgmt_post.windows
        mgmt_post.save_windows()  # saved to path 'm00/s00'

        # Delete windows, iterate step, retrieve fixed windows
        mgmt_post.windows = None
        mgmt_post.config.step_count += 1
        mgmt_post.window(fix_windows=True)

        # Just check some parameter for each window to make sure all goods
        for comp in mgmt_post.windows:
            for w, window in enumerate(mgmt_post.windows[comp]):
                for attr in ["left", "right", "cc_shift"]:
                    assert(getattr(window, attr) ==
                           getattr(saved_windows[comp][w], attr))

        # Delete windows, slightly change synthetic waveforms and check to make
        # sure that recalculated criteria are different
        mgmt_post.windows = None
        for tr in mgmt_post.st_syn:
            tr.data *= 2
        mgmt_post.window(fix_windows=True)

        for comp in mgmt_post.windows:
            for w, window in enumerate(mgmt_post.windows[comp]):
                # Amplitude ratios will be different since we multipled them
                assert (getattr(window, "dlnA") !=
                        getattr(saved_windows[comp][w], "dlnA"))


def test_save_adjsrcs(tmpdir, mgmt_post):
    """
    Checks that adjoint sources can be written to dataset and will match the 
    formatting required by Specfem3D
    """
    with ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5")) as ds:
        mgmt_post.ds = ds
        mgmt_post.save_adjsrcs()
        assert(hasattr(ds.auxiliary_data.AdjointSources.default, "NZ_BFZ_BXN"))


def test_format_windows(mgmt_post):
    """
    Basic check that format windows returns as formatted lists expected
    """
    adjoint_windows = mgmt_post._format_windows()
    assert(isinstance(adjoint_windows, dict))
    for key, window in adjoint_windows.items():
        assert(isinstance(key, str))
        assert(isinstance(window, list))
        for values in window:
            assert(isinstance(values, list))
            assert(len(values) == 2)
            for value in values:
                assert(isinstance(value, float))

def test_flow_multiband(mgmt_pre):
    """
    Test that the workflow for multiple period bands returns a single
    adjoint source
    """
    windows, adjsrcs = mgmt_pre.flow_multiband(
        periods=[(1, 10), (10, 30), (15, 40)]
    )
    # Just check that the expected values don't change
    assert(pytest.approx(adjsrcs["E"].max(), .001) == 8914.48)
    assert(pytest.approx(adjsrcs["N"].max(), .001) == 3173.05)
    assert(pytest.approx(adjsrcs["Z"].max(), .001) == 2749.96)


def test_resample_numerical_noise(mgmt_pre):
    """
    Raised in Issue #34 by Ridvan O. (rdno).
    
    Numerical noise can be introduced by resampling a trace that already has 
    the correct sampling rate, which leads to non-zero adjoint sources/misfit 
    when traces are identical due to very slight time shifts introduced from 
    the resampling method.

    This check ensures the fix (do not resample if sampling rates are the same)
    continues to work.
    """
    # Ensure that both traces are the same 
    mgmt_pre.st_obs = mgmt_pre.st_syn.copy()

    # Take some parameters from the Issue
    mgmt_pre.config.adj_src_type = "waveform"
    mgmt_pre.config.min_period = 20
    mgmt_pre.config.max_period = 100.
    mgmt_pre.config.st_obs_type = "syn"
    mgmt_pre.config._set_external_configs()

    mgmt_pre.standardize()
    mgmt_pre.preprocess()
    mgmt_pre.measure()  # skip windowing as we get the same result without

    # mgmt_pre.plot(choice="wav")  # creates the figure shown in Issue #34 

    # adjoint sources should be zero because these are the same trace
    for comp, adjsrc in mgmt_pre.adjsrcs.items():
        assert(adjsrc.adjoint_source.max() == 0)
        assert(adjsrc.adjoint_source.min() == 0)

  

