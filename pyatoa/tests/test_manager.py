"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import pyasdf
from IPython import embed
from pyatoa import Config, Manager
from obspy import read, read_events, read_inventory


@pytest.fixture
def obs():
    """Real data stream"""
    return read("./test_data/test_obs_data_NZ_BFZ_2018p130600.ascii")


@pytest.fixture
def syn():
    """Synthetic data stream"""
    return read("./test_data/test_syn_data_NZ_BFZ_2018p130600.ascii")


@pytest.fixture
def cat():
    """Event Catalog"""
    return read_events("./test_data/test_event_2018p130600.xml")


@pytest.fixture
def event(cat):
    """Event from Catalog"""
    return cat[0]


@pytest.fixture
def inv():
    """StationXML information"""
    return read_inventory("./test_data/test_inv.xml")


@pytest.fixture
def config():
    """Pyatoa Config object"""
    return Config(model_number="m00", event_id="2018p130600")


@pytest.fixture
def ds():
    """Pyasdf ASDFDataSet"""
    return pyasdf.ASDFDataSet("./test_data/test_empty_dataset.h5")


def test_init_empty(config):
    """Ensure that opening a Manager doesnt allow workflow"""
    mgmt = Manager(config=config, empty=True)
    assert mgmt


def test_standardize_resample(config, obs, syn):
    """Ensure that stream standardization resampling works"""
    mgmt = Manager(config=config, st_obs=obs, st_syn=syn)

    # Assert that streams are not standardized
    for tro, trs in zip(mgmt.st_obs, mgmt.st_syn):
        assert(tro.stats.sampling_rate != trs.stats.sampling_rate)

    mgmt.standardize()
    for tro, trs in zip(mgmt.st_obs, mgmt.st_syn):
        assert(tro.stats.sampling_rate == trs.stats.sampling_rate)
        pytest.set_trace()


def test_standardize_zero_pad(config, obs, syn):
    """Ensure Zero padding of data works"""
    mgmt = Manager(config=config, st_obs=obs, st_syn=syn)
    mgmt.config.zero_pad = 10
    mgmt.standardize()
    pytest.set_trace()

# def test_preprocess():
# def test_preprocess_overwrite():
# def test_window():
# def test_window_fixed():
#
# def test_load():
# def test_write():
