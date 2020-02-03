"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import pyasdf
from pyatoa import Config, Manager
from obspy import read, read_events, read_inventory


@pytest.fixture
def obs():
    """Real data stream"""
    return read("./test_data/test_obs_data _NZ_BFZ_2018p130600.ascii")


@pytest.fixture
def syn():
    """Synthetic data stream"""
    return read("./test_data/test_syn_data _NZ_BFZ_2018p130600.ascii")


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


def test_empty():
    """Ensure that opening an empty Manager doesnt work"""
    with pytest.raises(TypeError):
        Manager()


def test_init_empty(config):
    """Ensure that opening a Manager doesnt allow workflow"""
    mgmt = Manager(config=config, empty=True)
    assert(mgmt)


# def test_