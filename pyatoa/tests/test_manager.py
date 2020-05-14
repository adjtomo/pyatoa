"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import pyasdf
from IPython import embed
from pyatoa import Config, Manager
from obspy import read, read_events, read_inventory


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
    """Default Pyatoa Config object"""
    return Config(model="m00", step="s00", event_id="2018p130600")


def test_property_st_shows_correct_number_streams(st_obs, st_syn):
    """
    Ensure that the stream property correctly evaluates the internal streams
    """
    mgmt = Manager(st_obs=st_obs)
    assert(len(mgmt.st) == 3)
    mgmt.st_syn = st_syn
    assert(len(mgmt.st) == 6)


def test_internal_flag_checks(event, inv, st_obs, st_syn):
    """
    Ensure assertions in Manager _check catch incorrectly specified parameters
    """
    mgmt = Manager(event=event, st_obs=st_obs, st_syn=st_syn, inv=inv)

    # embed(colors="neutral")
    assert(mgmt._event_resource_id == "smi:nz.org.geonet/2018p130600")
    assert(mgmt._len_obs == 3)
    assert(mgmt._len_syn == 3)
    assert(mgmt._standardize_flag == False)
    assert(mgmt._inv_name == "NZ.BFZ")



