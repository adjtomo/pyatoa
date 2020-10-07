"""
Test the functionalities of the suite of processing functions in the utilities
"""
import pytest
import numpy as np
from pyatoa.utils import process
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
    """
    Default Pyatoa Config object
    """
    return Config(event_id="2018p130600", client="GEONET")


@pytest.fixture
def mgmt_pre(config, event, st_obs, st_syn, inv):
    """
    A manager filled with data but pre-workflow
    """
    return Manager(config=config, event=event, st_obs=st_obs, st_syn=st_syn,
                   inv=inv)


def test_filters(st_obs):
    """
    Make sure st_obs works with various input frequencies/ periods
    """
    # Only want to check one component really
    st_obs = st_obs.select(component="Z")

    min_period = 10
    max_period = 30
    outputs = [
        process.filters(st_obs, min_period=min_period, max_period=max_period),
        process.filters(st_obs, min_period=min_period, min_freq=1/max_period),
        process.filters(st_obs, max_freq=1/min_period, min_freq=1/max_period),
        process.filters(st_obs, max_freq=1 / min_period, max_period=max_period)
    ]
    check_result = outputs[0][0].data
    for output in outputs:
        assert np.array_equal(output[0].data, check_result)


def test_trim_streams(st_obs, st_syn):
    """
    Make sure that time offset tapering does what it's expected to
    :param st_obs:
    :return:
    """
    st_a = st_obs.copy()
    st_b = st_syn.copy()

    starttime_a = st_a[0].stats.starttime
    endtime_a = st_a[0].stats.endtime
    starttime_b = st_b[0].stats.starttime
    endtime_b = st_b[0].stats.endtime

    # Just make sure these values are different
    assert starttime_a != starttime_b
    assert endtime_a != endtime_b

    st_a, st_b = process.trim_streams(st_obs, st_syn, force="a")
    assert (st_a[0].stats.starttime != starttime_b)
    assert (st_a[0].stats.endtime != endtime_b)
    assert (st_b[0].stats.starttime == starttime_a)
    # assert (st_b[0].stats.endtime == endtime_a)  # endtime_b < endtime_a


    # This doesn't work because the starttimes differ by milliseconds?
    st_a, st_b = process.trim_streams(st_obs, st_syn, force="b")
    import ipdb;ipdb.set_trace()
    assert (st_a[0].stats.starttime == starttime_b)
    assert (st_a[0].stats.endtime == endtime_b)
    assert (st_b[0].stats.starttime != starttime_a)
    assert (st_b[0].stats.endtime != endtime_a)


def test_match_npts(st_obs, st_syn):
    """

    :param st_obs:
    :param st_syn:
    :return:
    """