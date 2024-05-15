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
        process.apply_filter(st_obs, min_period=min_period,
                             max_period=max_period),
        process.apply_filter(st_obs, min_period=min_period,
                             min_freq=1/max_period),
        process.apply_filter(st_obs, max_freq=1/min_period,
                             min_freq=1/max_period),
        process.apply_filter(st_obs, max_freq=1 / min_period,
                             max_period=max_period)
    ]
    check_result = outputs[0][0].data
    for output in outputs:
        assert np.array_equal(output[0].data, check_result)


def test_zero_pad(st_obs):
    """
    Ensure that zero padding adds the same number of data points each time
    """
    dt = st_obs[0].stats.delta  # time step
    npts = st_obs[0].stats.npts
    pad_length_in_seconds = 20
    pl_samples = pad_length_in_seconds / dt

    st_pad_before = process.zero_pad(st_obs.copy(), pad_length_in_seconds,
                                     before=True, after=False)
    st_pad_after = process.zero_pad(st_obs.copy(), pad_length_in_seconds,
                                    before=False, after=True)
    st_pad_both = process.zero_pad(st_obs.copy(), pad_length_in_seconds,
                                   before=True, after=True)

    for st, npts_check in zip(
            [st_pad_before, st_pad_after, st_pad_both],
            [pl_samples, pl_samples, pl_samples * 2]
    ):
        assert(st[0].stats.npts - npts == npts_check)


def test_trim_streams(st_obs, st_syn):
    """
    Test Stream trimming functionality
    """
    st_a = st_obs.copy()
    st_b = st_syn.copy()

    # Resample first so that we can actually trim the exact same start and 
    # end times
    st_a.resample(st_b[0].stats.sampling_rate)

    # Trim A so we know it is different from B
    st_a.trim(starttime=st_b[0].stats.starttime + 100,
              endtime=st_b[0].stats.endtime - 100)

    for tr_a, tr_b in zip(st_a, st_b):
        assert(tr_a.stats.starttime != tr_b.stats.starttime)
        assert(tr_a.stats.endtime != tr_b.stats.endtime)

    # Trim streams to the same time window
    st_at, st_bt = process.trim_streams(st_a=st_a, st_b=st_b, force="a")

    for tr_a, tr_b in zip(st_at, st_bt):
        assert(tr_a.stats.starttime == tr_b.stats.starttime)
        assert(tr_a.stats.endtime == tr_b.stats.endtime)

    # Make sure this throws an error if it is unable to trim streams
    with pytest.raises(AssertionError):
        process.trim_streams(st_a=st_a, st_b=st_b, force="b")


def test_match_npts(st_obs, st_syn):
    """
    Ensure that forcing number of points standardization works
    """
    st_a = st_obs.copy()
    st_b = st_syn.copy()

    # Resample A to B and trim A to B
    st_a.resample(st_b[0].stats.sampling_rate)

    # Trim A so we know it is different from B
    st_a.trim(starttime=st_b[0].stats.starttime + 100,
              endtime=st_b[0].stats.endtime - 100)
    assert(st_a[0].stats.npts != st_b[0].stats.npts)

    # Trim streams so that they are the same length
    st_a, st_b = process.trim_streams(st_a=st_a, st_b=st_b, force="a")

    # Force one stream to be slightly shorter
    st_a.trim(endtime=st_a[0].stats.endtime - 1)
    
    # Check they're different
    for tr_a, tr_b in zip(st_a, st_b):
        assert(tr_a.stats.npts != tr_b.stats.npts)

    # Match npts
    st_at, st_bt = process.match_npts(st_a, st_b, force="b")
    assert(st_at[0].stats.npts == st_bt[0].stats.npts)


def test_is_preprocessed(st_obs):
    """
    Test the check function that determines if a stream is preprocessed
    """
    # Check for a few different, commonly used processing functions
    st = st_obs.copy()
    assert(process.is_preprocessed(st) is False)
    st.filter("bandpass", freqmin=1, freqmax=5)
    assert(process.is_preprocessed(st) is True)

    st = st_obs.copy()
    assert(process.is_preprocessed(st) is False)
    st.detrend("demean")
    assert(process.is_preprocessed(st, filter_only=False) is True)
    assert(process.is_preprocessed(st, filter_only=True) is False)

    st = st_obs.copy()
    assert(process.is_preprocessed(st) is False)
    st.resample(st[0].stats.sampling_rate // 2 )
    assert(process.is_preprocessed(st, filter_only=False) is True)


def test_stf_convolve():
    """
    !!! TO DO
    """
    pass

