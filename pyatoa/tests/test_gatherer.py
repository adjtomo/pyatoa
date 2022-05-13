"""
Test the functionalities of the Pyatoa Gatherer class
"""
import pytest
from obspy import read_events
from pyasdf import ASDFDataSet
from pyatoa import Config
from pyatoa.core.gatherer import (Gatherer, append_focal_mechanism_to_event,
                                  get_gcmt_moment_tensors,
                                  get_usgs_moment_tensors)


@pytest.fixture
def code():
    """
    Example NZ station code
    """
    return "NZ.BFZ.??.HH*"


@pytest.fixture
def event_id():
    """
    Example NZ event identifier 
    """
    return "2018p130600"


@pytest.fixture
def dataset_fid():
    """
    The name for the dataset file id used for temporary storage
    :return:
    """
    return "./test_data/test_ASDFDataSet.h5"


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
def origintime(event):
    """
    The origin time of the example event
    :return:
    """
    return event.preferred_origin().time


@pytest.fixture
def config(event_id):
    """
    Default Pyatoa Config object
    """
    return Config(event_id=event_id, client="GEONET", iteration=1,
                  step_count=0, synthetics_only=False, save_to_ds=False)


@pytest.fixture
def gatherer(config, origintime):
    """
    The Gatherer which is responsible for gathering data.
    """
    return Gatherer(config=config, origintime=origintime)


def test_gather_no_client(gatherer, code):
    """
    Make sure that external getter functions always returns None if no Client
    is provided
    """
    gatherer.Client = None

    assert gatherer.get_event_from_fdsn() is None
    assert gatherer.get_inv_from_fdsn(code) is None
    assert gatherer.get_waveform_from_fdsn(code) is None


def test_event_get(gatherer):
    """
    Ensure that querying for an event works. Only tests GeoNet gathering.
    """
    # Get by event id from GeoNet
    assert gatherer.get_event_from_fdsn() is not None

    # Get by origin time from GeoNet
    gatherer.config.event_id = None
    assert gatherer.get_event_from_fdsn() is not None


def test_station_get(gatherer, code):
    """
    Ensure that querying for a station works. Only tests GeoNet gathering.
    """
    net, sta, loc, cha = code.split('.')

    # !!! This throws an ObsPy UserWarning that the StationXML file has version
    # !!! 1, but ObsPy only accepts 1.0 or 1.1. Acceptable warning so pass.
    with pytest.warns(UserWarning):
        inv = gatherer.get_inv_from_fdsn(code, station_level="response")

    assert inv[0].code == net
    assert inv[0][0].code == sta
    assert hasattr(inv[0][0][0], "response")


def test_obs_waveform_get(gatherer, code):
    """
    Ensure that querying for a waveforms works. Only tests GeoNet gathering.
    """
    net, sta, loc, cha = code.split('.')

    st = gatherer.get_waveform_from_fdsn(code)
    assert(len(st) == 3)

    stats = st.select(component="Z")[0].stats
    assert stats.network == net
    assert stats.station == sta


def test_asdf_event_fetch(gatherer, dataset_fid):
    """
    Get event from an ASDFDataSet.
    """
    with ASDFDataSet(dataset_fid) as ds:
        gatherer.ds = ds
        assert gatherer.fetch_event_from_dataset() is not None


def test_asdf_station_fetch(gatherer, dataset_fid, code):
    """
    Get station from an ASDFDataSet
    """
    with ASDFDataSet(dataset_fid) as ds:
        gatherer.ds = ds
        inv = gatherer.fetch_inv_from_dataset(code)
        assert len(inv[0][0]) == 3


def test_asdf_waveform_fetch(gatherer, dataset_fid, code, config):
    """
    Get waveforms from an ASDFDataSet
    """
    with ASDFDataSet(dataset_fid) as ds:
        gatherer.ds = ds
        for tag in [config.observed_tag, config.synthetic_tag]:
            st = gatherer.fetch_waveform_from_dataset(code, tag)
            assert len(st) == 3


def test_fetch_event_by_dir(gatherer, event_id):
    """
    Get event information based on given directory structure. Test the various
    types of input sources that are allowable by Pyatoa
    """
    # No Config path means fetching returns nada
    assert gatherer.fetch_event_by_dir(event_id) is None

    gatherer.config.paths["events"] = "./test_data/"

    prefixes = ["test_CMTSOLUTION_", "test_FORCESOLUTION_", "test_SOURCE_"]
    
    # Test each type of available source from SPECFEM2D and SPECFEM3D
    for prefix in prefixes:
        event = gatherer.fetch_event_by_dir(event_id=event_id, prefix=prefix)

        assert event is not None
        assert event_id in event.resource_id.id


def test_fetch_inv_by_dir(gatherer, code):
    """
    Get response based on given directory structure
    """
    # No Config path means fetching returns nada
    assert gatherer.fetch_inv_by_dir(code) is None

    gatherer.config.paths["responses"] = "./test_data/test_seed"
    inv = gatherer.fetch_inv_by_dir(code)
    assert inv is not None
    assert hasattr(inv[0][0][0], "response")


def test_fetch_observed_by_dir(gatherer, code):
    """
    Get waveforms based on given directory strucutre
    """
    assert gatherer.fetch_observed_by_dir(code) is None

    gatherer.config.paths["waveforms"] = "./test_data/test_mseeds"
    st = gatherer.fetch_observed_by_dir(code)
    assert st is not None
    assert len(st) == 3


def test_fetch_synthetic_by_dir(gatherer, code):
    """
    Get synthetics based on given directory strucutre
    """
    assert gatherer.fetch_synthetic_by_dir(code) is None

    gatherer.config.paths["synthetics"] = "./test_data/synthetics"
    st = gatherer.fetch_synthetic_by_dir(code)
    assert st is not None
    assert len(st) == 3

    # Testing out the synthetics only feature where obs data is searched for
    # using synthetic data structure
    assert gatherer.fetch_synthetic_by_dir(
        code, syn_cfgpath="waveforms") is None
    gatherer.config.paths["waveforms"] = "./test_data/synthetics"
    st = gatherer.fetch_synthetic_by_dir(code, syn_cfgpath="waveforms")
    assert st is not None
    assert len(st) == 3


def test_gather_event(gatherer, dataset_fid):
    """
    Ensure gatherer can get an event from the correct sources
    """
    assert gatherer.gather_event(append_focal_mechanism=False) is not None

    with ASDFDataSet(dataset_fid) as ds:
        gatherer.ds = ds
        gatherer.Client = None
        assert gatherer.gather_event(append_focal_mechanism=False) is not None


def test_append_focal_mechanism(event):
    """
    Try appending focal mechanism from GeoNet. GCMT doesn't have this regional
    event so we will need to test that separately.
    """
    del event.focal_mechanisms
    assert(len(event.focal_mechanisms) == 0)

    # Gather using the GeoNet client
    event = append_focal_mechanism_to_event(event, client="GEONET",
                                            overwrite_focmec=False,
                                            overwrite_event=False)
    assert(len(event.focal_mechanisms) != 0)
    m_rr = event.preferred_focal_mechanism().moment_tensor.tensor.m_rr
    assert(pytest.approx(m_rr, 1E-16) == -2.47938E16)


def test_get_gcmt_moment_tensor():
    """
    Just ensure that getting via GCMT works as intended using an example event
    """
    # Kaikoura Earthquake
    origintime = "2016-11-13T11:02:00"
    magnitude = 7.8

    cat = get_gcmt_moment_tensors(event=None, origintime=origintime,
                                    magnitude=magnitude)
    assert(len(cat) == 1)
    event = cat[0]
    assert hasattr(event, "focal_mechanisms")
    m_rr = event.preferred_focal_mechanism().moment_tensor.tensor.m_rr
    assert(pytest.approx(m_rr, 1E-20) == 3.56E20)

    return event


def test_get_usgs_moment_tensor():
    """
    Just ensure that getting via USGS works as intended using an example event
    """
    event = test_get_gcmt_moment_tensor()
    del event.focal_mechanisms

    cat = get_usgs_moment_tensors(event=event)
    assert(len(cat) == 1)
    event = cat[0]
    assert hasattr(event, "focal_mechanisms")

    m_rr = event.preferred_focal_mechanism().moment_tensor.tensor.m_rr
    assert(pytest.approx(m_rr, 1E-20) == 4.81E20)
