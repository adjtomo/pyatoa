"""
Test the functionalities of the Pyatoa Gatherer class
"""
import pytest
from obspy import read_events
from pyasdf import ASDFDataSet
from pyatoa import Config
from pyatoa.core.gatherer import Gatherer, GathererNoDataException


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
    return Config(event_id=event_id, iteration=1, step_count=0,
                  save_to_ds=False)


@pytest.fixture
def gatherer(config, origintime):
    """
    The Gatherer which is responsible for gathering data.
    """
    return Gatherer(config=config, origintime=origintime)


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
    with pytest.raises(GathererNoDataException):
        gatherer.gather_event()

    with ASDFDataSet(dataset_fid) as ds:
        gatherer.ds = ds
        gatherer.Client = None
        assert gatherer.gather_event() is not None

