"""
Test the functionalities of the Pyatoa Gatherer class
"""
import pytest
from obspy import read_events
from obspy.clients.fdsn import Client
from pyasdf import ASDFDataSet
from pyatoa import Config
from pyatoa.core.gatherer import (ExternalGetter, InternalFetcher, Gatherer,
                                  get_gcmt_moment_tensor, append_focal_mechanism
                                  )


@pytest.fixture
def code():
    """
    Example NZ station code
    """
    return "NZ.BFZ.??.HH?"


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
def config():
    """
    Default Pyatoa Config object
    """
    return Config(event_id="2018p130600", client="GEONET", iteration=1,
                  step_count=0, synthetics_only=False, save_to_ds=False)


@pytest.fixture
def internal_fetcher(config, origintime):
    """
    The fetcher class which searches internal directory structures and
    ASDFDataSets. Normally these functionalities are inherited directly via
    the Gatherer class but for testing purposes we initiate separately.
    """
    int_fetch = InternalFetcher()
    int_fetch.origintime = origintime
    int_fetch.ds = None
    int_fetch.config = config
    return int_fetch


@pytest.fixture
def external_getter(config, origintime):
    """
    The getter class which queries FDSN webservices for data.
    Normally these functionalities are inherited directly via
    the Gatherer class but for testing purposes we initiate separately.
    """
    ext_get = ExternalGetter()
    ext_get.origintime = origintime
    ext_get.ds = None
    ext_get.config = config
    ext_get.Client = Client("GEONET")
    return ext_get


@pytest.fixture
def gatherer(config, origintime):
    """
    The Gatherer which is responsible for gathering data.
    """
    return Gatherer(config=config, origintime=origintime)


def test_external_getter_no_client(external_getter, code):
    """
    Make sure that external getter functions always returns None if no Client
    is provided
    """
    external_getter.Client = None

    assert external_getter.event_get() is None
    assert external_getter.station_get(code) is None
    assert external_getter.obs_waveform_get(code) is None


def test_event_get(external_getter):
    """
    Ensure that querying for an event works. Only tests GeoNet gathering.
    """
    # Get by event id from GeoNet
    assert external_getter.event_get() is not None

    # Get by origin time from GeoNet
    external_getter.config.event_id = None
    assert external_getter.event_get() is not None


def test_station_get(external_getter, code):
    """
    Ensure that querying for a station works. Only tests GeoNet gathering.
    """
    net, sta, loc, cha = code.split('.')

    inv = external_getter.station_get(code, level="response")
    assert inv[0].code == net
    assert inv[0][0].code == sta
    assert hasattr(inv[0][0][0], "response")


def test_obs_waveform_get(external_getter, code):
    """
    Ensure that querying for a waveforms works. Only tests GeoNet gathering.
    """
    net, sta, loc, cha = code.split('.')

    st = external_getter.obs_waveform_get(code)
    assert(len(st) == 3)

    stats = st.select(component="Z")[0].stats
    assert stats.network == net
    assert stats.station == sta


def test_asdf_event_fetch(internal_fetcher, dataset_fid):
    """
    Get event from an ASDFDataSet.
    """
    with ASDFDataSet(dataset_fid) as ds:
        internal_fetcher.ds = ds
        internal_fetcher.asdf_event_fetch()


def test_asdf_station_fetch(internal_fetcher, dataset_fid, code):
    """
    Get station from an ASDFDataSet
    """
    with ASDFDataSet(dataset_fid) as ds:
        internal_fetcher.ds = ds
        inv = internal_fetcher.asdf_station_fetch(code)
        assert len(inv[0][0]) == 3


def test_asdf_waveform_fetch(internal_fetcher, dataset_fid, code, config):
    """
    Get waveforms from an ASDFDataSet
    """
    with ASDFDataSet(dataset_fid) as ds:
        internal_fetcher.ds = ds
        for tag in [config.observed_tag, config.synthetic_tag]:
            st = internal_fetcher.asdf_waveform_fetch(code, tag)
            assert len(st) == 3


def test_fetch_event_by_dir(internal_fetcher, event_id):
    """
    Get response based on given directory structure
    """
    # No Config path means fetching returns nada
    assert internal_fetcher.fetch_event_by_dir(event_id) is None

    internal_fetcher.config.paths["events"] = "./test_data/"
    event = internal_fetcher.fetch_event_by_dir(
                                  event_id, event_id_prefix="test_CMTSOLUTION_")

    assert event is not None
    assert event_id in event.resource_id.id

def test_fetch_resp_by_dir(internal_fetcher, code):
    """
    Get response based on given directory structure
    """
    # No Config path means fetching returns nada
    assert internal_fetcher.fetch_resp_by_dir(code) is None

    internal_fetcher.config.paths["responses"] = "./test_data/test_seed"
    inv = internal_fetcher.fetch_resp_by_dir(code)
    assert inv is not None
    assert hasattr(inv[0][0][0], "response")


def test_fetch_obs_by_dir(internal_fetcher, code):
    """
    Get waveforms based on given directory strucutre
    """
    assert internal_fetcher.fetch_obs_by_dir(code) is None

    internal_fetcher.config.paths["waveforms"] = "./test_data/test_mseeds"
    st = internal_fetcher.fetch_obs_by_dir(code)
    assert st is not None
    assert len(st) == 3


def test_fetch_syn_by_dir(internal_fetcher, code):
    """
    Get synthetics based on given directory strucutre
    """
    assert internal_fetcher.fetch_syn_by_dir(code) is None

    internal_fetcher.config.paths["synthetics"] = "./test_data/synthetics"
    st = internal_fetcher.fetch_syn_by_dir(code)
    assert st is not None
    assert len(st) == 3

    # Testing out the synthetics only feature where obs data is searched for
    # using synthetic data structure
    assert internal_fetcher.fetch_syn_by_dir(code,
                                             syn_cfgpath="waveforms") is None
    internal_fetcher.config.paths["waveforms"] = "./test_data/synthetics"
    st = internal_fetcher.fetch_syn_by_dir(code, syn_cfgpath="waveforms")
    assert st is not None
    assert len(st) == 3


def test_obs_waveform_fetch(internal_fetcher, dataset_fid, code):
    """
    Test the mid level fetching function which chooses whether to search via
    ASDFDataSet or directory structure
    """
    # Temporarily set config paths to check the correct order of funcs called
    internal_fetcher.config.paths["waveforms"] = "./test_data/test_mseeds"
    assert internal_fetcher.obs_waveform_fetch(code) is not None

    internal_fetcher.config.paths["waveforms"] = "./test_data/synthetics"
    internal_fetcher.config.synthetics_only = True
    assert internal_fetcher.obs_waveform_fetch(code) is not None

    internal_fetcher.config.paths["waveforms"] = None

    with ASDFDataSet(dataset_fid) as ds:
        internal_fetcher.ds = ds
        assert internal_fetcher.obs_waveform_fetch(code) is not None


def test_syn_waveform_fetch(internal_fetcher, dataset_fid, code):
    """
    Test the mid level fetching function which chooses whether to search via
    ASDFDataSet or directory structure
    """
    internal_fetcher.config.paths["synthetics"] = "./test_data/synthetics"
    assert internal_fetcher.syn_waveform_fetch(code) is not None
    internal_fetcher.config.paths["synthetics"] = None

    with ASDFDataSet(dataset_fid) as ds:
        internal_fetcher.ds = ds
        assert internal_fetcher.syn_waveform_fetch(code) is not None


def test_station_fetch(internal_fetcher, dataset_fid, code):
    """
    Test the mid level fetching function which chooses whether to search via
    ASDFDataSet or directory structure
    """
    internal_fetcher.config.paths["responses"] = "./test_data/test_seed"
    assert internal_fetcher.station_fetch(code) is not None
    internal_fetcher.config.paths["responses"] = None

    with ASDFDataSet(dataset_fid) as ds:
        internal_fetcher.ds = ds
        assert internal_fetcher.station_fetch(code) is not None


def test_gather_event(gatherer, dataset_fid):
    """
    Ensure gatherer can get an event from the correct sources
    """
    assert gatherer.gather_event(try_fm=False) is not None

    with ASDFDataSet(dataset_fid) as ds:
        gatherer.ds = ds
        gatherer.Client = None
        assert gatherer.gather_event(try_fm=False) is not None


def test_append_focal_mechanism(gatherer, event):
    """
    Try appending focal mechanism from GeoNet. GCMT doesn't have this regional
    event so we will need to test that separately.
    """
    del event.focal_mechanisms
    assert(len(event.focal_mechanisms) == 0)

    # Gather using the GeoNet client
    event = append_focal_mechanism(event, client="GEONET", overwrite=True)
    assert(len(event.focal_mechanisms) != 0)


def test_get_gcmt_moment_tensor():
    """
    Just ensure that getting via GCMT works as intended using an example event
    """
    # Kaikoura Earthquake
    origintime = "2016-11-13T11:02:00"
    magnitude = 7.8

    event = get_gcmt_moment_tensor(origintime, magnitude)
    assert hasattr(event, "focal_mechanisms")
