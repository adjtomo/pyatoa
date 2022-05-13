"""
Test the functionalities of the Executive class which runs many Managers
"""
import pytest
from pyatoa import Config, Executive, logger


# Turn off the logger for tests
logger.propogate = False
logger.setLevel("CRITICAL")


@pytest.fixture
def events():
    """
    A list of GeoNet event ids used for gathering metadata from GeoNet client
    """
    event_ids = ["2018p130600", "2012p242656", "2017p015402", "2228901"]
    return event_ids


@pytest.fixture
def stations():
    """
    A list of GeoNet station codes used for gathering metadata and waveforms
    """
    location = "*"
    channel = "HH?"
    codes = ["NZ.BFZ", "NZ.KNZ", "NZ.COVZ", "NZ.PXZ", "NZ.WEL"]
    codes = [f"{c}.{location}.{channel}" for c in codes]

    return codes


@pytest.fixture
def config(events):
    """
    A preset Config object that specifies where to grab data from, which
    already exists in the test data directory
    """
    syn_path = "./test_data/test_executive/{}"
    synthetics = [syn_path.format(_) for _ in events]

    cfg = Config(iteration=1, step_count=0, min_period=10, max_peropd=30,
                 client="GEONET", pyflex_preset="default",
                 adj_src_type="cc_traveltime_misfit",
                 paths={"synthetics": synthetics})
    return cfg


def test_single_event_single_station_no_concurrent(config, events, stations):
    """
    Attempt a single event single stationprocessing
    """
    exc = Executive(event_ids=events[0], station_codes=stations[0],
                    config=config)
    misfit = exc.process_station(f"{events[0]}-{stations[0]}")
    assert(pytest.approx(misfit, .001) == 1.6696)


def test_executive_single_event_single_station(config, events, stations):
    """
    Attempt a single event single stationprocessing
    """
    exc = Executive(event_ids=events[0], station_codes=stations[0],
                    config=config)
    misfits = exc.process()

    pytest.set_trace()
