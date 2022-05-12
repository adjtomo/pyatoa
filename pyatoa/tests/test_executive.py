"""
Test the functionalities of the Executive class which runs many Managers
"""
import pytest
from pyatoa import Config, Executive, logger


# Turn off the logger for tests
logger.propogate = False
logger.setLevel("DEBUG")


@pytest.fixture
def config():
    """
    A preset Config object that specifies where to grab data from
    """
    cfg = Config(iteration=1, step_count=0, min_period=10, max_peropd=30,
                 client="GEONET", pyflex_preset="default",
                 adj_src_type="cc_traveltime_misfit")
    return cfg


@pytest.fixture
def events():
    """
    A list of GeoNet event ids used for gathering metadata from GeoNet client
    """
    event_ids = ["2018p130600", "2012p242656", "2017p015402", "2228901",
                 "2013p617227", "2016p275188", "2017p087060", "2019p927023"
                 ]
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


def test_single_event_single_station_no_concurrent(config, events, stations):
    """
    Attempt a single event single stationprocessing
    """
    exc = Executive(event_ids=events[0], station_codes=stations[0],
                    config=config)
    misfit = exc.process_station(f"{events[0]}-{stations[0]}")

    pytest.set_trace()


# def test_executive_single_event_single_station(config, events, stations):
#     """
#     Attempt a single event single stationprocessing
#     """
#     exc = Executive(event_ids=events[0], station_codes=stations[0])
#     misfits = exc.process()
#
#     pytest.set_trace()
