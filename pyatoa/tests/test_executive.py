"""
Test the functionalities of the Executive class which runs many Managers
"""
import os
import pytest
import numpy as np
from pyatoa import Config, Executive, logger


# Turn off the logger for tests
logger.propogate = False
logger.setLevel("CRITICAL")


pytest.skip(allow_module_level=True)

@pytest.fixture
def events():
    """
    A list of GeoNet event ids used for gathering metadata from GeoNet client
    """
    event_ids = ["2018p130600", "2012p242656"]
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
    # syn_path = "./test_data/test_executive/{}"
    # synthetics = [os.path.abspath(syn_path.format(_)) for _ in events]

    syn_path = [os.path.abspath("./test_data/test_executive/")]

    cfg = Config(iteration=1, step_count=0, min_period=10, max_peropd=30,
                 client="GEONET", pyflex_preset="default",
                 adj_src_type="cc_traveltime_misfit",
                 paths={"synthetics": syn_path})
    return cfg


def test_executive_single_event_single_station_no_concurrent(tmpdir, config,
                                                             events, stations):
    """
    Attempt a single event single station processing without using concurrency
    """
    exc = Executive(event_ids=events[0], station_codes=stations[0],
                    config=config, cwd=tmpdir.strpath)
    misfit = exc.process_station(f"{events[0]}{exc.cat}{stations[0]}")
    assert(pytest.approx(misfit, .001) == 1.6696)


def test_executive_single_event_single_station(tmpdir, config, events,
                                               stations):
    """
    Attempt a single event single station processing with concurrency
    """
    exc = Executive(event_ids=events[0], station_codes=stations[0],
                    config=config, cwd=tmpdir.strpath)
    misfits = exc.process()
    misfit = misfits[events[0]][stations[0]]
    assert(pytest.approx(misfit, .001) == 1.6696)


def test_executive_single_event_multi_station(tmpdir, config, events,
                                               stations):
    """
    Attempt a single event multi station processing with concurrency
    """
    exc = Executive(event_ids=events[0], station_codes=stations,
                    config=config, cwd=tmpdir.strpath, max_events=1,
                    max_stations=os.cpu_count())
    misfits = exc.process()
    assert(len(misfits) == 1)
    assert(len(misfits[events[0]]) == len(stations))
    misfit = misfits[events[0]][stations[4]]
    assert(pytest.approx(misfit, .001) == 0.76983)


# def test_executive_multi_event_multi_station(tmpdir, config, events,
#                                                stations):
#     """
#     !!! This test is causing my computer to crash, not sure why, must rework
#
#     Attempt a single event multi station processing with concurrency.
#     Only do 2 events and 2 stations max to avoid crashing out system
#     """
#     exc = Executive(event_ids=events, station_codes=stations,
#                     config=config, cwd=tmpdir.strpath, max_events=2,
#                     max_stations=2)
#     misfits = exc.process()
#     assert(len(misfits) == 2)
#     assert(len(misfits[events[1]]) == len(stations))
#     misfit = misfits[events[1]][stations[2]]
#     assert(pytest.approx(misfit, .001) == 12.242)
