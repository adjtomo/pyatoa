"""
Test the Inspector class and its ability to generate dataframes for bulk
analyses of an inversion
"""
import pytest
import numpy as np
from pyasdf import ASDFDataSet
from pyatoa import Inspector, logger

# Turn off logger for tests
logger.propagate = False
logger.setLevel("CRITICAL")


@pytest.fixture
def test_data():
    return "./test_data"


@pytest.fixture
def asdf_dataset_fid():
    return "./test_data/2018p130600.h5"


@pytest.fixture
def seisflows_inspector(test_data):
    """
    An Inspector filled in by a SeisFlows workflow used to test some of the
    more complex functions
    """
    tag = "test_inspector"
    insp = Inspector()
    insp.read(path=test_data, fmt="csv", tag=tag)
    return insp


@pytest.fixture
def inspector(asdf_dataset_fid):
    insp = Inspector()
    insp.append(asdf_dataset_fid)
    return insp


def test_append(asdf_dataset_fid):
    """
    Make sure that the Inspector can append a single event and put the correct
    data in the correct place
    """
    insp = Inspector()
    insp.append(asdf_dataset_fid)
    assert(insp.events == "2018p130600")
    assert(insp.stations == "BFZ")
    assert(insp.iterations == "i01")
    assert(insp.evaluations == 1)


def test_discover(test_data):
    """
    Make sure Inspector can find HDF5 files generally and read them in.
    """
    insp = Inspector()
    insp.discover(path=test_data, ignore_symlinks=True)
    assert(insp.events == "2018p130600")
    assert(insp.stations == "BFZ")
    assert(insp.iterations == "i01")
    assert(insp.evaluations == 1)


def test_extend(inspector):
    """
    Make sure that you can extend the current inspector with the widnows of
    another. In this example we just use the same inspector twice
    Also tests the copy function to make sure that the extension doesn't
    affect both inspectors
    """
    insp_a = inspector.copy()
    insp_b = inspector.copy()

    insp_a.extend(insp_b.windows)
    assert(insp_a.evaluations == 2)
    assert(insp_b.evaluations == 1)


def test_read_write_csv(tmpdir, inspector):
    """
    Test the read and write functions with both an empty and a filled inspector
    """
    filled_insp = inspector.copy()
    empty_insp = Inspector()

    empty_insp.save(path=tmpdir, fmt="csv", tag="empty_inspector")
    filled_insp.save(path=tmpdir, fmt="csv", tag="filled_inspector")

    check_insp = Inspector()
    check_insp.read(path=tmpdir, fmt="csv", tag="filled_inspector")

    # just a simple check to make sure we read in the same
    assert(filled_insp.srcrcv.distance_km[0] ==
           check_insp.srcrcv.distance_km[0])


def test_isolate(seisflows_inspector):
    """
    Test the isolate function to grab specific data from a filled inspector.
    Test by checking number of windows at each isolation call
    """
    insp = seisflows_inspector.copy()
    assert(len(insp.windows) == 714)
    assert(len(insp.isolate(iteration="i01", step_count="s00")) == 206)
    assert(len(insp.isolate(station="BKZ")) == 21)
    assert(len(insp.isolate(channel="BXE")) == 234)
    assert(len(insp.isolate(component="Z")) == 240)
    assert(np.shape(insp.isolate(iteration="i01",
                                 step_count="s00",
                                 keys="cc_shift_in_seconds")) == (206, 1)
           )
    assert(np.shape(insp.isolate(iteration="i01",
                                 step_count="s00",
                                 exclude="cc_shift_in_seconds")) == (206, 16)
           )
    assert(np.shape(insp.isolate(unique_key="cc_shift_in_seconds")) == (714, 7))


def test_nwin(seisflows_inspector):
    """
    Test the number of windows function
    """
    check_list = [206, 204, 207, 97]
    insp = seisflows_inspector.copy()
    assert(insp.nwin().nwin.to_list() == check_list)


def test_misfit(seisflows_inspector):
    """
    Test the misfit calculation function
    """
    check_list = [6.624891459456695, 6.435833032352672, 6.453905909311119,
                  7.09944815562882]
    insp = seisflows_inspector.copy()
    for i, misfit in enumerate(insp.misfit().misfit.to_list()):
        assert(misfit == pytest.approx(check_list[i], .00001))


def test_stats(seisflows_inspector):
    """
    Test the per-level stats calculations
    """
    insp = seisflows_inspector.copy()
    misfit = insp.stats(level="event").misfit.i01.s03["2013p617227"]
    assert(misfit == pytest.approx(17.1739208, .000001))
    misfit = insp.stats(level="station").misfit.i01.s03["BFZ"]
    assert(misfit == pytest.approx(26.143711, .0001))


def test_minmax(seisflows_inspector):
    """
    Test the minmax printing function
    """
    insp = seisflows_inspector.copy()
    minmax = insp.minmax(pprint=False)
    assert(minmax["dlnA_max"] == pytest.approx(1.03947, .00001))


def test_compare_events(seisflows_inspector):
    """
    Test inter-event comparisons
    """
    insp = seisflows_inspector.copy()
    assert(insp.compare_events().diff_nwin["2013p617227"] == -110)


def test_get_models(seisflows_inspector):
    """
    Test the model state tracker
    """
    states = [0, 1, -1, -1]
    insp = seisflows_inspector.copy()
    for i, state in enumerate(insp.models.state):
        assert(state == states[i])


def test_get_unique_models(seisflows_inspector):
    """
    Test convenience function that finds accepted models only
    """
    insp = seisflows_inspector.copy()
    assert(np.shape(insp.get_unique_models()) == (2, 6))
