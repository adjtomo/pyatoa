"""
Test the Inspector class and its ability to generate dataframes for bulk
analyses of an inversion
"""
import pytest
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
def inspector(asdf_dataset_fid):
    insp = Inspector()
    insp.append(asdf_dataset_fid)
    return insp


def test_append(asdf_dataset_fid):
    """
    Inspector with data
    """
    insp = Inspector()
    insp.append(asdf_dataset_fid)


def test_discover(test_data):
    """
    Make sure Inspector can find HDF5 files and read them in
    """
    insp = Inspector()
    insp.discover(path=test_data)




