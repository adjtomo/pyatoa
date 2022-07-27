"""
Test the workflow management class Pyaflowa and its supporting classes
PathStructure and IO
"""
import os
import glob
import shutil
import pytest
import yaml
from pyasdf import ASDFDataSet
from pyatoa import Manager
from pyatoa.core.pyaflowa import Pyaflowa


class Dict(dict):
    """
    Dictionary with gettable attributes to mimic SeisFlows PAR and PATH variable
    """
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)


@pytest.fixture
def source_name():
    """
    Example path structure for Pyaflowa
    """
    return "2018p130600"


@pytest.fixture
def station_name():
    """
    Example path structure for Pyaflowa
    """
    return "NZ.BFZ.??.HH?"


@pytest.fixture
def parameter_file():
    """
    Example path structure for Pyaflowa
    """
    return os.path.join("./test_data", "test_seisflows_parameters.yaml")


@pytest.fixture
def seisflows_workdir():
    """
    Example path structure for Pyaflowa
    """
    return os.path.join("./test_data", "test_seisflows_workdir", "scratch")


@pytest.fixture
def seed_data():
    """
    Example path structure for Pyaflowa
    """
    return os.path.join("./test_data", "test_seed")


@pytest.fixture
def PAR(parameter_file):
    """
    An example PAR variable that will be passed in from SeisFlows to Pyaflowa.
    Mimics the Dict objects defined in SeisFlows
    """
    return Dict(yaml.safe_load(open(parameter_file)))


@pytest.fixture
def PATH(tmpdir, PAR):
    """
    An example PATH variable that will be passed in from SeisFlows to Pyaflowa
    Mimics the Dict objects defined in SeisFlows
    """
    path = Dict(PAR.PATHS)

    # These are the paths that are auto set by the System module in SeisFlows
    # and required by Pyaflowa path structure class
    path.WORKDIR = tmpdir.strpath
    path.SCRATCH = os.path.join(path.WORKDIR, "scratch")
    path.PREPROCESS = os.path.join(path.SCRATCH, "preprocess")
    path.SOLVER = os.path.join(path.SCRATCH, "solver")

    return path


def test_pyaflowa_setup(source_name, PAR, PATH):
    """
    Test the one-time setup of Pyaflowa which creates the IO object
    """
    pyaflowa = Pyaflowa(sfpath=PATH, sfpar=PAR)

    # Requirement that STATION file exists for Pyaflowa to run setup
    pyaflowa.paths = pyaflowa.paths.format(source_name=source_name)
    open(os.path.join(PATH.SOLVER, source_name, "DATA", "STATIONS"), "w")

    # SeisFlows usually takes care of placing source files into the data
    # directory, so we need to do it manually here
    src = os.path.join("test_data", "test_CMTSOLUTION_2018p130600")
    dst = os.path.join(PATH.SOLVER, source_name, "DATA", "CMTSOLUTION")
    shutil.copy(src, dst)

    # Initiate Pyaflowa which will create directory structure, read in source
    # file and create an ASDFDataSet
    pyaflowa.setup(source_name=source_name, iteration=1, step_count=0)

    # Simple check to make sure event id is set correctly and event reading
    # machinery is working
    assert(pyaflowa.config.event_id == source_name)
    assert(os.path.exists(pyaflowa.paths.ds_file))
    with ASDFDataSet(pyaflowa.paths.ds_file) as ds:
        assert(source_name in ds.events[0].resource_id.id)


def test_pyaflowa_process(tmpdir, seisflows_workdir, seed_data,
                                 source_name, station_name, PAR, PATH):
    """
    Test processing in serial with a single source receiver
    """
    # Turn off client to avoid searching FDSN, force local data search
    PAR.CLIENT = None
    PATH.DATA = tmpdir.strpath
    pyaflowa = Pyaflowa(sfpath=PATH, sfpar=PAR)

    # Copy working directory to tmpdir to avoid creating unnecessary files
    shutil.copytree(src=seisflows_workdir, dst=os.path.join(tmpdir, "scratch"))
    shutil.copytree(src=seed_data, dst=os.path.join(tmpdir, "seed"))

    # Set up the same machinery as process_event()
    pyaflowa.setup(source_name, iteration=1, step_count=0)
    misfit = pyaflowa.process(nproc=1)
    assert(pytest.approx(misfit, .01) == 184.38)

