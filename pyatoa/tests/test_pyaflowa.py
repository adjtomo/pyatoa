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
from pyatoa.core.pyaflowa import IO, PathStructure, Pyaflowa


class Dict(dict):
    """
    Dictionary with gettable attributes to mimic SeisFlows PAR and PATH variable
    """
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


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
    path.WORKDIR = tmpdir
    path.SCRATCH = os.path.join(path.WORKDIR, "scratch")
    path.PREPROCESS = os.path.join(path.SCRATCH, "preprocess")
    path.SOLVER = os.path.join(path.SCRATCH, "solver")

    return path


def test_path_structure_standalone(tmpdir, source_name):
    """
    Test that the path structure standalone creates the expected
    directory structure. Just make sure things work.
    """
    path_structure = PathStructure(structure="standalone", 
                                   workdir=tmpdir.strpath
                                   )
    # Missing files will throw a FileNotFoundError when trying to set up 
    # the standalone path structure
    with pytest.raises(FileNotFoundError):
        path_structure.format(source_name=source_name)

    # Make the required directory and file to allow the automated pathing
    os.mkdir(os.path.join(tmpdir, source_name))
    open(os.path.join(tmpdir, source_name, "STATIONS") ,"w") 
    path_structure.format(source_name=source_name)

    # Ensure that the correct file structure has been generated 
    for fid in ["datasets", "figures", "input", "logs"]:
        assert(os.path.exists(os.path.join(tmpdir, fid)))


def test_path_structure_seisflows(tmpdir, source_name, PAR, PATH):
    """
    Test that seisflows directory structure matches what's expected
    """
    path_structure = PathStructure(structure="seisflows", sfpaths=PATH)
    with pytest.raises(FileNotFoundError):
        path_structure.format(source_name=source_name)

    # Make required directory and file structure for auto pathing
    open(os.path.join(PATH.SOLVER, source_name, "DATA", "STATIONS"), "w")
    path_structure.format(source_name=source_name)

    # Ensure that directories are made where they need to be 
    for fid in ["datasets", "figures", "logs"]:
        assert(os.path.exists(os.path.join(PATH.PREPROCESS, fid)))


def test_pyaflowa_init(PAR, PATH):
    """
    Test that Pyaflowa can instantiate correctly
    """
    pyaflowa = Pyaflowa(structure="seisflows", sfpaths=PATH, sfpar=PAR)

    # Ensure that the initiated config takes values from the parameter file
    assert(pyaflowa.config.min_period == PAR.MIN_PERIOD)
    assert(pyaflowa.config.max_period == PAR.MAX_PERIOD)
    assert(pyaflowa.config.client == PAR.CLIENT)


def test_pyaflowa_setup(source_name, PAR, PATH):
    """
    Test the one-time setup of Pyaflowa which creates the IO object
    """
    pyaflowa = Pyaflowa(structure="seisflows", sfpaths=PATH, sfpar=PAR)

    # Requirement that STATION file exists for Pyaflowa to run setup
    with pytest.raises(FileNotFoundError):
        pyaflowa.path_structure.format(source_name=source_name)
    open(os.path.join(PATH.SOLVER, source_name, "DATA", "STATIONS"), "w")

    # SeisFlows usually takes care of placing source files into the data
    # directory, so we need to do it manually here
    src = os.path.join("test_data", "test_CMTSOLUTION_2018p130600")
    dst = os.path.join(PATH.SOLVER, source_name, "DATA", "CMTSOLUTION")
    shutil.copy(src, dst)

    # Initiate Pyaflowa which will create directory structure, read in source
    # file and create an ASDFDataSet 
    io = pyaflowa.setup(source_name=source_name)

    # Simple check to make sure event id is set correctly and event reading
    # machinery is working
    assert(io.config.event_id == source_name)
    assert(os.path.exists(io.paths.ds_file))
    with ASDFDataSet(io.paths.ds_file) as ds:
        assert(source_name in ds.events[0].resource_id.id)


def test_pyaflowa_process_station(tmpdir, seisflows_workdir, seed_data,
                                  source_name, station_name, PAR, PATH):
    """
    Test the single station processing function 
    """
    # Turn off client to avoid searching FDSN, force local data search
    PAR.CLIENT = None
    PATH.DATA = tmpdir.strpath
    pyaflowa = Pyaflowa(structure="seisflows", sfpaths=PATH, sfpar=PAR,
                        iteration=1, step_count=0)

    # Copy working directory to tmpdir to avoid creating unnecessary files
    shutil.copytree(src=seisflows_workdir, dst=os.path.join(tmpdir, "scratch"))
    shutil.copytree(src=seed_data, dst=os.path.join(tmpdir, "seed"))

    # Set up the same machinery as process_event()
    io = pyaflowa.setup(source_name)
    with ASDFDataSet(io.paths.ds_file) as ds:
        mgmt = Manager(ds=ds, config=io.config)
        mgmt, io = pyaflowa.process_station(mgmt=mgmt, code="NZ.BFZ.??.???", 
                                            io=io)

    assert(io.nwin == mgmt.stats.nwin == 3)
    assert(io.misfit == pytest.approx(65.39037, .001))


def test_pyaflowa_process_event(tmpdir, seisflows_workdir, seed_data,
                                source_name, PAR, PATH):
    """
    Test the actual machinery of processing a full station using Pyaflowa.
    Don't use FDSN gathering for simplicity, simply place all the requisite
    files in their requisite places and let Pyaflowa find things on its own.

    .. note::
        This is meant to mimic SeisFlows3.preprocess.pyatoa.prepare_eval_grad()
    """
    # Turn off client to avoid searching FDSN, force local data search
    PAR.CLIENT = None
    PATH.DATA = tmpdir.strpath
    pyaflowa = Pyaflowa(structure="seisflows", sfpaths=PATH, sfpar=PAR,
                        iteration=1, step_count=0)

    # Copy working directory to tmpdir to avoid creating unnecessary files
    shutil.copytree(src=seisflows_workdir, dst=os.path.join(tmpdir, "scratch"))
    shutil.copytree(src=seed_data, dst=os.path.join(tmpdir, "seed"))
   
    misfit = pyaflowa.process_event(source_name=source_name, 
                                    fix_windows=PAR.FIX_WINDOWS,
                                    event_id_prefix=PAR.SOURCE_PREFIX)

    # Make sure misfit value comes out right to a certain rounding error
    assert(misfit == pytest.approx(10.898, .001))

    # Make sure finalization files are created, including adj srcs and sta file
    paths = pyaflowa.path_structure.format(source_name=source_name)
    assert(len(glob.glob(os.path.join(paths.adjsrcs, "*"))) == 3)
    assert(os.path.exists(os.path.join(paths.data, "STATIONS_ADJOINT")))

