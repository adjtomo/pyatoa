"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import os
import json
import logging
import numpy as np
from IPython import embed
from pyasdf import ASDFDataSet
from pyatoa import logger, Config, Manager
from pyatoa.utils.write import write_stream_sem
from pyatoa.core.seisflows.pyaflowa import Pyaflowa
from obspy import read, read_events, read_inventory

logger.setLevel("DEBUG")


@pytest.fixture
def st_obs():
    """
    Raw observed waveforms from station NZ.BFZ.HH? for New Zealand event
    2018p130600 (GeoNet event id)
    """
    return read("./test_data/test_obs_data_NZ_BFZ_2018p130600.ascii")


@pytest.fixture
def st_syn():
    """
    Synthetic data stream generated using Specfem3D and focal mechanism for
    2018p130600. Minimum resolved period roughly 10s.
    """
    return read("./test_data/test_syn_data_NZ_BFZ_2018p130600.ascii")


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
def inv():
    """
    StationXML information for station NZ.BFZ.HH?
    """
    return read_inventory("./test_data/test_dataless_NZ_BFZ.xml")


@pytest.fixture
def config():
    """
    Default Pyatoa Config object
    """
    return Config(event_id="2018p130600", client="GEONET")


@pytest.fixture
def mgmt(config, event, st_obs, st_syn, inv):
    """
    A manager filled with data but pre-workflow
    """
    mgmt = Manager(config=config, event=event, st_obs=st_obs, st_syn=st_syn, 
                   inv=inv)
    mgmt.config.model = 0
    mgmt.config.step = 0
    mgmt.standardize()
    mgmt.preprocess()
    mgmt.window()
    mgmt.measure()

    return mgmt


@pytest.fixture
def ds(tmpdir, mgmt):
    """
    Fill an ASDFDataSet to check fixed windows
    """
    ds = ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5"))
    mgmt.ds = ds
    mgmt.write()
    return ds


@pytest.fixture
def pyaflowa(tmpdir):
    """
    Initate an empty Pyaflowa class
    """
    pars = json.load(open("./test_data/test_seisflows_parameters.json"))
    paths = {"PYATOA_IO": os.path.join(tmpdir, "pyatoa.io"),
             "WORKDIR": "./test_data", "SPECFEM_DATA": "./test_data"
             }
    return Pyaflowa(pars, paths)


def test_check_for_fixed_windows(pyaflowa, ds):
    """
    Ensure parameter checks work as advertised
    """
    assert pyaflowa._check_for_fixed_windows(ds)
    pyaflowa.iteration = 2
    assert(not pyaflowa._check_for_fixed_windows(ds))

def test_prepare_event(tmpdir, pyaflowa):
    """
    Check that the Config parameter and paths are set up correctly
    """
    cfg, paths = pyaflowa.prepare_event(cwd=tmpdir, event_id="2018p130600")
    assert(os.path.join(tmpdir, "traces", "syn") in cfg.cfgpaths["synthetics"])

def test_process_event(tmpdir, pyaflowa, st_obs, st_syn, cat, inv):
    """
    Run the event processing
    """
    cfg, paths = pyaflowa.prepare_event(cwd=tmpdir, event_id="2018p130600")

    # create a STATIONS file in the tmpdir as Pyaflowa will be looking for it
    os.makedirs(os.path.join(tmpdir, "DATA"), exist_ok=True)
    with open(os.path.join(tmpdir, "DATA", "STATIONS"), "w") as f:
        f.write("BFZ NZ -40.6796 176.2462 0.0 0.0")

    # write synthetic streams to tmpdir in Specfem3D format
    # hardcode known time offset in the synthetic data
    os.makedirs(cfg.cfgpaths["synthetics"][0], exist_ok=True)
    write_stream_sem(st_syn, "d", path=cfg.cfgpaths["synthetics"][0], 
                     time_offset=-20)

    # write observed waveforms into a SEED formatted directory structure
    for tr in st_obs:
        s = tr.stats
        t = tr.stats.starttime
        path = f"{t.year}/{s.network}/{s.station}/{s.channel}.D/"
        path = os.path.join(cfg.cfgpaths["waveforms"][0], path)
        os.makedirs(path, exist_ok=True)

        fid = f"{tr.id}.D.{t.year}.{t.julday:0>3}"
        tr.write(os.path.join(path, fid), format="MSEED")

    # write StationXML into a SEED formatted directory structure
    resp = os.path.join(tmpdir, "RESPONSE")
    cfg.cfgpaths["responses"].append(resp)
    path = os.path.join(resp, f"{s.station}.{s.network}")
    os.makedirs(path, exist_ok=True)
    for channel in inv.get_contents()["channels"]:
        cha = channel.split(".")[-1]
        fid = f"RESP.{channel}"
        inv.select(channel=cha).write(os.path.join(path, fid), 
                                      format="STATIONXML")

    # Run event processing and ensure data can be found and processed
    with ASDFDataSet(os.path.join(tmpdir, "process_event_dataset.h5")) as ds:
        # Set the event so that FDSN query not required
        ds.events = cat
        status = pyaflowa.process_event(ds, cfg, paths)
        pyaflowa.prepare_eval_grad(ds, paths)

    # Quick check that prepare_eval_grad() worked
    assert(os.path.exists(os.path.join(tmpdir, "DATA", "STATIONS_ADJOINT")))
    assert(os.path.exists(
                os.path.join(tmpdir, "traces", "adj", "NZ.BFZ.BXE.adj"))
    )


