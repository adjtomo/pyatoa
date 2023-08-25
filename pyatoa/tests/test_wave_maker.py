"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import os
import json
import logging
import numpy as np
import matplotlib.pyplot as plt
from pyasdf import ASDFDataSet
from obspy import read, read_events, read_inventory

from pyatoa import Config, Manager, logger
from pyatoa.visuals.wave_maker import WaveMaker


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
    return Config(event_id="2018p130600")


@pytest.fixture
def mgmt(config, event, st_obs, st_syn, inv):
    """
    A manager filled with data and that has progressed through the workflow
    """
    mgmt = Manager(config=config, event=event, st_obs=st_obs, st_syn=st_syn, 
                   inv=inv)
    mgmt.standardize()
    mgmt.preprocess(remove_response=True, output="DISP")
    mgmt.window()
    mgmt.measure()
    return mgmt


@pytest.fixture
def wm(mgmt):
    """
    WaveMaker object to create waveform figures
    """
    wm = WaveMaker(mgmt=mgmt)
    wm.setup_plot(dpi=100, figsize=(140, 60))
    return wm


def test_setup_plot(wm):
    """
    Make sure the correct number of axes are established
    """
    assert(len(wm.axes) == len(wm.twaxes) == 3)


def test_plot_waveforms(wm):
    """
    Plot waveforms by themselves
    """
    for i, obs in enumerate(wm.st_obs):
        comp = obs.stats.channel[-1]
        syn = wm.st_syn.select(component=comp)[0]
        wm.plot_waveforms(ax=wm.axes[i], obs=obs, syn=syn)
        wm.plot_amplitude_threshold(ax=wm.axes[i], obs=obs)
    plt.close()


def test_plot_stalta(wm):
    """
    Plot STALTA waveforms by themselves
    """
    for i, stalta in enumerate(wm.staltas.values()):
        wm.plot_stalta(ax=wm.axes[i], stalta=stalta)
    plt.close()


def test_plot_windows(wm):
    """
    Test plotting windows
    """
    for i, windows in enumerate(wm.windows.values()):
        wm.plot_windows(ax=wm.axes[i], windows=windows)
    plt.close() 


def test_plot_adjsrcs(wm):
    """
    Test plotting adjoint sources
    """
    for i, adjsrc in enumerate(wm.adjsrcs.values()):
        wm.plot_adjsrcs(ax=wm.twaxes[i], adjsrc=adjsrc)
    plt.close()


def test_plot_rejected_windows(wm):
    """
    Test plotting rejected windows 
    """
    for i, rejwins in enumerate(wm.rejwins.values()):
        wm.plot_rejected_windows(ax=wm.axes[i], rejwin=rejwins)
    plt.close()


def test_plot(wm):
    """
    Test the main plotting function
    """
    wm.plot(show=False, save=False)


