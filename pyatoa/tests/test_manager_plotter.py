"""
Test the functionalities of the Pyaflowa Manager class
"""
import pytest
import os
import json
import logging
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
from pyasdf import ASDFDataSet
from obspy import read, read_events, read_inventory

from pyatoa import Config, Manager, logger
from pyatoa.visuals.manager_plotter import ManagerPlotter


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
    A manager filled with data and that has progressed through the workflow
    """
    mgmt = Manager(config=config, event=event, st_obs=st_obs, st_syn=st_syn, 
                   inv=inv)
    mgmt.standardize()
    mgmt.preprocess()
    mgmt.window()
    mgmt.measure()
    return mgmt


def test_setup_plot(mgmt):
    """
    Make sure the correct number of axes are established
    """
    mp = ManagerPlotter(mgmt=mgmt)
    f, axes, twaxes = mp.setup_plot()
    assert(len(axes) == len(twaxes) == 3)


def test_plot_waveforms(mgmt):
    """
    Plot waveforms by themselves
    """
    mp = ManagerPlotter(mgmt=mgmt)
    _, axes, _ = mp.setup_plot()
    for i, obs in enumerate(mp.st_obs):
        comp = obs.stats.channel[-1]
        syn = mp.st_syn.select(component=comp)[0]
        mp.plot_waveforms(ax=axes[i], obs=obs, syn=syn)
        mp.plot_amplitude_threshold(ax=axes[i], obs=obs)
    plt.show()
    plt.close()


def test_plot_stalta(mgmt):
    """
    Plot STALTA waveforms by themselves
    """
    mp = ManagerPlotter(mgmt=mgmt)
    _, axes, _ = mp.setup_plot()
    for i, stalta in enumerate(mp.staltas.values()):
        mp.plot_stalta(ax=axes[i], stalta=stalta)
    plt.show()
    plt.close()


def test_plot_windows(mgmt):
    """
    Plot STALTA waveforms by themselves
    """
    mp = ManagerPlotter(mgmt=mgmt)
    _, axes, _ = mp.setup_plot()
    for i, windows in enumerate(mp.windows.values()):
        mp.plot_windows(ax=axes[i], windows=windows)
    plt.show()
    plt.close() 


def test_plot_adjsrcs(mgmt):
    """
    Plot STALTA waveforms by themselves
    """
    mp = ManagerPlotter(mgmt=mgmt)
    _, _, twaxes = mp.setup_plot()
    for i, adjsrc in enumerate(mp.adjsrcs.values()):
        mp.plot_adjsrcs(ax=twaxes[i], adjsrc=adjsrc)
    plt.show()
    plt.close()


def test_plot_rejected_windows(mgmt):
    """
    Plot STALTA waveforms by themselves
    """
    mp = ManagerPlotter(mgmt=mgmt)
    _, axes, _ = mp.setup_plot()
    for i, rejwins in enumerate(mp.rejected_windows.values()):
        mp.plot_rejected_windows(ax=axes[i], windows=rejwins)
    plt.show()
    plt.close()


def test_plot(mgmt):
    """
    Plot STALTA waveforms by themselves
    """
    mp = ManagerPlotter(mgmt=mgmt, show=True, save=None)
    mp.plot()


