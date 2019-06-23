#!/usr/bin/env python3
"""
Undoubtedly, mapping is very dependent on region and scale, so
things like fault lines, and landmarks can not be hardcoded into
the package. This function serves as an auxiliary entry point for the
User to implement their own mapping goodies.
"""
import pkg_resources

import numpy as np
import matplotlib as mpl

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def plot_hikurangi_trench(m):
    """
    Trace the hikurangi trench from a coordinate file

    :type m: Basemap
    :param m: basemap object
    :type path_: str
    :param path_: pathway to hikurangi trench coordinates
    """
    trenchcoords = np.load(
        pkg_resources.resource_filename(
            __name__, "fault_coordinates/hikurangi_trench.npz"))
    lats = trenchcoords['LAT']
    lons = trenchcoords['LON']
    x, y = m(lons, lats)

    # interpolate points to make a smoother curve
    xprime = np.flip(x, axis=0)
    yprime = np.flip(y, axis=0)
    xprimenew = np.linspace(x.min(), x.max(), 100)
    yprimenew = np.interp(xprimenew, xprime, yprime)

    m.plot(xprimenew, yprimenew, ':', linewidth=2.25, color='k', alpha=0.5,
           zorder=2)


def plot_geonet_active_faults(m):
    """
    Plot onshore and offshore fault coordinate files taken from GeoNet

    :type m: Basemap
    :param m: basemap object
    :type path_: str
    :param path_: pathway to hikurangi trench coordinates
    """
    int_fid = "fault_coordinates/north_island_550_641_{}.npz"
    for tag in ["onshore", "offshore"]:
        fid = pkg_resources.resource_filename(__name__, int_fid.format(tag))
        active_faults = np.load(fid)
        lats = active_faults['LAT']
        lons = active_faults['LON']
        faults = active_faults['FAULT']

        for i in range(faults.min(), faults.max()+1, 1):
            indices = np.where(faults == i)
            x, y = m(lons[indices], lats[indices])
            m.plot(x, y, '-.', linewidth=1, color='k', zorder=2, alpha=0.25)
