"""
Auxiliary criteria for weighting misfit windows to enhance or suppress certain
measurements in order to fairly assess misfit criteria. This is meant to
compliment the functionalities of Pyflex without having to directly edit the
Pyflex source code.

Functions should work in place on a Manager class to avoid having to pass in
all the different arguments from the Manager.
"""
from pyatoa import logger
from pyatoa.utils.calculate import abs_max


def compose_geographical_weights(cat, inv):
    """
    Create weights based on source-receiver distribution

    :type cat: obspy.event.Catalog
    :param cat: Catalog of events that should contain event locations
    :type inv: obspy.core.inventory.Inventory
    :param inv: Inventory of stations that should contain station locations
    :return:
    """
    raise NotImplementedError


def geographical_weights(mgmt, window, comp):
    """
    Up-weight stations and receivers that are sparsely distributed,
    down-weight stations and receivers that are densely distributed.
    Decreases the weight that dense distributions carry and keeps the inversion
    updates fair.

    :type mgmt: pyatoa.core.manager.Manager
    :param mgmt: Manager object that should already contain waveforms
    :type window: list of pyflex.window.Window
    :param window: list of window objects to check
    :type comp: str
    :param comp: component of the waveform to check the windows for
    :rtype: list of pyflex.window.Window
    :return: list of windows that have been suppressed
    """
    raise NotImplementedError


def category_weights(mgmt, window, comp):
    """
    Weight the measurement based on data category, e.g. direct arrival,
    Love wave, Rayleigh wave.

    :type mgmt: pyatoa.core.manager.Manager
    :param mgmt: Manager object that should already contain waveforms
    :type window: list of pyflex.window.Window
    :param window: list of window objects to check
    :type comp: str
    :param comp: component of the waveform to check the windows for
    :rtype: list of pyflex.window.Window
    :return: list of windows that have been suppressed
    """
    raise NotImplementedError


def window_by_amplitude(mgmt, window, comp):
    """
    Remove windows where amplitudes of the observed waveform fall below some
    threshold criteria. This was created in order to suppress direct arrival
    windows for long-period bandpasses, where e.g. P-wave direct arrivals are
    low-amplitude and elongated in time, sometimes leading to poor adjoint
    sources. This ensures that only the higher amplitude signals, e.g. surface
    waves, are included into the misfit metric.

    :type mgmt: pyatoa.core.manager.Manager
    :param mgmt: Manager object that should already contain waveforms
    :type window: list of pyflex.window.Window
    :param window: list of window objects to check
    :type comp: str
    :param comp: component of the waveform to check the windows for
    :rtype: list of pyflex.window.Window
    :return: list of windows that have been suppressed
    """
    suppressed_windows = []
    for win_ in window:
        waveform_peak = abs_max(mgmt.st_syn.select(component=comp)[0].data)
        window_peak = abs_max(mgmt.st_syn.select(
            component=comp)[0].data[win_.left:win_.right]
                              )
        # Check the waveform amplitudes
        if (abs(window_peak / waveform_peak) >
                mgmt.config.window_amplitude_ratio):
            suppressed_windows.append(win_)
        else:
            logger.info(
                "removing window due to global amplitude ratio: "
                f"{ abs(window_peak / waveform_peak):.2f} < "
                f"{mgmt.config.window_amplitude_ratio:.2f}")
            continue

    return suppressed_windows
