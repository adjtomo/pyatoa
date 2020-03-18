"""
Auxiliary functionality for windowing criteria to enhance or replace the
functionalities of Pyflex without having to edit the Pyflex source code.

Functions should work in place on a Manager class to avoid having to pass in
all the different arguments from the Manager.
"""
from pyatoa import logger
from pyatoa.utils.calculate import abs_max


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
                f"{ abs(window_peak / waveform_peak)} < "
                f"{mgmt.config.window_amplitude_ratio}")
            continue

    return suppressed_windows
