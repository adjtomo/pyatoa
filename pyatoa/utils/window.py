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


def geographical_weights(windows):
    """
    Up-weight stations and receivers that are sparsely distributed,
    down-weight stations and receivers that are densely distributed.
    Decreases the weight that dense distributions carry and keeps the inversion
    updates fair.

    :type windows: list of pyflex.window.Window
    :param windows: list of window objects to check
    :rtype: list of pyflex.window.Window
    :return: list of windows that have been suppressed
    """
    raise NotImplementedError


def category_weights(windows):
    """
    Weight the measurement based on data category, e.g. direct arrival,
    Love wave, Rayleigh wave.

    :type window: list of pyflex.window.Window
    :param window: list of window objects to check
    :rtype: list of pyflex.window.Window
    :return: list of windows that have been suppressed
    """
    raise NotImplementedError


def reject_on_global_amplitude_ratio(data, windows, ratio=0.2):
    """
    Reject windows where peak amplitude falls below some threshold value. 

    This was created in order to suppress windows containing long period direct 
    arrivals, which were creating high-frequency adjoint sources.

    :type array: np.ndarray
    :param array: data array to query amplitude values from
    :type windows: list of pyflex.window.Window
    :param windows: list of window objects to check
    :rtype: tuple of lists of pyflex.window.Window
    :return: lists of accepted and rejected windows
    """
    accepted_windows, rejected_windows = [], []
    for win in windows:
        waveform_peak = abs_max(data)
        window_peak = abs_max(data[win.left:win.right])
        # Check the waveform amplitudes
        if abs(window_peak / waveform_peak) > ratio:
            accepted_windows.append(win)
        else:
            rejected_windows.append(win)
            
    logger.info("rejection based on global amplitude ratio removed "
                f"{len(rejected_windows)} windows"
                )

    return accepted_windows, rejected_windows
