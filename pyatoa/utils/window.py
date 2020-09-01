"""
Auxiliary criteria for weighting misfit windows to enhance or suppress certain
measurements in order to fairly assess misfit criteria. This is meant to
compliment the functionalities of Pyflex without having to directly edit the
Pyflex source code.

Functions should work in place on a Manager class to avoid having to pass in
all the different arguments from the Manager.
"""
import numpy as np
from pyatoa import logger
from pyatoa.utils.calculate import abs_max


def zero_pad_then_window(ws, pad_by_fraction_of_npts=.2):
    """
    To address Pyflex throwing ValueErrors when source-receiver distances are

    .. note::
        Sept 1, 2020
        Work in progress, may not actually want to do this to avoid any
        near-source effects?

    :type ws: pyflex.WindowSelector
    :param ws: an already-filled window selector object that should
        be passed in from the Manager object
    :rtype: list of pyflex.Window
    :return: a list of Window objects, or an empty list if no windows found or
        the zero padding didnt work
    """
    raise NotImplementedError

    logger.warning("Pyflex has thrown a ValueError, most likely due to a small"
                   "source-receiver distance. Attempting to zero-pad waveforms"
                   "and re-run window selection")

    # We assume that these traces have already been standardized. These values
    # will be used to ensure that we can undo the zero-padding
    original_origintime = ws.observed.stats.starttime
    original_endtime = ws.observed.stats.endtime
    original_npts = ws.observed.stats.npts

    # Pad by a fraction of the trace length
    pad_width = int(original_npts * pad_by_fraction_of_npts)

    # Pad only the front of the data
    ws.observed.data = np.pad(ws.observed.data, (pad_width,), mode="constant")
    ws.observed.stats.starttime -= pad_width * ws.observed.stats.delta

    ws.observed.data = np.pad(ws.observed.data, (pad_width,), mode="constant")
    ws.observed.stats.starttime -= pad_width * ws.observed.stats.delta

    ws.select_windows()


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

    :type data: np.ndarray
    :param data: data array to query amplitude values from
    :type windows: list of pyflex.window.Window
    :param windows: list of window objects to check
    :type ratio: float
    :param ratio: percentage threshold of the peak value within a given window
        and the global peak value in the data array. Defaults to 0.2
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
