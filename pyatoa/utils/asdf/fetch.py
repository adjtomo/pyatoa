"""
Functions for extracting information from a Pyasdf ASDFDataSet object
"""
from pyatoa import logger
from pyatoa.utils.form import format_iter, format_step
from pyflex.window import Window
from obspy import UTCDateTime


def windows_from_dataset(ds, net, sta, iteration, step_count, 
                         return_previous_step=False):
    """
    Returns misfit windows from an ASDFDataSet for a given iteration, step,
    network and station, as well as a count of windows returned.

    If given iteration and step are not present in dataset (e.g. during line 
    search, new step), will try to search the previous step, which may or 
    may not be contained in the previous iteration. 

    Returns windows as Pyflex Window objects which can be used in Pyadjoint or
    in the Pyatoa workflow.

    Note:
        Expects that windows are saved into the dataset at each iteration and 
        step such that there is a coherent structure within the dataset

    To do:
        If windows are calculated for a given iteration/step but e.g. something 
        fails and I change parameters and retry, those windows will still be 
        there, and will be re-retrieved, which is not ideal. Might have to 
        clean dataset before rerunning? or add some choice variable.

    :type ds: pyasdf.ASDFDataSet
    :param ds: ASDF dataset containing MisfitWindows subgroup
    :type net: str
    :param net: network code used to find the name of the misfit window
    :type sta: str
    :param sta: station code used to find the name of the misfit window
    :type iteration: int or str
    :param iteration: current iteration, will be formatted by the function
    :type step_count: int or str
    :param step_count: step count, will be formatted by the function
    :type return_previous_step: bool
    :param return_previous_step: search the dataset for available windows from
        the previous iteration/step given the current iteration/step
    :rtype window_dict: dict
    :return window_dict: dictionary containing misfit windows, in a format
        expected by Pyatoa Manager class
    """
    # Ensure the tags are properly formatted
    iteration = format_iter(iteration)
    step_count = format_step(step_count)
    windows = ds.auxiliary_data.MisfitWindows

    window_dict = {}    
    if return_previous_step:
        # Retrieve windows from previous iter/step
        prev_windows = return_windows_from_previous_step(windows=windows,
                                                         iteration=iteration,
                                                         step_count=step_count
                                                         )
        window_dict = dataset_windows_to_pyflex_windows(windows=prev_windows,
                                                        network=net, station=sta
                                                        )
    else:
        if hasattr(windows, iteration) and \
                            hasattr(windows[iteration], step_count):
            # Attempt to retrieve windows from the given iter/step
            logger.debug(f"searching for windows in {iteration}{step_count}")
            window_dict = dataset_windows_to_pyflex_windows(
                windows=windows[iteration][step_count], network=net, 
                station=sta
                )

    return window_dict


def dataset_windows_to_pyflex_windows(windows, network, station):
    """
    Convert the parameter dictionary of an ASDFDataSet MisfitWindow into a 
    dictionary of Pyflex Window objects, in the same format as Manager.windows

    Returns empty dict and 0 if no windows are found

    :type windows: pyasdf.utils.AuxiliaryDataAccessor
    :param windows: ds.auxiliary_data.MisfitWindows[iter][step]
    :type network: str
    :param network: network of the station related to the windows
    :type station: str
    :param station: station related to the windows
    :rtype window_dict: dict
    :return window_dict: dictionary of window attributes in the same format
        that Pyflex outputs
    :rtype num_windows: int
    :return num_windows: number of windows for a given iter, step, net, sta
    """
    window_dict, _num_windows = {}, 0
    for window_name in windows.list():
        net, sta, comp, n = window_name.split("_")

        # Check the title of the misfit window to see if applicable
        if (net == network) and (sta == station):
            par = windows[window_name].parameters

            # Create a Pyflex Window object
            window = Window(
                left=par["left_index"], right=par["right_index"],
                center=par["center_index"], dt=par["dt"],
                time_of_first_sample=UTCDateTime(par["time_of_first_sample"]),
                min_period=par["min_period"], channel_id=par["channel_id"]
            )

            # We cant initiate these parameters so set them after the fact
            # If data changed, should recalculate with Window._calc_criteria()
            setattr(window, "dlnA", par["dlnA"])
            setattr(window, "cc_shift", par["cc_shift_in_samples"])
            setattr(window, "max_cc_value", par["max_cc_value"])

            # Save windows into the dictionary labelled by component
            if comp in window_dict.keys():
                # Either append to existing entry
                window_dict[comp] += [window]
            else:
                # Or create the first entry
                window_dict[comp] = [window]
            _num_windows += 1

    logger.debug(f"{_num_windows} window(s) found in dataset for "
                 f"{network}.{station}")
    return window_dict


def return_windows_from_previous_step(windows, iteration, step_count,
                                       continuous_search=False):
    """
    Given an iteration and step count, find windows from the previous step
    count. If none are found for the given iteration, return the most recently
    available windows.

    Note: Assumes that windows are saved at each iteration! Even if fixed 
        windows are used.

    :type iteration: int or str
    :param iteration: the current iteration 
    :type step_count: int or str
    :param step_count: the current step count
    :type continuous_search: bool
    :param continuous_search: allow the search to wind back through all
        available iterations/step_counts until an acceptable set of windows is
        found. Preferable to not allow this, but I wrote it so might as well
        keep it around.
    :rtype: pyasdf.utils.AuxiliaryDataAccessor
    :return: ds.auxiliary_data.MisfitWindows
    """
    # Ensure we're working with integer values for indexing, e.g. 's00' -> 0
    if isinstance(iteration, str):
        iteration = int(iteration[1:])
    if isinstance(step_count, str):
        step_count = int(step_count[1:])

    # Get a flattened list of iters and steps as unique tuples of integers
    iters = []
    steps = {i: windows[i].list() for i in windows.list()}
    for i, s in steps.items():
        for s_ in s:
            iters.append((int(i[1:]), int(s_[1:])))

    current = (iteration, step_count)
    if current in iters:
        # Windows have been found in the previous iteration/step_count
        prev_iter, prev_step = iters[iters.index(current) - 1]
    elif continuous_search:
        # Wind back the step to see if any windows in this given iteration
        while step_count >= 0:
            if (iteration, step_count) in iters:
                prev_iter, prev_step = iteration, step_count
                break
            step_count -= 1
        else:
            # If nothing is found return the most recent windows available
            prev_iter, prev_step = iters[-1]
    else:
        # Dont allow continous search, simply return empty dictionary
        logger.debug(f"no previous windows found w.r.t "
                     f"{format_iter(iteration)}{format_step(step_count)}")
        return {}

    prev_iter = format_iter(prev_iter)
    prev_step = format_step(prev_step)
    logger.debug(f"searching for windows in {prev_iter}{prev_step}")

    return windows[prev_iter][prev_step]

