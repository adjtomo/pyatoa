"""
Functions for extracting information from a Pyasdf ASDFDataSet object
"""
from pyatoa import logger
from obspy import UTCDateTime
from fnmatch import filter as fnf
from pyflex.window import Window
from pyadjoint.adjoint_source import AdjointSource
from pyatoa.utils.form import format_iter, format_step


def load_windows(ds, net, sta, iteration, step_count, return_previous=False):
    """
    Returns misfit windows from an ASDFDataSet for a given iteration, step,
    network and station, as well as a count of windows returned.

    If given iteration and step are not present in dataset (e.g. during line 
    search, new step), will try to search the previous step, which may or 
    may not be contained in the previous iteration. 

    Returns windows as Pyflex Window objects which can be used in Pyadjoint or
    in the Pyatoa workflow.

    .. note::
        Expects that windows are saved into the dataset at each iteration and 
        step such that there is a coherent structure within the dataset

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
    :type return_previous: bool
    :param return_previous: search the dataset for available windows
        from the previous iteration/step given the current iteration/step
    :rtype window_dict: dict
    :return window_dict: dictionary containing misfit windows, in a format
        expected by Pyatoa Manager class
    """
    # Ensure the tags are properly formatted
    iteration = format_iter(iteration)
    step_count = format_step(step_count)
    windows = ds.auxiliary_data.MisfitWindows

    window_dict = {}    
    if return_previous:
        # Retrieve windows from previous iter/step
        prev_windows = previous_windows(windows=windows, iteration=iteration,
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


def load_adjsrcs(ds, net, sta, iteration, step_count):
    """
    Load adjoint sources from a pyasdf ASDFDataSet and return in the format
    expected by the Manager class, that is a dictionary of adjoint sources

    :type ds: pyasdf.ASDFDataSet
    :param ds: ASDF dataset containing MisfitWindows subgroup
    :type net: str
    :param net: network code used to find the name of the adjoint source
    :type sta: str
    :param sta: station code used to find the name of the adjoint source
    :type iteration: int or str
    :param iteration: current iteration, will be formatted by the function
    :type step_count: int or str
    :param step_count: step count, will be formatted by the function
    :rtype: dict
    :return: dictionary containing adjoint sources, in a format expected by
        Pyatoa Manager class
    """
    # Ensure the tags are properly formatted before using them for access
    iteration = format_iter(iteration)
    step_count = format_step(step_count)
    adjsrcs = ds.auxiliary_data.AdjointSources[iteration][step_count]

    adjsrc_dict = {}
    # Use fnmatch filter to find all adjoint sources that match net/sta code
    for adjsrc_tag in fnf(adjsrcs.list(), f"{net}_{sta}_*"):
        component = adjsrc_tag[-1].upper()  # e.g. 'Z'

        # Build the adjoint source based on the parameters that were parsed in
        parameters = adjsrcs[adjsrc_tag].parameters
        assert(component == parameters["component"][-1]), (
            "AdjointSource tag does not match the component listed in the "
            "parameter dictionary when it should.")

        # Adjoint sources are time-reversed when saved into the dataset, so
        # reverse them back when returning to Manager. Also remove time axis.
        parameters["adjoint_source"] = adjsrcs[adjsrc_tag].data[()][:, 1][::-1]

        # Convert back from str to UTCDateTime object
        parameters["starttime"] = UTCDateTime(parameters["starttime"])

        # The parameter dicionary will have all the keywords necessary
        adjsrc_dict[component] = AdjointSource(**parameters)

    return adjsrc_dict


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
    :rtype: dict
    :return: dictionary of window attributes in the same format that Pyflex 
        outputs
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


def previous_windows(windows, iteration, step_count):
    """
    Given an iteration and step count, find windows from the previous step
    count. If none are found for the given iteration, return the most recently
    available windows.

    .. note:: 
        Assumes that windows are saved at each iteration.

    :type windows: pyasdf.utils.AuxiliaryDataAccessor
    :param windows: ds.auxiliary_data.MisfitWindows[iter][step]
    :type iteration: int or str
    :param iteration: the current iteration
    :type step_count: int or str
    :param step_count: the current step count
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
        # If windows have already been added to the auxiliary data
        prev_iter, prev_step = iters[iters.index(current) - 1]
    else:
        # Wind back the step to see if there are any windows for this iteration
        while step_count >= 0:
            if (iteration, step_count) in iters:
                prev_iter, prev_step = iteration, step_count
                break
            step_count -= 1
        else:
            # If nothing is found return the most recent windows available
            prev_iter, prev_step = iters[-1]

    # Format back into strings for accessing auxiliary data
    prev_iter = format_iter(prev_iter)
    prev_step = format_step(prev_step)

    logger.debug(f"most recent windows: {prev_iter}{prev_step}")

    return windows[prev_iter][prev_step]
