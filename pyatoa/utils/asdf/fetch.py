"""
Functions for extracting information from a Pyasdf ASDFDataSet object
"""
from pyatoa import logger
from pyatoa.utils.form import format_model_number, format_step_count
from pyflex.window import Window
from obspy import UTCDateTime


def windows_from_dataset(ds, net, sta, model, step, previous_step=False):
    """
    Returns misfit windows from an ASDFDataSet for a given model, step,
    network and station, as well as a count of windows returned.

    If given model and step are not present in dataset (e.g. during line search,
    new step), will try to search the previous step, which may or may not be
    contained in the previous model. 

    Returns windows as Pyflex Window objects which can be used in Pyadjoint or
    in the Pyatoa workflow.

    Note:
        Expects that windows are saved into the dataset at each model and step 
        such that there is a coherent structure within the dataset

    To do:
        If windows are calculated for a given model/step but e.g. something 
        fails and I change parameters and retry, those windows will still be 
        there, and will be re-retrieved, which is not ideal. Might have to 
        clean dataset before rerunning? or add some choice variable.

    :type ds: pyasdf.ASDFDataSet
    :param ds: ASDF dataset containing MisfitWindows subgroup
    :type net: str
    :param net: network code used to find the name of the misfit window
    :type sta: str
    :param sta: station code used to find the name of the misfit window
    :type model: int or str
    :param model: model number, will be formatted by the function
    :type step: int or str
    :param step: step count, will be formatted by the function
    :type check_previous: bool
    :param check_previous: if no windows are found for the given model, step,
        search the dataset for available windows from the previous step
    :rtype window_dict: dict
    :return window_dict: dictionary containing misfit windows, in a format
        expected by Pyatoa Manager class
    """
    # Ensure the tags are properly formatted
    model_number = format_model_number(model)
    step_count = format_step_count(step)
    windows = ds.auxiliary_data.MisfitWindows

    window_dict = {}    
    if previous_step:
        # Attempt to retrieve windows from previous model/step
        prev_windows = _return_windows_from_previous_step(windows=windows, 
                                                          model=model, 
                                                          step=step
                                                          )
        window_dict = dataset_windows_to_pyflex_windows(
            windows=prev_windows, network=net, station=sta
            )  
    else:
        if hasattr(windows, model_number) and \
                            hasattr(windows[model_number], step_count):
            # Attempt to retrieve windows from the given model/step
            logger.debug(f"searching for windows in {model_number}{step_count}")
            window_dict = dataset_windows_to_pyflex_windows(
                windows=windows[model_number][step_count], network=net, 
                station=sta
                )

    return window_dict


def dataset_windows_to_pyflex_windows(windows, network, station):
    """
    Convert the parameter dictionary of an ASDFDataSet MisfitWindow into a 
    dictionary of Pyflex Window objects, in the same format as Manager.windows

    Returns empty dict and 0 if no windows are found

    :type windows: pyasdf.utils.AuxiliaryDataAccessor
    :param windows: ds.auxiliary_data.MisfitWindows[model][step]
    :rtype window_dict: dict
    :return window_dict: dictionary of window attributes in the same format
        that Pyflex outputs
    :rtype num_windows: int
    :return num_windows: number of windows for a given model, step, net, sta
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

def _return_windows_from_previous_step(windows, model, step):
    """
    Given a model number and step count, find windows from the previous step
    count. If none are found for the given model, return the most recently
    available windows.

    Note: Assumes that windows are saved at each iteration! Even if fixed 
        windows are used.

    :type model: int or str
    :param model: the current model 
    :type step: int or str
    :param step: the current step
    :rtype: pyasdf.utils.AuxiliaryDataAccessor
    :return: ds.auxiliary_data.MisfitWindows
    """
    # Ensure were working with integer values for indexing
    if isinstance(model, str):
        model = int(model[1:])
    if isinstance(step, str):
        step = int(step[1:])

    # Get a flattened list of models and steps as unique tuples of integers
    iters = []
    steps = {m: windows[m].list() for m in windows.list()}
    for m, s in steps.items():
        for s_ in s:
            iters.append((int(m[1:]), int(s_[1:])))

    current = (model, step)
    if current in iters:
        prev_model, prev_step = iters[iters.index(current) - 1]
    else:
        # Wind back the step to see if there are any windows in this given model
        while step >= 0:
            if (model, step) in iters:
                prev_model, prev_step = model, step
                break
            step -= 1
        else:
            # If nothing is found return the most recent windows available
            prev_model, prev_step = iters[-1]

    prev_model = format_model_number(prev_model)
    prev_step = format_step_count(prev_step)
    logger.debug(f"searching for windows in {prev_model}{prev_step}")

    return windows[prev_model][prev_step]

