"""
Functions for extracting information from a Pyasdf ASDFDataSet object
"""
from pyatoa import logger
from pyatoa.utils.form import format_model_number, format_step_count
from pyflex.window import Window
from obspy import UTCDateTime


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
    steps = {m: misfit_windows[m].list() for m in misfit_windows.list()}
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

def windows_from_dataset(ds, net, sta, model, step, check_previous=True):
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
    if hasattr(windows, model_number) and \
                        hasattr(windows[model_number], step_count):
        # Attempt to retrieve windows from the given model/step
        logger.debug(f"searching for windows in {model_number}{step_count}")
        window_dict = dataset_windows_to_pyflex_windows(
            windows=windows[model_number][step_count], network=net, station=sta
            )
    else:
        if check_previous:
            # Attempt to retrieve windows from previous model/step
            prev_windows = _return_windows_from_previous_step(windows=windows, 
                                                              model=model, 
                                                              step=step
                                                              )
            window_dict = dataset_windows_to_pyflex_windows(
                windows=prev_windows, network=net, station=sta
                )  

    return window_dict


def sum_misfits(ds, model, step):
    """
    Misfits are stored in adjoint trace dictionaries, and are needed for
    Seisflows. This will sum the misfits and place them into a specific
    filepath.
    
    As per Tape (2010) Eq. 6, misfit for a single earthquake is given as:
    F^T_s(m) = (1/2) * (1/N_s) * sum[i=1:N_s] (F^T_i(m))
    
    N_s = total number of measurement windows for a source s
    ith window identified by a 
    source, station, component, period range, local window index

    The total misfit function F^T is given by Eq. 7
    F^T(m) = (1/S) * sum[s=1:S] (F^T_s(m))

    where S is the number of sources. The total misfit is calculated in
    the seisflows function: seisflows.workflow.inversion.write_misfit()    

    :type ds: pyasdf.ASDFDataSet
    :param ds: asdf dataset which must contain AdjointSource auxiliary data
    :type model: str 
    :param model: model number, e.g "m00"
    :type step: str
    :param step: step count, e.g. 's00', if None, defaults to latest step of
        given model. if False, assumes no step information present, only model
    :rtype summed_misfit: float
    :return summed_misfit: total misfit for a given dataset and model
    """
    misfits = []
    adjoint_sources = ds.auxiliary_data.AdjointSources[model][step]
    windows = ds.auxiliary_data.MisfitWindows[model][step]
    number_windows = len(windows)

    # collect the total misfit calculated by Pyadjoint
    for adjsrc in adjoint_sources.list():
        misfits.append(adjoint_sources[adjsrc].parameters["misfit_value"])

    return 0.5 * sum(misfits) / number_windows


def count_misfit_windows(ds, model, step, count_by_stations=False):
    """
    Figure out which stations contain which windows, return a dictionary
    which lists available components.
    Used by plotting scripts when annotating the number of misfit windows 

    :type ds: pyasdf.ADSFDataSet
    :param ds: dataset to count misfit windows from
    :type model: str
    :param model: model number e.g. 'm00'
    :type step: str or None or False
    :param step: step count, e.g. 's00', if None, defaults to latest step of
        given model. if False, assumes no step information present, only model
    :type count_by_stations: bool
    :param count_by_stations: if not, count by components
    :rtype counted_windows: dict
    :return counted_windows: station name as key, number of windows as value
    """
    assert(hasattr(ds.auxiliary_data.MisfitWindows, model)), \
        f"dataset has no misfit windows for model {model}"
    assert (hasattr(ds.auxiliary_data.MisfitWindows[model], step)), \
        f"dataset has no misfit windows for step {step} in model {model} "

    windows = ds.auxiliary_data.MisfitWindows[model][step]

    # Count up windows for each channel
    window_list = []
    for window in windows.list():
        window_list.append(windows[window].parameters['channel_id'])

    counted_windows = {}
    # Count by station, not by component
    if count_by_stations:
        stations = []
        # Remove the component tag from the name
        for comp in window_list:
            stations.append("{net}.{sta}".format(net=comp.split(".")[0],
                                                 sta=comp.split(".")[1])
                            )
        uniqueid = set(stations)
        quantity_out = stations
    # Count by component
    else:
        uniqueid = set(window_list)
        quantity_out = window_list

    # Count up the number of unique instances
    for id_ in uniqueid:
        counted_windows[id_] = quantity_out.count(id_)

    return counted_windows


def histogram_data(ds, model, data_type, step=None):
    """
    Returns values to be used in a statistics histogram. For example, a list
    containing all traveltime differences for every window, or amplitude
    differenes for a given dataset.
    
    :type ds: pyasdf.ADSFDataSet
    :param ds: dataset to count misfit windows from
    :type model: str
    :param model: model number e.g. 'm00'
    :type step: str or None or False
    :param step: step count, e.g. 's00', if None, defaults to latest step of
        given model. if False, assumes no step information present, only model
    :type data_type: str
    :param data_type: chosen data type to return histogram data for, available:
        'cc_shift_in_samples', 'cc_shift_in_seconds', 
    """
    from pyatoa.utils.calculate import amplitude_anomaly

    data_types = ['cc_shift_in_samples', 'cc_shift_in_seconds', 'amplitude']
    if data_type not in data_types:
        print(f"data_type must be in {data_types}")
        return
                         
    return_list = []
    # time shift statistics are kept in the misfit windows
    if "cc" in data_type:
        if hasattr(ds.auxiliary_data.MisfitWindows, model):
            windows = ds.auxiliary_data.MisfitWindows[model]
            if (step is not None) and (
                    hasattr(ds.auxiliary_data.MisfitWindows, step)):
                windows = windows[step]
            for window in windows:
                return_list.append(window.parameters[data_type])
    # amplitude information is recovered from waveforms
    elif data_type == "amplitude":
        for sta in ds.waveforms:
            # collect data for observed and synthetic traces
            st_obs = sta.observed 
            st_syn = None
            for tag in sta.get_waveform_tags():
                if model in tag:
                    st_syn = sta[tag]
            # get amplitude anomaly for each component
            if st_syn:
                for i in range(st_obs.count()):
                    # ensure you're comparing the correct components
                    obs = st_obs[i].data
                    syn = st_syn.select(
                            component=st_obs[i].get_id()[-1])[0].data
                    return_list.append(amplitude_anomaly(
                                       a=obs, b=syn, dt=st_obs[i].stats.delta)
                                       )
            else:
                continue 

    return return_list


def _count_waveforms(ds_waveforms, model):
    """
    Count the number of waveforms for each station for a given model, return
    a list with number of waveforms corresponding to the dataset station list
    
    :type ds_waveforms: pyasdf.ASDFDataSet.waveforms
    :param ds_waveforms: asdf dataset waveforms subdirectory
    :type model: str 
    :param model: model number, e.g "m00"
    :rtype syn_n: list of int
    :return syn_n: number of synthetic waveforms for a given station 
    :rtype obs_n: list of int
    :return obs_n: number of observed  waveforms for a given station
    """
    syn_n, obs_n = [], []
    for sta in ds_waveforms.list():
        try:
            syn_n.append(len(ds_waveforms[sta][f"synthetic_{model}"]))
        except KeyError:
            syn_n.append(0)
        try:
            obs_n.append(len(ds_waveforms[sta]["observed"]))
        except KeyError:
            obs_n.append(0)

    return syn_n, obs_n


def _misfit_info(ds_adjsrc):
    """
    Find the mean misfit for all adjoint sources in a given event
    Also find the min and max misfits and corresponding stations
    
    :type ds_adjsrc: pyasdf.ASDFDataSet.auxiliary_data.AdjointSources
    :param ds_adjsrc: asdf dataset adjoint sources subdirectory
    :rtype misfits: list of floats
    :return misfits: misfit values for a given dataset
    :rtype max_misfit: tuple of floats
    :return max_misfit: (max misfit, name of component with max misfit)
    :rtype min_minsfit: tuple of floats
    :return min_misfit: (min misfit, name of component with min misfit)
    """
    misfits = []
    for comp in ds_adjsrc.list():
        misfits.append(ds_adjsrc[comp].parameters["misfit_value"])

    # gather information from the list of misfits   
    # TO DO check this works
    max_misfit = (max(misfits), ds_adjsrc.list()[misfits.index(max(misfits))])
    min_misfit = (min(misfits), ds_adjsrc.list()[misfits.index(min(misfits))])
     
    return misfits, max_misfit, min_misfit

    
def _window_info(ds_windows):
    """
    Collect information from each of the misfit windows in a dataset
    
    :type ds_windows: pyasdf.ASDFDataSet.auxiliary_data.MisfitWindows
    :param ds_windows: asdf dataset misfit windows subdirectory
    :rtype cc_shifts_secs: list of floats
    :return cc_shifts_secs: individual values of cross correlation time shifts
        produced by Pyflex
    :rtype max_cc_values: list of floats
    :return max_cc_values: individual values of maximum cross correlation
        produced by Pyflex
    """
    cc_shift_secs = []
    max_cc_values = [] 
    for comp in ds_windows.list():
        cc_shift_secs.append(ds_windows[comp].parameters["cc_shift_in_seconds"])
        max_cc_values.append(ds_windows[comp].parameters["max_cc_value"])
    
    return cc_shift_secs, max_cc_values
    



