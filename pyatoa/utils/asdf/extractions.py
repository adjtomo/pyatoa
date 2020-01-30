"""
Functions for extracting information from a Pyasdf ASDFDataSet object
"""


def windows_from_ds(ds, net, sta, model=None):
    """
    If misfit windows are to be fixed, then window information needs to be saved
    between calls of Pyatoa. Fortunately, Pyatoa stores misfit window
    information into MisfitWindow auxiliary data subgroups. This function
    accesses this information and creates window dictionaries in the format
    expected by the Pyatoa Manager class.
    These windows only contain enough information so that Pyflex cooperates,
    Pyadjoint can use them, and Pyatoa can generate plots.

    :raises AttributeError: if no windows are found for the given parameters

    :type ds: pyasdf.ASDFDataSet
    :param ds: ASDF dataset containing MisfitWindows subgroup
    :type net: str
    :param net: network code used to find the name of the misfit window
    :type sta: str
    :param sta: station code used to find the name of the misfit window
    :type model: str or None
    :param model: model number, e.g. "m00". If None, most recent windows used
    :rtype window_dict: dict
    :return window_dict: dictionary containing misfit windows, in a format
        expected by Pyatoa Manager class
    """
    from pyflex.window import Window

    if not model:
        model = ds.auxiliary_data.MisfitWindows.list()[-1]
    misfit_windows = ds.auxiliary_data.MisfitWindows[model]

    # Pyatoa expects the Manager class windows as a dictionary with keys
    # corresponding to components, each item is then a list, containing
    # Pyflex Window objects
    window_dict = {}
    for window_name in misfit_windows.list():
        net_, sta_, comp_, n_ = window_name.split("_")
        # Check the title of the misfit window to see if applicable
        if (net == net_) and (sta == sta_):
            par = misfit_windows[window_name].parameters

            # Create the misfit window
            window_ = Window(
                left=par["left_index"], right=par["right_index"],
                center=par["center_index"], dt=par["dt"],
                time_of_first_sample=par["time_of_first_sample"],
                min_period=par["min_period"], channel_id=par["channel_id"]
            )

            # Pyflex is weird that it doesn't allow one to set all values
            # in __init__, so we need to manual set some values after init
            window_.dlnA = par["dlnA"]
            window_.cc_shift = par["cc_shift_in_samples"]
            window_.max_cc_value = par["max_cc_value"]

            # Either append to existing entry
            if comp_ in window_dict.keys():
                window_dict[comp_] += [window_]
            # Or create the first entry
            else:
                window_dict[comp_] = [window_]
    
    if not window_dict:
        raise AttributeError(f"No windows found for {model}.{net}.{sta}")

    return window_dict


def sum_misfits(ds, model, station=None):
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
    :type station: str
    :param station: station code to check, e.g. 'NZ.BFZ'
    :rtype summed_misfit: float
    :return summed_misfit: total misfit for a given dataset and model
    """
    misfits = []
    adjoint_sources = ds.auxiliary_data.AdjointSources[model]

    # Sum misfits only for a given station
    if station:
        try:
            number_windows = count_misfit_windows(ds, model, True)[station]
        except KeyError:
            return
        for adjsrc in adjoint_sources.list():
            if adjoint_sources[adjsrc].parameters['station_id'] == station:
                misfits.append(
                    adjoint_sources[adjsrc].parameters['misfit_value']
                )

    # Sum misfits for all given adjoint sources
    else:
        number_windows = len(ds.auxiliary_data.MisfitWindows[model])
        # collect the total misfit calculated by Pyadjoint
        for adjsrc in adjoint_sources.list():
            misfits.append(adjoint_sources[adjsrc].parameters["misfit_value"])

    return 0.5 * sum(misfits)/number_windows


def misfit_stats(ds, model, include_lists=False):
    """
    Extract misfit statistics from a dataset for a given model.
    Return information as a dictionary object for easy access.    
    Used primarily to write misfit stats into ASDF Dataset auxiliary data.
 
    :type ds: pyasdf.ASDFDataSet
    :param ds: asdf dataset which must contain AdjointSource auxiliary data
    :type model: str 
    :param model: model number, e.g "m00"
    :type include_lists: bool
    :param include_lists: include list information containing statistics for
        individual stations, components etc. Not always necessary
    :rtype stats: dict
    :return stats: a dictionary of statistics values 
    """
    # collect relevant information
    ds_adjsrc = ds.auxiliary_data.AdjointSources[model]
    if not hasattr(ds.auxiliary_data.MisfitWindows, model):
        latest_model = ds.auxiliary_data.MisfitWindows.list()[-1]
        ds_windows = ds.auxiliary_data.MisfitWindows[latest_model]
    else:
        ds_windows = ds.auxiliary_data.MisfitWindows[model]
    
    syn_n, obs_n = _count_waveforms(ds.waveforms, model) 
    misfits, max_misfit, min_misfit = _misfit_info(ds_adjsrc) 
    cc_shift_secs, max_cc_values = _window_info(ds_windows) 
    
    # save stats in a dictionary
    # commented out sections are lists that I don't think are necessary 
    stats = {
        "number_stations": len(ds.waveforms.list()),
        "number_syn_waveforms": sum(syn_n),
        "number_obs_waveforms": sum(obs_n),
        "number_adjoint_sources": len(ds_adjsrc),
        "min_misfit": min_misfit[0],
        "min_misfit_component": min_misfit[1],
        "max_misfit": max_misfit[0],
        "max_misfit_component": max_misfit[1],
        "average_misfit": sum(misfits)/len(misfits),
        "number_misfit_windows": len(ds_windows),
            }

    # lists of information that may or may not be relevant
    if include_lists:
        lists = {
            "stations": ds.waveforms.list(),
            "syn_waveforms": syn_n,
            "obs_waveforms": obs_n,
            "adjoint_sources": ds_adjsrc.list(),
            "misfits" : misfits,
            "cc_shift_in_seconds": cc_shift_secs,
            "misfit_windows": ds_windows.list(),
            "max_cc_values": max_cc_values
                }
        for item in lists.keys():
            stats[item] = lists[item]

    return stats


def count_misfit_windows(ds, model, count_by_stations=False):
    """
    Figure out which stations contain which windows, return a dictionary
    which lists available components.
    Used by plotting scripts when annotating the number of misfit windows 

    :type ds: pyasdf.ADSFDataSet
    :param ds: dataset to count misfit windows from
    :type model: str
    :param model: model number e.g. 'm00'
    :type count_by_stations: bool
    :param count_by_stations: if not, count by components
    :rtype counted_windows: dict
    :return counted_windows: station name as key, number of windows as value
    """
    components = []
    # Determine where to count windows from
    if hasattr(ds.auxiliary_data.MisfitWindows, model):
        windows = ds.auxiliary_data.MisfitWindows[model]
    else:
        latest_model = ds.auxiliary_data.MisfitWindows.list()[-1] 
        windows = ds.auxiliary_data.MisfitWindows[latest_model]

    # Count up windows for each channel
    for window in windows.list():
        components.append(windows[model][window].parameters['channel_id'])

    counted_windows = {}
    # Count by component
    if count_by_stations:
        stations = []
        # Remove the component tag from the name
        for comp in components:
            stations.append("{net}.{sta}".format(net=comp.split(".")[0],
                                                 sta=comp.split(".")[1])
                            )
        uniqueid = set(stations)
        quantity_out = stations
    # Count by station
    else:
        uniqueid = set(components)
        quantity_out = components

    # Count up the number of unique instances
    for id_ in uniqueid:
        counted_windows[id_] = quantity_out.count(id_)

    return counted_windows


def histogram_data(ds, model, data_type):
    """
    Returns values to be used in a statistics histogram. For example, a list
    containing all traveltime differences for every window, or amplitude
    differenes for a given dataset.
    
    :type ds: pyasdf.ADSFDataSet
    :param ds: dataset to count misfit windows from
    :type model: str
    :param model: model number e.g. 'm00'
    :type data_type: str
    :param data_type: chosen data type to return histogram data for, available:
        'cc_shift_in_samples', 'cc_shift_in_seconds', 
    """
    from pyatoa.utils.tools.calculate import amplitude_anomaly

    data_types = ['cc_shift_in_samples', 'cc_shift_in_seconds', 'amplitude']
    if data_type not in data_types:
        print(f"data_type must be in {data_types}")
        return
                         
    return_list = []
    # time shift statistics are kept in the misfit windows
    if "cc" in data_type:
        if hasattr(ds.auxiliary_data.MisfitWindows, model):
            for window in ds.auxiliary_data.MisfitWindows[model]:
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
    

def _count_by_station(ds, model):
    """
    TO DO: write a function that returns lists corresponding to each station 
    that count: number of misfit windows, number of adjoint sources
    """
    number_windows, number_adjsrcs = [],[]
    stations = ds.waveforms.list()
    
    raise NotImplementedError


