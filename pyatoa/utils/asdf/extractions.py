"""
Functions for extracting information from a Pyasdf ASDFDataSet object
"""


def sum_misfits(ds, model):
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
    :rtype summed_misfit: float
    :return summed_misfit: total misfit for a given dataset and model
    """
    adjoint_sources = ds.auxiliary_data.AdjointSources[model]
    number_windows = len(ds.auxiliary_data.MisfitWindows[model])
    misfits = []

    # collect the total misfit calculated by Pyadjoint
    for adjsrc in adjoint_sources.list():
        misfits.append(adjoint_sources[adjsrc].parameters["misfit_value"])
    
    summed_misfits = 0.5 * sum(misfits)/number_windows

    return summed_misfits


def misfit_stats(ds, model, include_lists=False):
    """
    Extract misfit statistics from a dataset for a given model.
    Return information as a dictionary object for easy access.    
 
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
    ds_windows = ds.auxiliary_data.MisfitWindows[model]
    
    syn_n, obs_n = _count_waveforms(ds.waveforms, model) 
    misfits, max_misfit, min_misfit = _misfit_info(ds_adjsrc) 
    cc_shift_secs, max_cc_values = _window_info(ds_windows) 
    
    # save stats in a dictionary
    # commented out sections are lists that I don't think are necessary 
    stats = {
        "number_stations": len(ds.waveforms.list()) ,
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
    from pyasdf import WaveformNotInFileException
    syn_n, obs_n = [],[]
    for sta in ds_waveforms.list():
        try:
            syn_n.append(len(ds_waveforms[sta]["synthetic_{}".format(model)]))
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
