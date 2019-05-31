"""
Convenience function for editing Pyasdf ASDFDataSet objects. All functions
work with the dataset as an input and act in-place on the dataset so no returns
"""


def clean_ds(ds):
    """
    removes synthetic waveforms and auxiliary data so that only observation
    data remains for new model runs
    :param ds:
    :return:
    """
    del_synthetic_waveforms(ds=ds, model=None)
    del_auxiliary_data(ds=ds)


def del_synthetic_waveforms(ds, model=None):
    """
    Removing only synthetic tagged waveforms from an asdf dataset.
    e.g. if a dataset has been used for processing and you want to retain
    the observation waveforms and start over with synthetics
    :param ds: pyasdf.ASDFDataSet
    """
    for sta in ds.waveforms.list():
        for stream in ds.waveforms[sta].list():
            if "synthetic" in stream:
                if (model is not None) and model in stream:
                    del ds.waveforms[sta][stream]
                elif model is None:
                    del ds.waveforms[sta][stream]


def del_auxiliary_data(ds):
    """
    delete all items in auxiliary data
    :param ds:
    :param model:
    :return:
    """
    for aux in ds.auxiliary_data.list():
        del ds.auxiliary_data[aux]


def del_adjoint_sources(ds, model=None):
    """
    delete adjoint sources from auxiliary data
    :param ds: pyasdf.ASDFDataSet
    :param model: str
    :return:
    """
    if model:
        for comp in ds.auxiliary_data.AdjointSources[model].list():
            del ds.auxiliary_data.AdjointSources[model][comp]
    else:
        for model in ds.auxiliary.data.AdjointSources.list():
            for comp in ds.auxiliary_data.AdjointSources[model].list():
                del ds.auxiliary_data.AdjointSources[model][comp]


def del_misfit_windows(ds, model=None):
    """
    delete misfit windows from auxiliary data
    :param ds: pyasdf.ASDFDataSet
    :param model: str
    :return:
    """
    if model:
        for comp in ds.auxiliary_data.MisfitWindows[model].list():
            del ds.auxiliary_data.MisfitWindows[model][comp]
    else:
        for model in ds.auxiliary.data.MisfitWindows.list():
            for comp in ds.auxiliary_data.MisfitWindows[model].list():
                del ds.auxiliary_data.MisfitWindows[model][comp]


def del_configs(ds, model=None):
    """
    delete config files from auxiliary data
    :param ds: pyasdf.ASDFDataSet
    :param model: str
    :return:
    """
    if model:
        for comp in ds.auxiliary_data.Configs[model].list():
            del ds.auxiliary_data.Configs[model][comp]
    else:
        for model in ds.auxiliary.data.Configs.list():
            for comp in ds.auxiliary_data.Configs[model].list():
                del ds.auxiliary_data.Configs[model][comp]


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

    :param ds:
    :param model_number:
    :param filepath:
    :return:
    """
    adjoint_sources = ds.auxiliary_data.AdjointSources[model]
    misfits = []
    for adjsrc in adjoint_sources.list():
        misfits.append(adjoint_sources[adjsrc].parameters["misfit_value"])
   
    summed_misfits = 0.5 * sum(misfits)/len(misfits)

    return summed_misfits

