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


