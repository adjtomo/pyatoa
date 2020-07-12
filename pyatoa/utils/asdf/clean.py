"""
Convenience functions for removing data from Pyasdf ASDFDataSet objects. 
All functions work with the dataset as an input and act in-place on the 
dataset so no returns
"""


def clean_ds(ds, iteration=None, step_count=None, fix_windows=False):
    """
    Removes synthetic waveforms and auxiliary data so that only observation
    data remains for new iteration runs. If no iteration nmber is given, will wipe
    all non-observation data and all auxiliary data
    
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset to be cleaned
    :type fix_windows: bool
    :param fix_windows: don't delete MisfitWindows
    :type iteration: str
    :param iteration: iteration number, e.g. "m00"
        :type step_count: str
    :param step_count: step_count number, e.g. "s00"
    """
    del_synthetic_waveforms(ds=ds, iteration=iteration, step_count=step_count)
    retain = []
    if fix_windows:
        retain.append("MisfitWindows")

    del_auxiliary_data(ds=ds, iteration=iteration, step_count=step_count, 
                       retain=retain)


def del_synthetic_waveforms(ds, iteration=None, step_count=None):
    """
    Remove "synthetic_{iteration}" tagged waveforms from an asdf dataset.
    If no iteration number given, wipes all synthetic data from dataset.   
 
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset to be cleaned
    :type iteration: str
    :param iteration: iteration number, e.g. "m00"
    :type step_count: str
    :param iteration: step_count count, e.g. "s00"
    """
    for sta in ds.waveforms.list():
        for stream in ds.waveforms[sta].list():
            if "synthetic" in stream:
                if (iteration is not None) and iteration in stream:
                    if (step_count is not None) and step_count in stream:
                        del ds.waveforms[sta][stream]
                    elif step_count is None:
                        del ds.waveforms[sta][stream]
                elif iteration is None:
                    del ds.waveforms[sta][stream]


def del_auxiliary_data(ds, iteration=None, step_count=None, retain=None,
                       only=None):
    """
    Delete all items in auxiliary data for a given iteration, if iteration not given,
    wipes all auxiliary data.

    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset to be cleaned
    :type iteration: str
    :param iteration: iteration number, e.g. "m00"
    :type step_count: str
    :param step_count: step_count number, e.g. "s00"
    :type retain: list of str
    :param retain: list of auxiliary data tags to retain, that is: delete all 
        auxiliary data EXCEPT FOR the names given in this variable
    :type only: list of str
    :param only: list of auxiliary data tags to remove, that is: ONLY delete 
        auxiliary data that matches the names given in this variable. 
        Lower in priority than 'retain'
    """
    for aux in ds.auxiliary_data.list():
        if retain and aux in retain:
            continue
        if only and aux not in only:
            continue
        if (iteration is not None) and hasattr(ds.auxiliary_data[aux], 
                                               iteration):
            # If the aux data doesn't contain this iteration, nothing to clean
            if (step_count is not None) and (
                        hasattr(ds.auxiliary_data[aux][iteration], step_count)):
                del ds.auxiliary_data[aux][iteration][step_count]
            # If no 'step_count' given, simply delete the iteration subgroup 
            elif step_count is None:
                del ds.auxiliary_data[aux][iteration]
        elif iteration is None:
            del ds.auxiliary_data[aux]
            
