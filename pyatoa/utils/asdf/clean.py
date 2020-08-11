"""
Convenience functions for removing data from Pyasdf ASDFDataSet objects. 
All functions work with the dataset as an input and act in-place on the 
dataset so no returns
"""
from pyatoa.utils.form import format_iter, format_step


def clean_dataset(ds, iteration=None, step_count=None, fix_windows=False):
    """
    Removes synthetic waveforms and auxiliary data so that only observation
    data remains for new iterations. If no iteration is given, will 
    wipe all non-observation data and all auxiliary data
    
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset to be cleaned
    :type fix_windows: bool
    :param fix_windows: don't delete MisfitWindows
    :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
    :type step_count: str or int
    :param step_count: step count e.g. "s00". Will be formatted so int ok.
    """
    del_synthetic_waveforms(ds=ds, iteration=iteration, step_count=step_count)

    # Retain misfit windows if windows to  be fixed
    if fix_windows:
        retain = ["MisfitWindows"]
    else:
        retain = []

    del_auxiliary_data(ds=ds, iteration=iteration, step_count=step_count, 
                       retain=retain)


def del_synthetic_waveforms(ds, iteration=None, step_count=None):
    """
    Remove "synthetic_{iter_tag}{step_tag}" tagged waveforms from an asdf 
    dataset. If no iter_tag number given, wipes all synthetic data from dataset.   
 
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset to be cleaned
    :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
    :type step_count: str or int
    :param step_count: step count e.g. "s00". Will be formatted so int ok.
    """
    iter_tag = format_iter(iteration)
    step_tag = format_step(step_count)

    for sta in ds.waveforms.list():
        for stream in ds.waveforms[sta].list():
            # stream is e.g. 'synthetic_i00s00'
            if "synthetic" in stream:
                if (iter_tag is not None) and iter_tag in stream:
                    if (step_tag is not None) and step_tag in stream:
                        del ds.waveforms[sta][stream]
                    elif step_tag is None:
                        del ds.waveforms[sta][stream]
                elif iter_tag is None:
                    del ds.waveforms[sta][stream]


def del_auxiliary_data(ds, iteration=None, step_count=None, retain=None,
                       only=None):
    """
    Delete all items in auxiliary data for a given iter_tag, if iter_tag not
    given, wipes all auxiliary data.

    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset to be cleaned
    :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
    :type step_count: str or int
    :param step_count: step count e.g. "s00". Will be formatted so int ok.
    :type retain: list of str
    :param retain: list of auxiliary data tags to retain, that is: delete all 
        auxiliary data EXCEPT FOR the names given in this variable
    :type only: list of str
    :param only: list of auxiliary data tags to remove, that is: ONLY delete 
        auxiliary data that matches the names given in this variable. 
        Lower in priority than 'retain'
    """
    iter_tag = format_iter(iteration)
    step_tag = format_step(step_count)
    
    for aux in ds.auxiliary_data.list():
        # Check if auxiliary data tag matches optional lists
        if retain and aux in retain:
            continue
        if only and aux not in only:
            continue

        if (iter_tag is not None) and hasattr(ds.auxiliary_data[aux], iter_tag):
            # If the aux data doesn't contain this iter_tag, nothing to clean
            if (step_tag is not None) and (
                        hasattr(ds.auxiliary_data[aux][iter_tag], step_tag)):
                del ds.auxiliary_data[aux][iter_tag][step_tag]
            # If no 'step_tag' given, simply delete the iter_tag subgroup 
            elif step_tag is None:
                del ds.auxiliary_data[aux][iter_tag]
        elif iter_tag is None:
            del ds.auxiliary_data[aux]

