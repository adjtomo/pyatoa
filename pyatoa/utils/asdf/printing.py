"""
Pretty printing certain parts of an asdf dataset
"""


def print_statistics(ds, model):
    """
    print the information from statistics in a neat organized fashion
    """
    try:
        stats = ds.auxiliary_data.Statistics[model].parameters
    except AttributeError:
        return
    
    format_str = "{sta:>5s} {obs:>2d} {syn:>2s} "
    for i, sta in enumerate(stats["stations"]):
        print("hello")

    raise NotImplementedError
              
     
