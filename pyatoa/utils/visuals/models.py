"""
Plot various model visualizations using Numpy arrays
"""
import numpy as np


class Model:
    """
    A class for dealing with model plotting
    """
    def __init__():


    # Global dictionary that sets the indices of the array with coordinates
    md = {"x": 0, "y": 1, "z": 2}


    def clip(model, bound, normal="z", verbose=True):
        """
        Clip the model, retaining only the values from value, in the direction of 
        the normal
        """
        idx = md[normal[-1]]
        # negative normal means we face the other direction
        if "-" in normal:
            indices = np.where(model[:, idx] <= bound)[0]
        else:
            indices = np.where(model[:, idx] >= bound)[0]

        if verbose:
            print(f"clipping model at {bound} for {normal}: "
                  f"{len(model)} -> {len(indices)} points")

        return model[indices]



    def surface(model):
        """
        Plot surface projection of a model
        """

