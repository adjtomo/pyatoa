"""
Simple scatterplot to show misfit convergence based on the misfits.json
output that is created during a Seisflows run
"""
import os
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 14
mpl.rcParams['lines.linewidth'] = 2.
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['axes.linewidth'] = 3


def parse_json(fid):
    """
    Parse the misfits json file

    :type fid: str
    :param fid: File identifier for misfits.json file
    """
    misfits = json.load(open(fid))

    misfit_array, window_array, model = [], [], []
    for i, key in enumerate(misfits):
        try:
            misfit = misfits[key]['s00']['misfit']
            windows = misfits[key]['s00']['windows']
            # if s00 is available but is the same as the previous misfit,
            # force the exception to get s01
            if misfit_array and misfit == misfit_array[-1]:
                raise KeyError
        # if thrifty inversion, s00 is taken from the previous iteration
        except KeyError:
            misfit = misfits[key]['s01']['misfit']
            windows = misfits[key]['s01']['windows']
        
        misfit_array.append(misfit)
        window_array.append(windows)
        model.append(i)
    
    return model, np.array(misfit_array), window_array


if __name__ == "__main__":
    # Set pathanames here
    plot_windows = True
    basepath = "./"
    fid_list = [os.path.join(basepath, "30cc/pyatoa.io/misfits.json"),
                os.path.join(basepath, "30mtm/pyatoa.io/misfits.json"),
                # os.path.join(basepath, "both/pyatoa.io/misfits.json")
                ]

    # Colors and labels for the various files, must match length fid_list
    color_list = ["darkorange", "mediumturquoise"]
    label_list = ["traveltime_cc", "multitaper misfit", "both"]

    # set up plots
    f, ax1 = plt.subplots(figsize=(8,6))
    if plot_windows:
        ax2 = ax1.twinx()
    else:
        ax2 = None

    # Plot the scatterplot
    for i, fid in enumerate(fid_list):
        models, misfits, windows = parse_json(fid)
        misfits /= misfits.max() 
        ax1.plot(models, misfits, 'o-', c=color_list[i], label=label_list[i])
        if ax2:
            ax2.plot(models, windows, 'd-.', c=color_list[i], alpha=0.5)

    # Finalize plotting attributes
    ax1.legend(loc="center right")
    plt.title('Checkerboard Convergence')
    ax1.set_xticks(models)
    ax1.set_xlabel('Model Number')
    ax1.set_ylabel('Total Misfit (solid)')
    if ax2:
        ax2.set_ylabel('Number of Windows (dashed)', rotation=270, labelpad=15.)
    ax1.grid(True, alpha=0.5, linestyle='--', linewidth=1.)
    
    plt.savefig('convergence.png') 
    plt.show()


