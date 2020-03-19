#!/usr/bin/env python3
"""
Reads, parses and plots the misfits.json file outputted by a Seisflows/Pyatoa
workflow, in order to show convergence behavior. Allows comparisons of multiple
inversions. This functionality is also provided in the Inspector class.
"""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 14
mpl.rcParams['lines.linewidth'] = 2.
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['axes.linewidth'] = 3


def parse_json(fid, choice="models", value="windows"):
    """
    Parse the misfits json file

    :type fid: str
    :param fid: File identifier for misfits.json file
    :type choice: str
    :param choice: 'models' to only collect information for each model update
                   'steps' to collect information about each step
    :type value: str
    :param value: choice of value to plot from the json file
        "windows" or ,"adj_srcs"
    """
    misfits = json.load(open(fid))

    misfit_array, window_array, models, iters = [], [], [], []
    if choice == "models":
        for i, model in enumerate(misfits):
            print(model)
            # Misfit defined by M00S00
            try:
                misfit = misfits[model]["s00"]["misfit"]
                windows = misfits[model]["s00"][value]
                print(f"\ts00: {misfit:.2E}")
            # If thrifty inversion, S00 doesn't exist and misfit needs to be 
            # retrieved from the minimum misfit in previous line search
            except KeyError:
                smallest_misfit, step_check = 1E9, 0
                for step in list(misfits[previous_model].keys()):
                    temp_misfit = misfits[previous_model][step]["misfit"] 
                    if temp_misfit < smallest_misfit:
                        misfit = misfits[previous_model][step]["misfit"] 
                        windows = misfits[previous_model][step][value]
                        step_check = step
                    else:
                        smallest_misfit = temp_misfit
                print(f"\t{previous_model}{step_check}: {misfit:.2E}")
            
            misfit_array.append(misfit)
            window_array.append(windows)
            previous_model = model
            models.append(i)
    elif choice == "steps":
        k = 0
        for i, model in enumerate(misfits):
            iters.append(k)
            for j, step in enumerate(misfits[model]):
                misfit_array.append(misfits[model][step]["misfit"])
                window_array.append(misfits[model][step][value])
                models.append(k)
                k += 1
    
    return models, np.array(misfit_array), window_array, iters


if __name__ == "__main__":
    # vv SET PARAMETERS HERE vv
    # Figure parameters
    choice = "models"  # 'models' or 'iter'
    plot_windows = True  # plot number of misfit windows alongside
    basepath = "./"  # if common base location
    figure_size = (8, 6)
    title = ""  # Figure title
    fid_out = "convergence.png"  # name of the output file
    show = True  # show the figure

    # Choice of files to plot, lists must have the same lengths
    fid_list = ["path/to/misfits.json_1", "path/to/misfits.json_2"]
    color_list = ["red", "black"]
    label_list = ["label_1", "label_2"]
    # ^^ SET PARAMETERS HERE ^^

    # Set up plots
    f, ax1 = plt.subplots(figsize=figure_size)
    if plot_windows:
        ax2 = ax1.twinx()
    else:
        ax2 = None

    # Plot the scatterplot
    for i, fid in enumerate(fid_list):
        print(fid)
        models, misfits, windows, iters = parse_json(fid, choice, "windows")
        misfits /= misfits.max() 
        ax1.plot(models, misfits, 'o-', c=color_list[i], label=label_list[i])
        if ax2:
            ax2.plot(models, windows, 'v--', c=color_list[i])

        # For step count plot iters
        if iters:
            for j in iters:
                plt.axvline(x=j, linewidth=3, c=color_list[i], linestyle=":", 
                            alpha=0.25)

    # Finalize plotting attributes
    ax1.legend(loc="center right")
    if title:
        plt.title(title)
    if choice == "models":
        ax1.set_xlabel('Model Number')
    elif choice == "steps":
        ax1.set_xlabel('Step Count')
    ax1.set_xticks(models)
    ax1.set_ylabel('total normalized misfit (solid)')
    if ax2:
        ax2.set_ylabel('number of measurements (dashed)', rotation=270, 
                       labelpad=15.)
    ax1.grid(True, alpha=0.5, linestyle='--', linewidth=1.)

    if fid_out:
        plt.savefig(fid_out)
    if show:
        plt.show()
