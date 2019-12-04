"""
Create misfit histograms to understand changes in model parameters
Inspired by the histograms from Tape et al. 2010
To do:
    This should be moved into Pyatoa's core functionalities at some point
"""
import os
import glob
import pyasdf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa.utils.asdf.extractions import histogram_data

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['axes.linewidth'] = 2

# set parameters
path_to_datasets = "./"
models_to_compare = ["m00", "m09"]
data_types = ["amplitude", "cc_shift_in_seconds"]
# data_types = ["amplitude"]
# data_types = ["cc_shift_in_seconds"]
colors = ["darkorange", "deepskyblue"]
label_dict = {"cc_shift_in_seconds": "Time Shift (s)",
              "amplitude": "dlnA=ln(A_obs/A_syn)"}

# initiate the dictionary to hold values
histo_data = {}
for dtype in data_types:
    histo_data[dtype] = {}

# initiate the plot
for fid in glob.glob(os.path.join(path_to_datasets, "*.h5")):
    print(fid)
    with pyasdf.ASDFDataSet(fid) as ds:
        for model in models_to_compare:
            for dtype in data_types:
                # initiate dictionary lists, fill them with data from dataset
                if model not in histo_data[dtype]:
                    histo_data[dtype][model] = []
                histo_data[dtype][model] += histogram_data(ds, model, dtype)

for data_type in histo_data.keys():
    f, ax = plt.subplots(1)
    for i, model in enumerate(histo_data[data_type].keys()):
        min_value = np.floor(min(histo_data[data_type][model]))
        max_value = np.ceil(max(histo_data[data_type][model]))
        bound = max(abs(min_value), abs(max_value))
        if data_type == "amplitude":
            binsize = 0.2
        else:
            binsize = 0.25
        
        # plot the main histogram in full color
        n, bins, patches = plt.hist(
                                x=histo_data[data_type][model], 
                                bins=len(np.arange(-1*bound, bound, binsize)),
                                color=colors[i], histtype="bar", 
                                edgecolor="black", linewidth=2.5, label=model,
                                alpha=1 - i*0.2, zorder=10
                                )

        # plot the overlying histogram that sits ontop of all the full histo
        if i < len(histo_data[data_type].keys()) - 1:
            n, bins, patches = plt.hist(
                                x=histo_data[data_type][model], 
                                bins=len(np.arange(-1*bound, bound, binsize)),
                                histtype="step", edgecolor=colors[i], linewidth=2., 
                                alpha=0.7, zorder=100
                                        )

    # amplitude information only relevant in these bounds 
    if data_type == "amplitude":
        plt.xlim([-1.25, 1.25])
    else:
        plt.xlim([-6, 6])

    plt.xlabel(label_dict[data_type])
    plt.ylabel("Count")
    plt.title("Misfit Histogram")
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.grid(linewidth=1, linestyle=":", which="both", zorder=1)
    plt.axvline(x=0, ymin=0, ymax=1, linewidth=1.5, c="k", zorder=2, alpha=0.75,
                linestyle='--')
    plt.legend()
    plt.savefig(f'misfithisto_{data_type}.png')
    plt.close()
             
        
        
        
    
