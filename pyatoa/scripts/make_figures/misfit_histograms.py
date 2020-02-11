"""
DEPRECATED in favor of pyatoa.plugins.seisflows.inspector.Inspector

Create misfit histograms to understand changes in model parameters.

Histograms available are:
1) "amplitude" showing the quantity dln(a) which is a
measure of amplitude differences between Observed and Synthetic
2) "cc_shift_in_seconds" showing the overall time shifts for all measured
windows in a dataset

Inspired by the histograms from Tape et al. 2010
"""
import os
import glob
import pyasdf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa.utils.asdf.extractions import histogram_data

# Set histogram look here
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['axes.linewidth'] = 2

# >>>> Set parameters here <<<<
path_to_datasets = "path/to/asdf_datasets"

# Choose the models you want to compare. Two is standard, more than two and the
# histogram becomes jumbled, but it is allowed. Colors correspond to the models,
# must include same number of colors as models to compare.
models_to_compare = ["m00", "m09"]
colors = ["darkorange", "deepskyblue"]

# Can choose "amplitude", "cc_shift_in_seconds" or both
data_types = ["amplitude", "cc_shift_in_seconds"]

output_figure_template = "misfithisto_{data_type}.png"
# >>>> Set parameters here <<<<


# Set the x-axis label for each of the data types
label_dict = {"cc_shift_in_seconds": "Time Shift (s)",
              "amplitude": "dlnA=ln(A_obs/A_syn)"}

# initiate the dictionary to hold values
histo_data = {}
for dtype in data_types:
    histo_data[dtype] = {}

# Read in the ASDF Dataset and collect the relevant information as per data_type
for fid in glob.glob(os.path.join(path_to_datasets, "*.h5")):
    print(fid)
    with pyasdf.ASDFDataSet(fid) as ds:
        for model in models_to_compare:
            for dtype in data_types:
                # initiate dictionary lists, fill them with data from dataset
                if model not in histo_data[dtype]:
                    histo_data[dtype][model] = []
                histo_data[dtype][model] += histogram_data(ds, model, dtype)

# Initiate the plot and construct the histograms
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
        
        # Plot the main histogram in full color
        n, bins, patches = plt.hist(
                                x=histo_data[data_type][model], 
                                bins=len(np.arange(-1*bound, bound, binsize)),
                                color=colors[i], histtype="bar", 
                                edgecolor="black", linewidth=2.5, label=model,
                                alpha=1 - i*0.2, zorder=10
                                )

        # Plot the overlying histogram that sits ontop of all the full histo,
        # this allows for visualization of histograms that overlap.
        if i < len(histo_data[data_type].keys()) - 1:
            n, bins, patches = plt.hist(
                                x=histo_data[data_type][model], 
                                bins=len(np.arange(-1*bound, bound, binsize)),
                                histtype="step", edgecolor=colors[i], linewidth=2., 
                                alpha=0.7, zorder=100
                                        )

    # dln(A) information only relevant in these bounds
    if data_type == "amplitude":
        plt.xlim([-1.25, 1.25])

    # Finalize plot details
    plt.xlabel(label_dict[data_type])
    plt.ylabel("Count")
    plt.title("Misfit Histogram")
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.grid(linewidth=1, linestyle=":", which="both", zorder=1)
    plt.axvline(x=0, ymin=0, ymax=1, linewidth=1.5, c="k", zorder=2, alpha=0.75,
                linestyle='--')
    plt.legend()
    plt.savefig(output_figure_template.format(data_type))
    plt.close()
