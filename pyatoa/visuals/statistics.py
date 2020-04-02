"""
Plots of statistical information for use in misfit analysis
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def plot_output_optim(path_to_optim, show=False, save=''):
    """
    Sesiflows specific function

    Read the text file output.optim, which contains
    misfit values, step length values and iteration numbers.
    Parse the values and plot them on a scatter plot to show how the misfit
    evolves with subsequent iterations. Put some auxiliary plot information on
    for a more informative overall plot

    :type path_to_optim: str
    :param path_to_optim: path to the 'output.optim' text file
    :type show: bool
    :param show: show the plot
    :type save: str
    :param save: output filename to save the figure
    """
    from pyatoa.utils.write import parse_output_optim

    # read in the values of the file
    iterations, steplens, misfits = parse_output_optim(path_to_optim)

    # Plot the misfit values, and steplengths on a twin axis
    f, ax = plt.subplots(1)
    twax = ax.twinx()

    # Set the normalized colormap to differentiate iterations
    colormap = cm.nipy_spectral
    norm = mpl.colors.Normalize(vmin=min(iterations), vmax=max(iterations))

    # Offset allows for cumulative plotting despite restarting scatter plot
    offset = 0
    for itr in np.unique(iterations):
        # Determine unique arrays for each iteration
        indices = np.where(iterations == itr)
        misfits_ = misfits[indices]
        steplens_ = steplens[indices]

        # Set a specific color for each iteration and plot as a scatterplot
        color = [colormap(norm(itr))] * len(misfits_)

        # Plot the misfit values
        ax.scatter(range(offset, offset + len(misfits_)), misfits_, c=color,
                   zorder=5, s=60)
        ax.annotate(s=f"iter {itr}", xytext=(offset+0.2, misfits_[0]),
                    xy=(offset, misfits_[0]), color=colormap(norm(itr)),
                    fontsize=9, zorder=6
                    )

        # Plot the steplengths, remove values of 0 steplength
        steplens_[steplens_ == 0] = 'nan'
        twax.scatter(range(offset, offset + len(steplens_)), steplens_, c=color,
                     zorder=3, marker='d', alpha=0.75, s=50, edgecolor='w')

        # Connect scatter plots with a line in the back
        ax.plot(range(offset, offset + len(misfits_)), misfits_,
                c=colormap(norm(itr)))

        # Make sure scatter points overlap because the starting misfit of the
        # next iteration is the same misfit as the last iteration
        offset += len(misfits_) - 1

    # Scatter points outside the plot for legend
    ax.scatter(-1, min(misfits), c='w', edgecolor='k', marker='o',
               label='misfits')
    ax.scatter(-1, min(misfits), c='w', edgecolor='k', marker='d',
               label='step lengths')

    # Put a line at y=1 for step lengths to show where L-BFGS is working
    twax.axhline(y=1, xmin=0, xmax=1, linestyle='--', linewidth=1.5, c='k',
                 alpha=0.5, zorder=2)

    # Because the scatter plots overlap, number of points must be defined 
    num_points = len(iterations) - len(np.unique(iterations))
    steplens_min_notzero = steplens[np.where(steplens > 0)].min()
    # Set plot attributes
    if max(misfits)/min(misfits) > 10:
        ax.set_yscale('log')
    twax.set_yscale('log')
    twax.set_ylim(bottom=steplens_min_notzero)
    plt.title("seisflows output.optim")
    ax.set_xlim([-0.5, num_points + 0.5])
    ax.set_xlabel("step count")
    ax.set_ylabel("pyatoa misfits")
    ax.legend()
    twax.set_ylabel("step lengths", rotation=270, labelpad=15)

    # Set the tick labelling based on the number of iterations
    plt.xticks(range(0, num_points, num_points//10 or 1))

    # Set the format of the plot to match pretty_grids formatting
    for axis in [ax, twax]:
        axis.tick_params(which='both', direction='in', top=True, right=True)
    for axis_ in ['major', 'minor']:
        ax.grid(which=axis_, linestyle=':', linewidth='0.5',
                color='k', alpha=0.25)

    if show:
        plt.show()
    if save:
        plt.savefig(save)


def plot_misfit_histogram(path_to_datasets, choice, show=False, save=False):
    """
    Make histograms of misfit values for model m_a, can compare with model m_b

    :type path_to_datasets: str
    :param path_to_datasets: path to the pyasdf datasets
    :type choice: str
    :param choice: choice for misfit value, available
        "cc_shift_in_seconds", "amplitude", "both"
    :type show: bool
    :param show: show the plot after making it
    :type save: str
    :param save: fid to save the figure
    """
    import os
    import pyasdf
    from glob import glob
    from pyatoa.utils.asdf.extractions import histogram_data

    # Pre set available data types
    data_types = ["cc_shift_in_seconds", "amplitude", "both"]
    colors = ["darkorange", "deepskyblue"]
    assert(choice in data_types), f"choice must be in {data_types}"

    # Pre set the x-axis label for each of the data types
    label_dict = {"cc_shift_in_seconds": "Time Shift (s)",
                  "amplitude": "dlnA=ln(A_obs/A_syn)"}

    # initiate the dictionary to hold values
    histo_data = {}
    for dtype in data_types:
        histo_data[dtype] = {}

    # Read in the ASDF Dataset and collect relevant information per data_type
    for fid in glob(os.path.join(path_to_datasets, "*.h5")):
        with pyasdf.ASDFDataSet(fid) as ds:
            # Models comprise the first and final models
            models = [ds.auxiliary_data.Statistics.list()[0],
                      ds.auxiliary_data.Statistics.list()[-1]
                      ]
            for model in models:
                for dtype in data_types:
                    # initiate dictionary lists, fill them with data from ds
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
            plt.hist(x=histo_data[data_type][model],
                     bins=len(np.arange(-1*bound, bound, binsize)),
                     color=colors[i], histtype="bar",
                     edgecolor="black", linewidth=2.5, label=model,
                     alpha=1 - i * 0.2, zorder=10
                     )

            # Plot the overlying histogram that sits ontop of all the full histo
            # this allows for visualization of histograms that overlap.
            if i < len(histo_data[data_type].keys()) - 1:
                plt.hist(x=histo_data[data_type][model],
                         bins=len(np.arange(-1*bound, bound, binsize)),
                         histtype="step", edgecolor=colors[i], linewidth=2.,
                         alpha=0.7, zorder=100
                         )

        # dln(A) information only relevant between -1 and 1
        if data_type == "amplitude":
            plt.xlim([-1.25, 1.25])

        # Finalize plot details
        plt.xlabel(label_dict[data_type])
        plt.ylabel("Count")
        plt.title("Misfit Histogram")
        plt.tick_params(which='both', direction='in', top=True, right=True)
        plt.grid(linewidth=1, linestyle=":", which="both", zorder=1)
        plt.axvline(x=0, ymin=0, ymax=1, linewidth=1.5, c="k", zorder=2,
                    alpha=0.75, linestyle='--')
        plt.legend()

        if save:
            plt.savefig(save)
        if show:
            plt.show()
        plt.close()



