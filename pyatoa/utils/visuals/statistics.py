"""
Plots of statistical information for use in misfit analysis
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from pyatoa.utils.tools.calculate import myround

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def plot_output_optim(path_to_optim, show=False, save=''):
    """
    Sesiflows specific
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
    from pyatoa.utils.tools.io import parse_output_optim

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
        ax.annotate(s="iter {}".format(itr), xytext=(offset+0.2, misfits_[0]),
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


def plot_misfit_histogram(misfit_values, config, binsize=0.1):
    """
    Make histograms of misfit values for model m_a, can compare with model m_b
    :param misfit_values: dict
    :param binsize:
    :return:
    """
    for model in misfit_values.keys():
        misfits = np.fromiter(misfit_values[model].values(), dtype="float")
        maxmisfit = myround(misfits.max(), base=1, choice="up")
        n, bins, patches = plt.hist(x=misfits,
                                    bins=len(np.arange(0, maxmisfit, binsize)),
                                    range=(0, maxmisfit), color="orange",
                                    histtype="bar", edgecolor="black",
                                    linewidth=1.5, zorder=10, label=model
                                    )
    plt.xlabel("Misfit Value ")
    plt.ylabel("Count (N={})".format(len(misfits)))
    plt.title("{eid} Misfits ".format(eid=config.event_id))
    plt.grid(linewidth=1.0, which='both', zorder=1)
    plt.legend()
    plt.xlim([-0.05, bins.max() + 0.05])
    plt.ylim([0, max(n) + 0.5])
    plt.show()


def plot_cc_time_shift_histogram(cc_time_shifts, config, binsize=0.1):
    """
    create a histogram of cross correlation time shifts
    """
    raise NotImplementedError

    for model in misfit_values.keys():
        ccts = np.fromiter(cc_time_shifts[model].values(), dtype="float")
        max_ccts = myround(max([abs(ccts.max()), abs(ccts.min())]), base=0.1,
                           choice="up"
                           )
        n, bins, patches = plt.hist(x=ccts,
                                    bins=len(np.arange(0, max_ccts, binsize)),
                                    range=(-max_ccts, max_ccts), color="orange",
                                    histtype="bar", edgecolor="black",
                                    linewidth=1.5, zorder=10, label=model
                                    )
        mean = np.mean(ccts)
        onesigma = np.std(ccts)
        mean_one_sigma = "$\mu$ + $\sigma$ = {m:.2f} $\pm$ {s:.2f}s".format(
            m=mean, s=onesigma
        )
        plt.text(x=-max_ccts + .5, y=n.max() - 1, s=mean_one_sigma,
                 bbox=dict(facecolor='w', alpha=0.5), zorder=11
                 )
        plt.axvline(x=mean, ymin=n.min(), ymax=n.max(), linestyle='-',
                    color='k')
        for i in [-1, 1]:
            plt.axvline(x=i * onesigma, ymin=n.min(), ymax=n.max(),
                        linestyle='--', color='k'
                        )
    # plot attributes
    plt.xlabel("CC Time Shift [s]")
    plt.ylabel("Count (N={})".format(len(ccts)))
    plt.title("{eid} CC Time Shift ".format(eid=config.event_id))
    plt.legend()
    plt.grid(linewidth=1.0, which='both', zorder=1)
    plt.xlim([-(max_ccts + 0.05), max_ccts + 0.05])
    plt.ylim([0, max(n) + 0.5])
    plt.show()
