"""
Plots of statistical information for use in misfit analysis
"""
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from pyatoa.utils.operations.calculations import myround

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def read_and_plot_output_optim(path_to_optim, show=False, save=False):
    """
    Sesiflows specific - Reading the text file output.optim, which contains
    misfit values
    :param path_to_optim:
    :return:
    """
    with open(path_to_optim, 'r') as f:
        lines = f.readlines()

    # Parse the file, skip the header, ignore tail
    iterations, steplengths, misfits = [], [], []
    for line in lines[2:]:
        line = line.strip().split()
        # Each iteration will have an integer to represent the iter number
        if len(line) == 3:
            iteration = line[0]
            iterations.append(int(iteration))
            steplengths.append(float(line[1]))
            misfits.append(float(line[2]))
        # Each trial step length will follow and will carry the same iteration
        elif len(line) == 2:
            iterations.append(int(iteration))
            steplengths.append(float(line[0]))
            misfits.append(float(line[1]))
        elif len(line) == 0:
            continue
        else:
            print(line)
            print("invalid line length encountered in output.optim")
            sys.exit()

    # Set the lists as numpy arrays to do some manipulations
    iterations = np.asarray(iterations)
    steplengths = np.asarray(steplengths)
    misfits = np.asarray(misfits)

    # Plot
    offset = 0
    norm = mpl.colors.Normalize(vmin=min(iterations), vmax=max(iterations))
    for itr in np.unique(iterations):
        # Determine unique arrays for each iteration
        indices = np.where(iterations == itr)
        misfits_ = misfits[indices]
        colormap = cm.nipy_spectral

        # Set a specific color for each iteration and plot as a scatterplot
        color = [colormap(norm(itr))] * len(misfits_)
        plt.scatter(range(offset, offset + len(misfits_)), misfits_, c=color,
                    label="Iter {}".format(str(itr)), zorder=5)
        plt.annotate(s="Iter {}".format(itr), xytext=(offset+0.5, misfits_[0]),
                     xy=(offset, misfits_[0]), c=colormap(norm(itr)),
                     fontsize=9, zorder=6
                     )
        plt.plot(range(offset, offset + len(misfits_)), misfits_,
                 c=colormap(norm(itr)))
        offset += len(misfits_) - 1

    # Set plot attributes
    if max(misfits)/min(misfits) > 10:
        plt.yscale('log')
    plt.title("Seisflows Output Optimization")
    plt.xlabel("Step Count")
    plt.ylabel("Pyatoa Misfits")
    # plt.legend(ncol=len(np.unique(iterations))//5)
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.xticks(range(0, offset + 1, 2))
    for axis_ in ['major', 'minor']:
        plt.grid(which=axis_, linestyle=':', linewidth='0.5',
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
