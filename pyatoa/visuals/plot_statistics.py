"""
Plots of statistical information for use in misfit analysis
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyatoa.utils.operations.calculations import myround
from pyatoa.visuals.plot_utils import align_yaxis, pretty_grids, format_axis

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


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
