"""
Here are some recipes for generating certain figures using the Inspector class
that didn't quite make it into the core functionality, but are still important.
Although they are specific to a problem the functions used may serve as a
template for accomplishing other tasks
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyatoa import Inspector
from pyatoa.utils.calculate import myround


def dlna_by_distance_bins():
    """
    Plot dlnA mean and standard deviation in bins organized by distance.
    """
    name = "birch"
    fid_insp = f"/Users/Chow/Documents/academic/vuw/forest/{name}/{name}"
    iteration = "i09"
    step_count = "s01"
    binsize_km = 25

    insp = Inspector(fid_insp)

    # Get the correct information into a separate dataframe
    df = insp.isolate(iteration, step_count)
    df = df.merge(insp.srcrcv, on=["event", "network", "station"])
    df = df[["dlnA", "distance_km"]]

    # Sort distances into orderly bins
    distance_bins = np.arange(
        0, myround(df.distance_km.max() + 1, base=binsize_km, choice="up"),
        binsize_km)
    df["dist_bin"] = pd.cut(df.distance_km, distance_bins)
    df.drop("distance_km", axis=1)

    # Find mean and std of dlnA based on bins
    mean_vals = df.groupby("dist_bin").mean()
    std_vals = df.groupby("dist_bin").std()

    # Get the center of the Categorical ranges that define the distance bins
    midpoints = pd.IntervalIndex(mean_vals.index).mid

    plt.errorbar(x=midpoints.to_numpy(), y=mean_vals.dlnA.to_numpy(),
                 yerr=std_vals.dlnA.to_numpy(), marker="o", markersize=4, c="r",
                 capsize=2)
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    yabsmax = max(abs(ymin), abs(ymax))

    # ax.set_ylim([-yabsmax, yabsmax])
    plt.xlabel("Distance [km]")
    plt.ylabel("dlnA")
    plt.title(f"{name} {iteration}{step_count}\n"
              f"Distance [{binsize_km}km bins] vs Amplitude Anomaly")
    plt.close()





