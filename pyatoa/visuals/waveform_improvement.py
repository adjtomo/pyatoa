#!/usr/bin/env python3
"""
Functions to create a figure showing progressive waveform improvement over
the course of a seismic inversion.

Show the changes in synthetic waveforms with progressive model updates. 
Each individual model gets its on row in the plot.
"""
import os
import sys
import glob
import pyasdf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
  
from pyatoa import Config, Manager
from pyatoa.utils.form import format_event_name
from pyatoa.utils.asdf.fetch import windows_from_dataset


def setup_plot(nrows, ncols, label_units=False):
    """
    Dynamically set up plots according to number_of given
    Returns a list of lists of axes objects
    e.g. axes[i][j] gives the ith column and the jth row

    :type nrows: int
    :param nrows: number of rows in the gridspec
    :type ncols: int
    :param ncols: number of columns in the gridspec
    :type label_units: bool
    :param label_units: label the y-axis units and tick marks. can get messy
        if there are multiple models plotted together, so usually best to leave
        it off.
    :rtype axes: matplotlib axes
    :return axes: axis objects
    """
    gs = mpl.gridspec.GridSpec(nrows, ncols, hspace=0, wspace=0.05, 
                               width_ratios=[1] * ncols,
                               height_ratios=[3] * nrows
                               )

    axes = []
    for row in range(0, gs.get_geometry()[0]):
        components = []
        for col in range(0, gs.get_geometry()[1]):
            if row == 0:
                ax = plt.subplot(gs[row, col])
            else:
                # Share the x axis with the same component in the column
                ax = plt.subplot(gs[row, col], sharex=axes[0][col],
                                 sharey=axes[0][col]
                                 )
            ax.set_axisbelow(True)
            ax.tick_params(which='both', direction='in', top=True,
                                 right=True)
            # Set the grids on
            ax.minorticks_on()
            components.append(ax)
        axes.append(components)

    # remove x-tick labels except for last axis
    for row in axes[:-1]:
        for col in row:
            plt.setp(col.get_xticklabels(), visible=False)

    # remove y-tick labels except for first axis
    if label_units:
        for i, row in enumerate(axes):
            for j, col in enumerate(row):
                if j != 0: 
                    plt.setp(col.get_yticklabels(), visible=False)
                else:    
                    col.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    else:
        for i, row in enumerate(axes):
            for col in row:
                plt.setp(col.get_yticklabels(), visible=False)

    return axes


def collect_data(ds, sta, synthetics_only=False):
    """
    Parse dataset for given station, gather observed and synthetic data, 
    preprocess data and return as stream objects.

    :type ds: pyasdf.ASDFDataSet
    :param dsfid: dataset containing waveform information
    :type synthetics_only: bool
    :param synthetics_only: synthetics only parameter to be passed to the
        Config. if False, instrument response will be removed from obs.
    """
    network, station = sta.split(".")
    models = ds.auxiliary_data.MisfitWindows.list()

    # Preprocess all traces using Pyatoa and store in dict
    st_obs, synthetics, windows = None, {}, {}
    for m, model in enumerate(models):
        steps = ds.auxiliary_data.MisfitWindows[model].list()
        for s, step in enumerate(steps):
            # If specific paths are given, only choose those.
            # Else only take m00s00, and last step each other model
            if m == 0 and s != 0:
                continue
            elif m != 0 and s != len(steps) - 1:
                continue
            path = f"{model}/{step}"
            print(path)
            cfg = Config(path=path, ds=ds, synthetics_only=synthetics_only)
            mgmt = Manager(config=cfg, ds=ds)   
            mgmt.load(sta)
            mgmt.standardize()
            mgmt.preprocess() 
            synthetics[path] = mgmt.st_syn.copy() 
            windows[path] = windows_from_dataset(ds=ds, model=model, 
                                                 step=step, net=network, 
                                                 sta=station)
            # Observed waveform will be the same
            if st_obs is None:
                st_obs = mgmt.st_obs.copy()

    return st_obs, synthetics, windows, mgmt.stats.time_offset_sec


def wavimprov(ds, sta, min_period, max_period, save=None, 
              synthetics_only=False, trace_length=None, show=True, **kwargs):

    """
    Plot waveforms iterative based on model updates

    :type ds: pyasdf.ASDFDataSet
    :param dsfid: dataset containing waveform information
    :type min_period: float
    :param min_period: minimum filter period for waveforms
    :type max_period: float
    :param max_period: maximum filter period for waveforms
    :type synthetics_only: bool
    :param synthetics_only: synthetics only parameter to be passed to the
        Config. if False, instrument response will be removed from obs.
    :type trace_length: list of floats
    :param trace_length: [trace_start, trace_end] will be used to set the x
        limit on the waveform data. If none, no xlim will be set
    :type show: bool
    :param show: Show the plot or do not
    """
    # Check that plotting can be performed
    if sta not in ds.waveforms.list():
        raise KeyError(f"Dataset has no station {sta}")

    # Plotting parameters
    dpi = kwargs.get("dpi", 100)
    color_list = kwargs.get("color_list", ["r", "b", "g"])
    figsize = kwargs.get("figsize", (10, 8))

    # Set some random info for plotting.
    z = 5  # for zorder

    st_obs, synthetics, windows, time_offset = collect_data(
                    ds, sta, synthetics_only=synthetics_only)

    # Instantiate the plotting object
    f = plt.figure()
    axes = setup_plot(nrows=len(synthetics.keys()), ncols=len(st_obs))
    middle_column = len(st_obs) // 2
    t = st_obs[0].times(
           reftime=st_obs[0].stats.starttime - time_offset)

    # Make sure the models are in order 
    synthetic_keys = list(synthetics.keys())
    synthetic_keys.sort()
    synthetic_init = synthetics[synthetic_keys[0]]
    component_list = [_.stats.component for _ in st_obs]

    # Plot each model on a different row
    ymin = ymax = None
    for row, syn_key in enumerate(synthetic_keys):
        ylab = syn_key.split('_')[-1]  # e.g. 'm00'

        # Plot each component in a different column
        for col, comp in enumerate(component_list):
            obs = st_obs.select(component=comp)[0]
            syn = synthetics[syn_key].select(component=comp)[0]

            # Plot waveforms
            a1, = axes[row][col].plot(t, obs.data, 'k', zorder=z,
                                      label="True", linewidth=3)
            a2, = axes[row][col].plot(t, syn.data, 
                                      color_list[col], zorder=z,
                                      label="Syn", linewidth=3)

            # set a 'global' min and max value for the y-axis

            # Plot windows if available for this given component
            if comp in windows[syn_key]:
                for w, win in enumerate(windows[syn_key][comp]):
                    ymin, ymax = axes[row][col].get_ylim()

                    tleft = win.left * win.dt + time_offset
                    tright = win.right * win.dt + time_offset
                    tshift = win.cc_shift * syn.stats.delta

                    axes[row][col].add_patch(mpl.patches.Rectangle(
                                  xy=(tleft, ymin * 2), 
                                  width=tright-tleft, edgecolor='k',
                                  height=2*(ymax + abs(ymin)), 
                                  facecolor='orange',
                                  alpha=(win.max_cc_value **2) / 4)
                                            )
                    # Annotate time shift value into window,
                    # Alternate height if multiple windows so no overlap
                    axes[row][col].text(s=f"{tshift:.2f}s", x=tleft, 
                                        y=(ymax-ymin)*[0.75, 0.075][w%2]+ymin,
                                        fontsize=12, zorder=z+1)

                    axes[row][col].set_ylim([ymin, ymax])

            if row == 0:
                # determine how long the traces should be
                # hardcode the trace length based on user params
                if isinstance(trace_length, list):
                    axes[row][col].set_xlim(trace_length)
                else:
                    axes[row][col].set_xlim([np.maximum(time_offset, -10), 
                                             t[-1]
                                             ])

                # Set titles for the first row
                if col == middle_column:
                    title = (f"{st_obs[0].stats.network}."
                             f"{st_obs[0].stats.station} "
                             f"{format_event_name(ds)}")
                    axes[row][col].set_title(title)
        
                # Append component to bottom right of subplot  
                axes[row][col].text(
                            x=0.95, y=0.15, s=comp.upper(),
                            horizontalalignment="center",
                            verticalalignment="center",
                            transform=axes[row][col].transAxes)
      
        # y-label after all the processing has occurred
        axes[row][0].set_ylabel(ylab, rotation="horizontal", ha="right")

    # Label the time axis on the bottom row
    axes[-1][middle_column].set_xlabel("time [sec]")

    # Save the generated figure
    if save:
        plt.savefig(save)
    if show:
        plt.show()



