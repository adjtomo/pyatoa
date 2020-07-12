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
  
from pyatoa import Config, Manager, logger
from pyatoa.utils.form import format_event_name
from pyatoa.visuals.manager_plotter import format_axis
from pyatoa.utils.asdf.fetch import windows_from_dataset
from pyflex import logger as pflogger

pflogger.setLevel("DEBUG")
logger.setLevel("INFO")

class WaveformImprovement:
    """
    A class to plot waveform improvement for a given ASDFDataSet

    Use case:
        ds = pyasdf.ASDFDataSet("dataset.h5")
        wi = WaveformImprovement(ds)
        wi.gather("NZ.BFZ", 10, 30)
        wi.plot()
        wi.plot("NZ.KNZ", 8, 30)
    """
    def __init__(self, ds):
        """
        Initiate empty objects and keep dataset as an internal attribute

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset containing waveform data and windows
        """
        self.ds = ds

        self.st_obs = None
        self.synthetics = None
        self.windows = None
        self.time_axis = None

    
    def sort_steps(self):
        """
        Figure out which step goes to which model
        """
        model_dict = {"m00": "m00/s00"}
        models = self.ds.auxiliary_data.MisfitWindows.list()
        for model in models:
            steps = self.ds.auxiliary_data.MisfitWindows[model]
            if model in model_dict:
                continue
            elif "s00" in steps.list():
                model_dict[model] = f"{model}/s00"
            else:
                model_dict[model] = f"{prev_model}/{prev_steps.list()[-1]}"
            prev_model = model
            prev_steps = steps

        # Get the last step
        if prev_steps.list()[-1] != "s00":
            final_model = f"m{int(prev_model.split('m')[-1]) + 1:0>2}"
            model_dict[final_model] = f"{prev_model}/{prev_steps.list()[-1]}"

        return model_dict


    def gather(self, sta, min_period, max_period, rotate_to_rtz=False,
               synthetics_only=False, fix_windows=False, pyflex_preset=False):
        """
        Parse dataset for given station, gather observed and synthetic data, 
        preprocess data and return as stream objects.

        :type synthetics_only: bool
        :param synthetics_only: synthetics only parameter to be passed to the
            Config. if False, instrument response will be removed from obs.
        """
        if min_period is None or max_period is None:
            raise TypeError("must specify 'min_period' and 'max_period'")

        assert(sta in self.ds.waveforms.list()), f"{sta} not in ASDFDataSet"

        network, station = sta.split(".")
        models = self.ds.auxiliary_data.MisfitWindows.list()

        model_dict = self.sort_steps()

        # Preprocess all traces using Pyatoa and store in dict
        st_obs, synthetics, windows = None, {}, {}
        for model, path in model_dict.items():
            # Gather synthetic data
            mgmt = Manager(ds=self.ds)   
            print(path)
            mgmt.load(sta, path)

            # Overwrite some config parameters
            mgmt.config.save_to_ds = False
            mgmt.config.min_period = min_period
            mgmt.config.max_period = max_period
            if rotate_to_rtz:
                mgmt.config.rotate_to_rtz = rotate_to_rtz
                mgmt.config.component_list = ["Z", "R", "T"]
            if pyflex_preset:
                mgmt.config.pyflex_preset = pyflex_preset
                mgmt.config._check()

            mgmt.standardize()
            mgmt.preprocess() 
            mod_, step_ = path.split("/")
            mgmt.window(fixed=fix_windows, model=mod_, step=step_)

            windows[model] = mgmt.windows
            synthetics[model] = mgmt.st_syn.copy() 

            # Observed waveform will be the same
            if st_obs is None:
                st_obs = mgmt.st_obs.copy()

        # Internally used by plotting function
        self.st_obs = st_obs
        self.synthetics = synthetics
        self.windows = windows
        self.time_axis = self.st_obs[0].times(
            reftime=st_obs[0].stats.starttime - mgmt.stats.time_offset_sec
            ) 

    def setup_plot(self, nrows, ncols, label_units=False, **kwargs):
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
            if there are multiple models plotted together, so usually best to 
            leave it off.
        :rtype axes: matplotlib axes
        :return axes: axis objects
        """
        figsize = kwargs.get("figsize", (6, 8))
        dpi = kwargs.get("dpi", 150)
        fontsize = kwargs.get("fontsize", 10)
        axis_linewidth = kwargs.get("axis_linewidth", 2)

        f = plt.figure(figsize=figsize, dpi=dpi)
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
                                     # sharey=axes[0][col]
                                     )
                ax.set_axisbelow(True)
                ax.minorticks_on()
                ax.tick_params(which='major', direction='in', top=True,
                               right=False, left=False, labelsize=fontsize, 
                               length=3, width=2*axis_linewidth/3)
                ax.tick_params(which='minor', direction='in', length=1.5, 
                               top=True, bottom=True, right=False, left=False, 
                               width=2*axis_linewidth/3)

                for axis in ["top", "bottom", "left", "right"]:
                    ax.spines[axis].set_linewidth(axis_linewidth)
                
                # Turn off the y axes because we wont show units
                ax.get_yaxis().set_ticks([])

                # Set the grids on
                components.append(ax)
            axes.append(components)

        # remove x-tick labels except for last axis
        for row in axes[:-1]:
            for col in row:
                plt.setp(col.get_xticklabels(), visible=False)
        
        # Deprecated, prefer to just set ticks to empty, will probably never
        # need to show y-values/labels if showing lots of waveforms
        # 
        # # remove y-tick labels except for first axis
        # if label_units:
        #     for i, row in enumerate(axes):
        #         for j, col in enumerate(row):
        #             if j != 0: 
        #                 plt.setp(col.get_yticklabels(), visible=False)
        #             else:    
        #                 col.ticklabel_format(style='sci', axis='y', 
        #                                      scilimits=(0,0))
        # else:
        #     for i, row in enumerate(axes):
        #         for col in row:
        #             plt.setp(col.get_yticklabels(), visible=False)

        return f, axes


    def plot(self, sta=None, min_period=None, max_period=None, 
             plot_windows=False, trace_length=None, show=True, save=False, 
             **kwargs):

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
        linewidth = kwargs.get("linewidth", 2)
        fontsize = kwargs.get("fontsize", 10)
        anno_fontsize = kwargs.get("anno_fontsize", 8)
        window_color = kwargs.get("window_color", "orange")
        percent_over = kwargs.get("percent_over", 0.125)
        anno_choice = kwargs.get("anno_choice", "all")

        # Allows for skipping the gather call and including it directly in plot
        if sta is not None:
            self.gather(sta, min_period, max_period)

        assert(self.st_obs), "must collect data for a station before plotting"

        # Instantiate the plotting object
        f, axes = self.setup_plot(nrows=len(self.synthetics.keys()), 
                                  ncols=len(self.st_obs), **kwargs)

        if not trace_length:
            trace_length = [self.time_axis[0], self.time_axis[-1]]

        # Plot each model on a different row
        synthetic_keys = list(self.synthetics.keys())
        synthetic_keys.sort()
        for row, syn_key in enumerate(synthetic_keys):
            ylab = syn_key.split('_')[-1]  # e.g. 'm00'

            # Plot each component in a different column
            component_list = [_.stats.channel[-1] for _ in self.st_obs]
            for col, comp in enumerate(component_list):
                obs = self.st_obs.select(component=comp)[0]
                syn = self.synthetics[syn_key].select(component=comp)[0]

                # Plot waveforms
                a1, = axes[row][col].plot(self.time_axis, obs.data, 'k', 
                                          zorder=10, label="Obs", 
                                          linewidth=linewidth)
                a2, = axes[row][col].plot(self.time_axis, syn.data, 
                                          ["r", "b", "g"][col], zorder=10,
                                          label="Syn", linewidth=linewidth)

                # Format the axes for a standardized look
                format_axis(axes[row][col], percent_over=percent_over)

                # Plot windows if available for this given component
                tshift_max = 0  # temporary
                if plot_windows:
                    windows = self.windows[syn_key].get(comp, [])
                    for w, win in enumerate(windows):
                        ymin, ymax = axes[row][col].get_ylim()

                        tleft = win.left * win.dt + self.time_axis[0]
                        tright = win.right * win.dt + self.time_axis[0]
                        tshift = win.cc_shift * win.dt

                        axes[row][col].add_patch(mpl.patches.Rectangle(
                                      xy=(tleft, ymin), width=tright-tleft, 
                                      ec='k', fc=window_color,
                                      height=(ymax + np.abs(ymin)), 
                                      alpha=(win.max_cc_value **2) / 4)
                        )
                        # Outline the rectangle with solid lines
                        for t_ in [tleft, tright]:
                            axes[row][col].axvline(x=t_, ymin=0, ymax=1, 
                                                   color="k", alpha=1., 
                                                   zorder=11)

                        # Annotate time shift value into window,
                        # Alternate height if multiple windows so no overlap
                        if anno_choice == "all":
                            axes[row][col].text(
                                s=f"{tshift:.2f}s", x=tleft, 
                                y=(ymax-ymin)*[0.7, 0.06][w%2]+ymin,
                                fontsize=anno_fontsize, zorder=11
                                )

                        # If only annotate the largest timeshift
                        if abs(tshift) > abs(tshift_max):
                            tshift_max = tshift
                            tleft_max = tleft
                            tright_max = tright

                    # If annotate largesttime shift value into window
                    if tshift_max and anno_choice == "max":
                        axes[row][col].text(
                            s=f"{tshift_max:.2f}s", x=tleft_max, 
                            y=(ymax-ymin)*0.06 + ymin,
                            fontsize=anno_fontsize, zorder=11
                            )

                if row == 0:
                    # determine how long the traces should be
                    # hardcode the trace length based on user params
                    if isinstance(trace_length, list):
                        axes[row][col].set_xlim(trace_length)
                    else:
                        axes[row][col].set_xlim([
                            np.maximum(self.time_axis[0], -10), t[-1]
                            ])

                    # Set titles for the first row, middle column
                    if col == len(self.st_obs) // 2:
                        title = (f"{self.st_obs[0].stats.network}."
                                 f"{self.st_obs[0].stats.station} "
                                 f"{format_event_name(self.ds)}")
                        axes[row][col].set_title(title, fontsize=fontsize)
            
                    # Append component to bottom right of subplot  
                    axes[row][col].text(
                                x=0.95, y=0.15, s=comp.upper(),
                                horizontalalignment="center",
                                verticalalignment="center",
                                transform=axes[row][col].transAxes)
          
            # y-label after all the processing has occurred
            axes[row][0].set_ylabel(ylab, rotation="horizontal", ha="right",
                                    fontsize=fontsize)

        # Label the time axis on the bottom row, middle column
        axes[-1][len(self.st_obs) // 2].set_xlabel("time [sec]", 
                                                   fontsize=fontsize)

        # Save the generated figure
        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, axes



