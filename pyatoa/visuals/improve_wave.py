#!/usr/bin/env python3
"""
Functions to create a figure showing progressive waveform improvement over
the course of a seismic inversion.

Show the changes in synthetic waveforms with progressive model updates. 
Each individual model gets its on row in the plot.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyasdf
from pyasdf import ASDFDataSet as asdf 
from pyatoa import Manager, logger
from pyatoa.utils.form import format_event_name
from pyatoa.visuals.wave_maker import format_axis
from pyflex import logger as pflogger

pflogger.setLevel("DEBUG")
logger.setLevel("INFO")


class ImproveWave:
    """
    A class to plot waveform improvement for a given ASDFDataSet

    .. code:: python

        ds = pyasdf.ASDFDataSet("dataset.h5")
        wi = WaveformImprovement(ds)
        wi.gather("NZ.BFZ", 10, 30)
        wi.plot()
        wi.plot("NZ.KNZ", 8, 30)
    """
    def __init__(self):
        """
        Initiate empty objects and keep dataset as an internal attribute

        :type ds: pyasdf.ASDFDataSet
        :param ds: dataset containing waveform data and windows
        """
        # self.ds = ds

        self.st_obs = None
        self.synthetics = None
        self.windows = None
        self.time_axis = None

    def get_models(self):
        """
        Figure out which step goes to which iteration to get model numbers
        """
        models = {"m00": "i01/s00"}
        iterations = self.ds.auxiliary_data.MisfitWindows.list()
        for iter_ in iterations:
            steps = self.ds.auxiliary_data.MisfitWindows[iter_]
            if iter_.replace("i", "m") in models:
                continue
            elif "s00" in steps.list():
                models[iter_.replace("i", "m")] = f"{iter_}/s00"
            else:
                models[iter_.replace("i", "m")] = \
                    f"{prev_iter}/{prev_steps.list()[-1]}"
            prev_iter = iter_
            prev_steps = steps

        # Get the last step
        if prev_steps.list()[-1] != "s00":
            final_model = f"m{int(prev_iter.split('i')[-1]) + 1:0>2}"
            models[final_model] = f"{prev_iter}/{prev_steps.list()[-1]}"

        return models

    def gather(self, sta, min_period, max_period, rotate_to_rtz=False,
               fix_windows=False, pyflex_preset=False):
        """
        Parse dataset for given station, gather observed and synthetic data, 
        preprocess data and return as stream objects.

        :type sta: str
        :param sta: station to gather data for
        :type min_period: float
        :param min_period: minimum filter period in seconds
        :type max_period: float
        :param max_period: maximum filter period in seconds
        :type rotate_to_rtz: bool
        :param rotate_to_rtz: rotate components from NEZ to RTZ
            Config. if False, instrument response will be removed from obs.
        :type fix_windows: bool
        :param fix_windows: dont recalculate windows when gathering
        :type pyflex_preset: str
        :param pyflex_preset: overwrite the pyflex preset provided in the
            Config object
        """
        if min_period is None or max_period is None:
            raise TypeError("must specify 'min_period' and 'max_period'")

        assert(sta in self.ds.waveforms.list()), f"{sta} not in ASDFDataSet"

        models = self.get_models()

        # Preprocess all traces using Pyatoa and store in dict
        st_obs, synthetics, windows = None, {}, {}
        for model, path in models.items():
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
            iter_, step_ = path.split("/")
            mgmt.window(fix_windows=fix_windows, iteration=iter_,
                        step_count=step_)

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

    def gather_simple(self, event, sta, min_period, max_period, path_dict=None,
                      component=None):
        """
        Manually set the model values based on inspection of the Inspector
        Don't return windows or anything, keep it simple
        """

        models = {"m00": ("i01/s00", "a"),
                  "m03": ("i03/s03", "a"),
                  "m09": ("i09/s02", "a"),
                  "m12": ("i12/s04", "a"),
                  "m17": ("i17/s01", "a"),
                  "m24": ("i07/s01", "b"),
                  "m28": ("i11/s03", "c"),
                  }
        st_obs, synthetics = None, {}
        windows = None
        for model, tup in models.items():
            path, tag = tup
            if path_dict:
                ds_fid = os.path.join(path_dict[tag], f"{event_id}.h5")
            else:
                ds_fid = f"{event_id}{tag}.h5"
            with asdf(ds_fid, mode="r") as ds:
                mgmt = Manager(ds=ds)
                mgmt.load(sta, path)
                mgmt.config.save_to_ds = False
                mgmt.config.min_period = min_period
                mgmt.config.max_period = max_period
                mgmt.standardize().preprocess()

                if component:
                    synthetics[model] = mgmt.st_syn.select(
                                                     component=component).copy()
                else:
                    synthetics[model] = mgmt.st_syn.copy()
                if st_obs is None:
                    if component:
                        st_obs = mgmt.st_obs.select(component=component).copy()
                    else:
                        st_obs = mgmt.st_obs.copy()

        self.st_obs = st_obs
        self.synthetics = synthetics

    def setup_plot(self, nrows, ncols, **kwargs):
        """
        Dynamically set up plots according to number_of given
        Returns a list of lists of axes objects
        e.g. axes[i][j] gives the ith column and the jth row

        :type nrows: int
        :param nrows: number of rows in the gridspec
        :type ncols: int
        :param ncols: number of columns in the gridspec
        :rtype axes: matplotlib axes
        :return axes: axis objects
        """
        dpi = kwargs.get("dpi", 150)
        figsize = kwargs.get("figsize", (500/dpi, 800/dpi))
        fontsize = kwargs.get("fontsize", 10)
        axis_linewidth = kwargs.get("axis_linewidth", 2)

        f = plt.figure(figsize=figsize, dpi=dpi)
        gs = mpl.gridspec.GridSpec(nrows, ncols, hspace=0, wspace=0.025, 
                                   width_ratios=[1] * ncols,
                                   height_ratios=[3] * nrows
                                   )

        axes = [[] for _ in range(nrows)]
        for row in range(0, gs.get_geometry()[0]):
            for col in range(0, gs.get_geometry()[1]):
                # Ensure axis sharing
                if col == 0:
                    sharey = None
                else:
                    sharey = axes[row][0]
                if row == 0 and col == 0:
                    sharex = None
                else:
                    sharex = axes[0][0]

                ax = plt.subplot(gs[row, col], sharey=sharey, sharex=sharex)
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

                axes[row].append(ax)


        # remove x-tick labels except for last axis
        for row in axes[:-1]:
            for col in row:
                plt.setp(col.get_xticklabels(), visible=False)

        return f, axes

    def plot(self, sta=None, event_id=None, min_period=None, max_period=None, 
             plot_windows=False, trace_length=None, show=True, save=False, 
             **kwargs):
        """
        Plot waveforms iterative based on model updates

        :type sta: str
        :param sta: station to gather data for, if None, skips gathering
            assuming data has already been gathered
        :type min_period: float
        :param min_period: minimum filter period for waveforms
        :type max_period: float
        :param max_period: maximum filter period for waveforms
        :type plot_windows: bool
        :param plot_windows: plot misfit windows above waveforms
        :type trace_length: list of floats
        :param trace_length: [trace_start, trace_end] will be used to set the x
            limit on the waveform data. If none, no xlim will be set
        :type show: bool
        :param show: Show the plot or do not
        :type save: str
        :param save: if given, save the figure to this path
        """
        linewidth = kwargs.get("linewidth", 1.)
        fontsize = kwargs.get("fontsize", 10)
        anno_fontsize = kwargs.get("anno_fontsize", 8)
        window_color = kwargs.get("window_color", "orange")
        percent_over = kwargs.get("percent_over", 0.125)
        anno_choice = kwargs.get("anno_choice", "all")

        # Allows for skipping the gather call and including it directly in plot
        if sta is not None:
            self.gather(sta, min_period, max_period)

        assert self.st_obs, "must collect data for a station before plotting"

        # Instantiate the plotting object
        f, axes = self.setup_plot(nrows=len(self.synthetics.keys()), 
                                  ncols=len(self.st_obs), **kwargs)

        # if not trace_length:
        #     trace_length = [self.time_axis[0], self.time_axis[-1]]

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
                a1, = axes[row][col].plot(obs.times(), obs.data, 'k', 
                                          zorder=10, label="Obs", 
                                          linewidth=linewidth)
                a2, = axes[row][col].plot(syn.times(), syn.data, 
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
                            # tright_max = tright

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
                    # else:
                    #     axes[row][col].set_xlim([
                    #         np.maximum(self.time_axis[0], -10), t[-1]
                    #         ])

                    # Set titles for the first row, middle column
                    if col == len(self.st_obs) // 2:
                        title = (f"{self.st_obs[0].stats.network}."
                                 f"{self.st_obs[0].stats.station} "
                                 f"{event_id} Z")
                        axes[row][col].set_title(title, fontsize=fontsize)
            
                    # Append component to bottom right of subplot  
                    if False:
                        axes[row][col].text(
                                    x=0.95, y=0.15, s=comp.upper(),
                                    horizontalalignment="center",
                                    verticalalignment="center",
                                    transform=axes[row][col].transAxes)
          
            # y-label after all the processing has occurred
            axes[row][0].set_ylabel(ylab, rotation="horizontal", ha="right",
                                    fontsize=fontsize)

        # Label the time axis on the bottom row, middle column
        axes[-1][len(self.st_obs) // 2].set_xlabel("time [s]", 
                                                   fontsize=fontsize)

        # Save the generated figure
        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, axes

    def gather_simple(self, models, event_id, sta, component,  min_period,
                      max_period):
        """
        Gather waveforms from manually input model values, usually determined
        by using the Inspector class

        :type models: dict of tuples
        :param models: model values as keys, (iter/step, tag) as tuple value.
            Tags allow multiple datasets to be used, e.g. if an inversion
            spans over multiple legs and more than one dataset is used to
            store waveform data
        :type event_id: str
        :param event_id: name of the event, used to access the ASDFDataSet
        ;type sta: str
        :param sta: station id to gather data for
        :type min_period: float
        :param min_period: period to filter data at
        :type max_period: float
        :param max_period: period to filter data at
        """
        st_obs, synthetics = None, {}

if __name__ == "__main__":
    pairs = [
         ("2013p617227", "NZ.TOZ", "Z"),
         # ("2014p952799", "NZ.NTVZ", "N"),
         # ("2016p105478", "NZ.PUZ", "Z"),
         # ("2016p881118", "NZ.MWZ", "E"),
         # ("2018p465580", "NZ.KHEZ", "E"),
         # ("2019p738432", "NZ.KHZ", "Z"),
         # ("2019p754447", "NZ.HIZ", "Z"),
         # ("2019p927023", "NZ.VRZ", "Z"),
         ]    
    
    path_dict = {"a": "../waveform_comparisons/aspen/",
                 "b": "/home/chowbr/current/birch/scratch/preprocess/datasets/i07_corrupted/",
                 "c": "../waveform_comparisons/birch/"}

    for event_id, sta, comp in pairs:
        wi = ImproveWave()
        wi.gather_simple(event_id, sta, 6, 30, path_dict=path_dict, component=comp)
        wi.plot(show=False, save=f"{event_id}_{sta}.png", event_id=event_id,
                trace_length=[70,290])
    a=1/0

    # MAIN
    event_id = "2019p927023"
    with asdf(f"{event_id}a.h5") as ds:
        stations = ds.waveforms.list()
        self.st_obs = st_obs
        self.synthetics = synthetics


if __name__ == "__main__":
    event_id = "2013p507880"
    component = "Z"
    models = {"m00": ("i01/s00", ""),
              "m01": ("i01/s04", ""),
              "m02": ("i02/s01", ""),
              "m02": ("i03/s00", ""),
              "m03": ("i03/s01", ""),
              "m03": ("i04/s00", ""),
              "m04": ("i04/s04", ""),
              "m04": ("i05/s00", ""),
              "m05": ("i05/s04", ""),
              "m05": ("i06/s00", ""),
              "m06": ("i06/s01", ""),
              "m07": ("i07/s01", ""),
              "m08": ("i08/s03", ""),
              "m09": ("i09/s01", ""),
              "m10": ("i10/s04", ""),
              "m10": ("i11/s00", ""),
              "m11": ("i11/s03", ""),
              "m12": ("i12/s01", ""),
              "m13": ("i13/s01", ""),}

    # with asdf(f"{event_id}.h5") as ds:
    #     stations = ds.waveforms.list()
    stations = ["NZ.TOZ"]

        
    for sta in stations:
        wi = ImproveWave()
        wi.gather_simple(models, event_id, sta, component, 4, 30)
        wi.plot()
