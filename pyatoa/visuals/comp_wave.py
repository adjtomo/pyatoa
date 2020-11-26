#!/usr/bin/env python3
"""
A stripped down (arguably better) version of the Waveform Improvement class,
used to simply compare two synthetic waveforms from a given PyASDF DataSet.
"""
import os
from glob import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa import Manager
from pyasdf import ASDFDataSet
from pyatoa.visuals.wave_maker import format_axis
from pyatoa.visuals.map_maker import MapMaker
from obspy.signal.filter import envelope



class CompWave:
    """
    A class to plot waveform improvement between two models for a given dataset
    """
    def __init__(self, dsfid, station, min_period, max_period,
                 dsfid_final=None):
        """
        Initiate empty objects and keep dataset as an internal attribute

        :type dsfid: str
        :param dsfid: name of the ASDF Dataset
        :type station: str
        :param station: Network and station code, e.g. NZ.BFZ
        :type min_period: float
        :param min_period: minimum filter corner
        :type max_period: float
        :param max_period: maximum filter corner
        """
        self.dsfid = dsfid
        self.dsfid_final = dsfid_final
        self.event_id = os.path.splitext(os.path.basename(dsfid))[0]
        self.station = station
        self.min_period = min_period
        self.max_period = max_period

        # Initialize empty attributes to be filled
        self._ds = None
        self._inv = None
        self._event = None
        self._st_obs = None
        self._st_init = None
        self._st_final = None
        self._m_init = None
        self._m_final = None

    def _gather_model_from_dataset(self, dsfid, model=None, init_or_final=None):
        """
        Gather data from an ASDFDataSet based on the given model (iter/step)

        :type dsfid: str
        :param dsfid: file identifier for the dataset
        :type model: str
        :param model: iteration/step, e.g. 'i01/s00'
        :type init_or_final: str
        :param init_or_final: for choosing default values if model is None
            * 'init': choose the first iteration/step for the initial model
            * 'final': choose final iteration/step for final model
        :return:
        """
        assert(init_or_final in ["init", "final"])

        with ASDFDataSet(dsfid, mode="r") as ds:
            if model is None:
                configs = ds.auxiliary_data.Configs
                if init_or_final is "init":
                    idx = 0
                elif init_or_final is "final":
                    idx = -1
                iter_ = configs.list()[idx]
                step_ = configs[iter_].list()[idx]
                model = f"{iter_}/{step_}"

            # Use the Manager class to load in waveform data
            mgmt = Manager(ds=ds)
            mgmt.load(code=self.station, path=model)
            # Overwrite the filter corners stored in the dataset
            mgmt.config.min_period = self.min_period
            mgmt.config.max_period = self.max_period
            mgmt.standardize().preprocess()

            # Store data in class attributes, obs waveforms will be the same
            if self._st_obs is None:
                self._st_obs = mgmt.st_obs
            if self._inv is None:
                self._inv = mgmt.inv
            if self._event is None:
                self._event = mgmt.event
            setattr(self, f"_st_{init_or_final}", mgmt.st_syn)
            setattr(self, f"_m_{init_or_final}", model)

    def gather(self, m_init=None, m_final=None):
        """
        Gather data from the correct dataset. If no m_init or m_final given,
        will gather the first and last models

        :type m_init: str
        :param m_init: initial iteration/step, e.g. 'i01/s00'
        :type m_final: str
        :param m_final: final iteration/step
        """
        self._gather_model_from_dataset(dsfid=self.dsfid, model=m_init,
                                        init_or_final="init")

        # Default to dsfid if separate final dataset not provided
        self._gather_model_from_dataset(dsfid=self.dsfid_final or self.dsfid,
                                        model=m_final, init_or_final="final")

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
        figure = kwargs.get("figure", None)
        subplot_spec = kwargs.get("subplot_spec", None)
        dpi = kwargs.get("dpi", 100)
        figsize = kwargs.get("figsize", (1400 / dpi, 600 / dpi))
        fontsize = kwargs.get("fontsize", 8)
        axis_linewidth = kwargs.get("axis_linewidth", 1.75)

        # Initiate the figure, allow for external figure objects
        if figure is None:
            f = plt.figure(figsize=figsize, dpi=dpi)
        else:
            f = figure

        # Initiate gridspec, allow for nested grids
        subplot_kwargs = {"hspace": 0, "wspace": 0.025, 
                          "width_ratios": [5] * ncols, 
                          "height_ratios": [1] * nrows}
        if subplot_spec is None:
            gs = mpl.gridspec.GridSpec(nrows, ncols, **subplot_kwargs)
        else:
            gs = mpl.gridspec.GridSpecFromSubplotSpec(nrows, ncols,
                                                      subplot_spec=subplot_spec,
                                                      **subplot_kwargs)

        axes = [[] for _ in range(nrows)]
        for row in range(0, gs.get_geometry()[0]):
            for col in range(0, gs.get_geometry()[1]):
                # First column can't share y-values
                if col == 0:
                    sharey = None
                else:
                    sharey = axes[row][0]
                # First entry can't share x-values
                if row == 0 and col == 0:
                    sharex = None
                else:
                    sharex = axes[0][0]

                ax = plt.subplot(gs[row, col], sharey=sharey, sharex=sharex)

                ax.set_axisbelow(True)
                ax.minorticks_on()
                ax.tick_params(which='major', direction='in', top=True,
                               right=False, left=True, labelsize=fontsize,
                               length=3, width=2*axis_linewidth/3)
                ax.tick_params(which='minor', direction='in', length=1.5, 
                               top=True, bottom=True, right=False, left=True,
                               width=2*axis_linewidth/3)
                if col == 0:
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                    # Make sure the scientific notation has the same fontsize
                    ax.yaxis.get_offset_text().set_fontsize(fontsize)

                for axis in ["top", "bottom", "left", "right"]:
                    ax.spines[axis].set_linewidth(axis_linewidth)

                # Set the grids on
                axes[row].append(ax)

        # remove x-tick labels except for last axis
        for row in axes[:-1]:
            for col in row:
                plt.setp(col.get_xticklabels(), visible=False)

        # Turn off the y-tick labels and sci notation for columns except first
        for row in axes:
            for col in row[1:]:
                plt.setp(col.get_yticklabels(), visible=False)
                col.yaxis.get_offset_text().set_visible(False)

        return f, axes

    def _xlim_from_envelope(self, data, dt):
        """
        Get rough bounds for the xlimits by looking at waveform envelopes
        :return:
        """
        env = envelope(data)
        idx = np.where(env >= env.std())[0] * dt
        return [np.floor(idx[0]) - 20 , np.ceil(idx[-1]) + 20]

    def plot(self, show=True, save=False, **kwargs):
        """
        Plot waveforms iterative based on model updates

        :type show: bool
        :param show: Show the plot or do not
        :type save: str
        :param save: if given, save the figure to this path
        """
        linewidth = kwargs.get("linewidth", 1.3)
        fontsize = kwargs.get("fontsize", 8)
        xlim = kwargs.get("xlim", None)
        percent_over = kwargs.get("percent_over", 0.125)
        color_init = kwargs.get("color_init", "orangered")
        color_final = kwargs.get("color_final", "mediumorchid")

        # One row per component, one column for init and final each
        f, axes = self.setup_plot(nrows=3, ncols=2, **kwargs)

        # Plot each component in a different column
        component_list = [_.stats.channel[-1] for _ in self._st_obs]
        for row, comp in enumerate(component_list):
            obs = self._st_obs.select(component=comp)[0]
            syn_init = self._st_init.select(component=comp)[0]
            syn_final = self._st_final.select(component=comp)[0]

            # Plot init synthetics to the left, final synthetics to the right
            axes[row][0].plot(syn_init.times(), syn_init.data, color_init,
                              zorder=10, label="INIT", linewidth=linewidth
                              )
            axes[row][1].plot(syn_final.times(), syn_final.data, color_final,
                              zorder=10, label="FINAL", linewidth=linewidth
                              )
            # Plot obs in both columns
            for col in [0, 1]:
                axes[row][col].plot(obs.times(), obs.data, "k", zorder=11,
                                    label="OBS", linewidth=linewidth)

            # Format the axes for a standardized look
            for col in range(len(axes[row])):
                format_axis(axes[row][col], percent_over=percent_over)

            # Component in the y-label
            axes[row][0].set_ylabel(comp.upper(), rotation="horizontal",
                                    ha="left", va="center")

            # Set ylim based on ymax from either waveform
            max_yvals = []
            for data in [obs.data, syn_init.data, syn_final.data]:
                max_yvals.append(max(abs(data)))
            ylim = [-1 * max(max_yvals), max(max_yvals)]
            axes[row][0].set_ylim(ylim)

        # Set xlim for master axis
        if xlim is None:
            xlim = [self._st_obs[0].times()[0], self._st_obs[0].times()[-1]]
        elif xlim == "auto":
            xlim = self._xlim_from_envelope(
                data=self._st_final.select(component="Z")[0].data,
                dt=self._st_final[0].stats.delta
            )
        axes[0][0].set_xlim(xlim)

        # Main title ab"ove all the subplots
        plt.suptitle(f"Waveform Improvement | "
                     f"{self.event_id} {self.station} "
                     f"[{self.min_period}, {self.max_period}]s",
                     fontsize=fontsize)

        # Common X and Y labels, manually decided values
        f.text(0.375, 0.05, "Time [s]", ha="center", fontsize=fontsize)
        f.text(0.09, 0.375, "Displacement [m]", ha="center", rotation=90,
               fontsize=fontsize)

        # Title each of the models
        # axes[0][0].set_title(f"INITIAL: {self._m_init}", fontsize=fontsize)
        # axes[0][1].set_title(f"FINAL: {self._m_final}", fontsize=fontsize)
        axes[0][0].set_title(f"Initial", fontsize=fontsize)
        axes[0][1].set_title(f"Final", fontsize=fontsize)

        # Save the generated figure
        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, axes

    def plot_with_map(self, corners=None, dpi=100, figsize=None, show=True,
                      save=False,**kwargs):
        """
        Similar to Manager plotter, plot the waveform comparisons next to a
        source receiver map. Wraps the internal plotting functionality with
        a gridspec
        """
        # Default figure size
        if figsize is None:
            figsize = (2400 / dpi, 600 / dpi)

        # Create an overlying GridSpec that will contain both plots
        gs = mpl.gridspec.GridSpec(1, 2, wspace=0.0, hspace=0.,
                                   width_ratios=[2, 1], height_ratios=[1]
                                   )
        fig = plt.figure(figsize=figsize, dpi=dpi)

        # Plot the waveform on the left
        self.plot(figure=fig, subplot_spec=gs[0],  show=False, save=False,
                  fontsize=12)

        # Plot the map on the right
        mm = MapMaker(inv=self._inv, cat=self._event, **kwargs)
        ax = fig.add_subplot(gs[1])
        mm.plot(corners=corners, figure=fig, ax=ax, show=False, save=False,
                **kwargs)

        if save:
            plt.savefig(save)
        if show:
            plt.show()
        else:
            plt.close()


def main(event_id=None, station=None):
    """
    Main call script to choose event and station based on what's available
    """
    # =========================================================================
    # PARAMETER CHOICE
    choice = "all"  # pick or all
    min_period = 6
    max_period = 30
    m_init = None
    m_final = "i10/s02"
    event_id = event_id
    dsfid = f"./aspen/{event_id}.h5"
    dsfid_final = f"./birch/{event_id}.h5"
    station = station
    show = False
    plot_with_map = True
    xlim = None
    # =========================================================================

    # Get station information prior to plotting
    with ASDFDataSet(dsfid, mode="r") as ds:
        stations = ds.waveforms.list()

    # Ask user to choose station to plot
    if choice == "pick":
        for s, sta in enumerate(stations):
            print(f"{s}: {sta}")
        while True:
            idx = input(f"Which station index or name?: ")
            try:
                stations = [stations[int(idx)]]
            except ValueError:
                stations = [idx]

    for sta in stations:
        if station is not None and sta != station:
            continue

        fid_out = f"./figures/{event_id}_{sta}.png"
        if os.path.exists(fid_out):
            print(f"{fid_out} exists")
            continue
        print(sta)
        try:
            cw = CompWave(dsfid=dsfid, dsfid_final=dsfid_final,
                          station=sta, min_period=min_period,
                          max_period=max_period)
            cw.gather(m_init, m_final)

            if plot_with_map:
                cw.plot_with_map(show=show, save=fid_out, xlim=xlim)
            else:
                cw.plot(show=show, save=fid_out, xlim=xlim)
            plt.close()
        except Exception as e:
            print(e)
            pass


if __name__ == "__main__":
    for fid in sorted(glob("./aspen/*.h5")):
        event_id = os.path.splitext(os.path.basename(fid))[0]
        print(event_id)
        main(event_id, None)
