#!/usr/bin/env python3
"""
A stripped down (arguably better) version of the Waveform Improvement class,
used to simply compare two synthetic waveforms from a given PyASDF DataSet.
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa import Manager
from pyasdf import ASDFDataSet
from pyatoa.visuals.wave_maker import format_axis


class CompWave:
    """
    A class to plot waveform improvement between two models for a given dataset
    """
    def __init__(self, dsfid, station, min_period, max_period):
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
        self.event_id = os.path.splitext(os.path.basename(dsfid))[0]
        self.station = station
        self.min_period = min_period
        self.max_period = max_period

        # Initialize empty attributes to be filled
        self._ds = None
        self._st_obs = None
        self._st_init = None
        self._st_final = None
        self._m_init = None
        self._m_final = None

    def gather(self, m_init=None, m_final=None):
        """
        Gather data from the correct dataset. If no m_init or m_final given,
        will gather the first and last models

        :type m_init: str
        :param m_init: initial iteration/step, e.g. 'i01/s00'
        :type m_final: str
        :param m_final: final iteration/step
        """
        with ASDFDataSet(self.dsfid) as ds:
            # Figure out initial and final model iterations and step counts
            configs = ds.auxiliary_data.Configs
            if m_init is None:
                iter_init = configs.list()[0]
                step_init = configs[iter_init].list()[0]
                m_init = f"{iter_init}/{step_init}"

            if m_final is None:
                iter_final = configs.list()[-1]
                step_final = configs[iter_final].list()[-1]
                m_final = f"{iter_final}/{step_final}"

            # Use the Manager class to load in waveform data
            mgmt = Manager(ds=ds)
            mgmt.config.min_period = self.min_period
            mgmt.config.max_period = self.max_period

            # Store data in class attributes
            for path, attr in zip([m_init, m_final], ["_st_init", "_st_final"]):
                mgmt.load(code=self.station, path=path)
                mgmt.standardize().preprocess()
                setattr(self, attr, mgmt.st_syn)
            # Obs waveforms will be the same both times
            self._st_obs = mgmt.st_obs

        self._m_init = m_init
        self._m_final = m_final

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
        figsize = kwargs.get("figsize", (8, 5))
        dpi = kwargs.get("dpi", 100)
        fontsize = kwargs.get("fontsize", 8)
        axis_linewidth = kwargs.get("axis_linewidth", 1.5)

        f = plt.figure(figsize=figsize, dpi=dpi)
        gs = mpl.gridspec.GridSpec(nrows, ncols, hspace=0, wspace=0.075,
                                   width_ratios=[3] * ncols,
                                   height_ratios=[1] * nrows
                                   )

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

    def plot(self, m_init=None, m_final=None, show=True, save=False, **kwargs):
        """
        Plot waveforms iterative based on model updates

        :type m_init: str
        :param m_init: initial iteration/step, e.g. 'i01/s00'
        :type m_final: str
        :param m_final: final iteration/step
        :type show: bool
        :param show: Show the plot or do not
        :type save: str
        :param save: if given, save the figure to this path
        """
        linewidth = kwargs.get("linewidth", 1.75)
        fontsize = kwargs.get("fontsize", 8)
        xlim = kwargs.get("xlim", None)
        percent_over = kwargs.get("percent_over", 0.125)
        color_init = kwargs.get("color_init", "red")
        color_final = kwargs.get("color_final", "blue")

        self.gather(m_init, m_final)

        # One row per component, one column for init and final each
        f, axes = self.setup_plot(nrows=3, ncols=3, **kwargs)

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
                axes[row][col].plot(obs.times(), obs.data, "k", zorder=9,
                                    label="OBS", linewidth=linewidth)

            # Plot both syns together in last column
            axes[row][2].plot(syn_init.times(), syn_init.data, color_init, 
                              zorder=10, label="INIT", linewidth=linewidth
                              )
            axes[row][2].plot(syn_final.times(), syn_final.data, color_final,
                              zorder=10, label="FINAL", linewidth=linewidth
                              )

            # Format the axes for a standardized look
            for col in range(len(axes[row])):
                format_axis(axes[row][col], percent_over=percent_over)

            # Component in the y-label
            axes[row][0].set_ylabel(comp.upper(), rotation="horizontal",
                                    ha="left", va="center")

        # Set xlim for master axis
        if xlim is None:
            xlim = [self._st_obs[0].times()[0], self._st_obs[0].times()[-1]]
        axes[0][0].set_xlim(xlim)

        # Main title above all the subplots
        plt.suptitle(f"{self.event_id} {self.station} "
                     f"[{self.min_period}, {self.max_period}]s",
                     fontsize=fontsize)

        # Common x-label location
        f.text(0.5, 0.04, "Time [s]", ha="center", fontsize=fontsize)

        # Title each of the models
        axes[0][0].set_title(f"INITIAL: {self._m_init}", fontsize=fontsize)
        axes[0][1].set_title(f"FINAL: {self._m_final}", fontsize=fontsize)
        axes[0][2].set_title(f"INITIAL / FINAL", fontsize=fontsize)

        # Save the generated figure
        if save:
            plt.savefig(save)
        if show:
            plt.show()

        return f, axes


def plot_from_working_dir():
    """
    Main call script to choose event and station based on what's available
    :return:
    """
    from glob import glob
    choice = "all"  # pick or all

    # Hardcoded
    min_period = 6
    max_period = 30
    m_init = None
    m_final = None

    available = glob("*.h5")
    if len(available) == 1:
        idx = 0
    else:
        for a, avail in enumerate(available):
            print(f"{a:0>3}: {avail}")
        idx = int(input("Which file id index?: "))

    dsfid = available[idx]

    with ASDFDataSet(dsfid) as ds:
        stations = ds.waveforms.list()
    for s, sta in enumerate(stations):
        print(f"{s:0>3}: {sta}")

    # Ask user to choose station to plot
    if choice == "pick":
        while True:
            idx = input(f"Which station index or name?: ")
            try:
                station = stations[int(idx)]
            except ValueError:
                station = idx

            cw = CompWave(dsfid=dsfid, station=station, min_period=min_period,
                         max_period=max_period)
            cw.plot(m_init, m_final, show=True, save=False, fontsize=14)

    elif choice == "all":
        for station in stations:
            print(station)
            try:
                cw = CompWave(dsfid=dsfid, station=station, 
                              min_period=min_period, max_period=max_period)
                cw.plot(m_init, m_final, show=True, save=f"{station}.png", 
                        fontsize=12)
                plt.close()
            except Exception as e:
                print(e)
                pass


if __name__ == "__main__":
    plot_from_working_dir()
