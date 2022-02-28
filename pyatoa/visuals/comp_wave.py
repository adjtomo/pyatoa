#!/usr/bin/env python3
"""
A stripped down (arguably better) version of the Waveform Improvement class,
used to simply compare two synthetic waveforms from a given PyASDF DataSet.
"""
import os
import sys
from glob import glob
import traceback
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pyatoa import Manager
from pyasdf import ASDFDataSet
from pyatoa.core.config import set_pyflex_config
from pyatoa.visuals.wave_maker import format_axis
from pyatoa.visuals.map_maker import MapMaker
from pyatoa.utils.calculate import vrl
from pyatoa.utils.process import match_npts
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
        self.windows = None
        self.init_windows = None
        self.final_windows = None

        # Initialize empty attributes to be filled
        self._ds = None
        self._inv = None
        self._event = None
        self._st_obs = None
        self._st_init = None
        self._st_final = None
        self._m_init = None
        self._m_final = None

    def _gather_model_from_dataset(self, dsfid, model=None, init_or_final=None,
                                   save_windows=False):
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
            mgmt.config.save_to_ds = False

            # !!! NZ Temp network skip remove response
            net, sta = self.station.split(".")
            remove_response = not bool(net in ["ZX", "Z8"])
                
            # Overwrite the filter corners stored in the dataset
            mgmt.config.min_period = self.min_period
            mgmt.config.max_period = self.max_period
            mgmt.standardize().preprocess(remove_response=remove_response)
            # See if windows can be picked for the given load setup
            if save_windows:
                pf_cfg = "nznorth_6-30s"
                mgmt.config.pyflex_preset = pf_cfg
                mgmt.config.pyflex_config, _ = set_pyflex_config(
                                                         mgmt.config.min_period,
                                                         mgmt.config.max_period,
                                                         pf_cfg)
                mgmt.window()
                print(f"{mgmt.stats.nwin} windows for {init_or_final}")
               
                if save_windows == "init":
                    self.init_windows = mgmt.windows
                elif save_windows == "final":
                    self.final_windows = mgmt.windows
                else: 
                    # Any check that returns no windows will set this False
                    self.windows = mgmt.windows

            # Store data in class attributes, obs waveforms will be the same
            if self._st_obs is None:
                self._st_obs = mgmt.st_obs
            if self._inv is None:
                self._inv = mgmt.inv
            if self._event is None:
                self._event = mgmt.event
            setattr(self, f"_st_{init_or_final}", mgmt.st_syn)
            setattr(self, f"_m_{init_or_final}", model)

    def gather(self, m_init=None, m_final=None, save_windows=False):
        """
        Gather data from the correct dataset. If no m_init or m_final given,
        will gather the first and last models

        :type m_init: str
        :param m_init: initial iteration/step, e.g. 'i01/s00'
        :type m_final: str
        :param m_final: final iteration/step
        """
        self._gather_model_from_dataset(dsfid=self.dsfid, model=m_init,
                                        init_or_final="init", 
                                        # save_windows=bool(save_windows=="init"))
                                        save_windows="init")

        # Default to dsfid if separate final dataset not provided
        self._gather_model_from_dataset(dsfid=self.dsfid_final or self.dsfid,
                                        model=m_final, init_or_final="final",
                                        # save_windows=bool(
                                        #                 save_windows=="final"))
                                        save_windows="final")

    def calculate_vrl(self, init_or_final):
        """
        Caclulate the logarithmic variance reduction to look at how waveforms
        imrpove from m_init to m_final. Following Eq. 8 of Tape et al. (2010).
        """
        st_d = self._st_obs.copy()
        st_1 = self._st_init.copy()
        st_2 = self._st_final.copy()

        # Final synthetics dictate samp. rate because it will be higher
        if init_or_final == "init":
            sr = st_1[0].stats.sampling_rate
            dt = st_1[0].stats.delta
        elif init_or_final == "final":
            sr = st_2[0].stats.sampling_rate
            dt = st_2[0].stats.delta

        st_d.resample(sr)
        st_1.resample(sr)
        st_2.resample(sr)

        with open("vrl.txt", "a+") as f:
            for tr_d in st_d:
                comp = tr_d.stats.channel[-1] 
                tr_1 = st_1.select(component=comp)[0]
                tr_2 = st_2.select(component=comp)[0]
                # If no windows, calculate VRL for the entire trace
                if self.windows is None:
                    vrl_ = vrl(d=tr_d.data, s1=tr_1.data, s2=tr_2.data)
                    f.write(
                       f"{self.event_id}, {self.station}, {comp}, {vrl_:.4f}\n")
                # If windows, calculate VRL for each window
                else:
                    if comp in self.windows:
                        for i, win in enumerate(self.windows[comp]):
                            # Delta should match final synthetics
                            assert(dt == win.dt), "delta dont match"
                            vrl_ = vrl(d=tr_d.data[win.left:win.right], 
                                       s1=tr_1.data[win.left:win.right],
                                       s2=tr_2.data[win.left:win.right])
                            f.write(f"{self.event_id}, {self.station}, {comp}, "
                                    f"{i}, {vrl_:.4f}\n") 
                        

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
        dpi = kwargs.get("dpi", 150)
        figsize = kwargs.get("figsize", (1200 / dpi, (240 * nrows) / dpi))
        fontsize = kwargs.get("fontsize", 12)
        axis_linewidth = kwargs.get("axis_linewidth", 2)
        xticks = kwargs.get("xticks", True)
        yticks = kwargs.get("yticks", True)
        minor_ticks = kwargs.get("minor_ticks", False)

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
                
                left = bool(yticks)
                bottom = bool(xticks) 

                ax.tick_params(which='major', direction='in', top=False,
                               right=False, left=left, bottom=bottom, 
                               labelsize=fontsize,
                               length=5, width=2*axis_linewidth/3)
                if minor_ticks:
                    ax.minorticks_on()
                    ax.tick_params(which='minor', direction='in', length=3, 
                                   top=False, bottom=bottom, right=False, 
                                   left=left, width=2*axis_linewidth/3)
                if col == 0:
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                    # Make sure the scientific notation has the same4fontsize
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
            # for col in row[1:]:
            for col in row[:]:
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

    def plot(self, component_list=None, show=True, save=False, **kwargs):
        """
        Plot waveforms iterative based on model updates

        :type show: bool
        :param show: Show the plot or do not
        :type save: str
        :param save: if given, save the figure to this path
        """
        linewidth = kwargs.get("linewidth", 1.3)
        syn_linewidth = kwargs.get("syn_linewidth", 1.7)
        fontsize = kwargs.get("fontsize", 12)
        xlim = kwargs.get("xlim", None)
        percent_over = kwargs.get("percent_over", 0.125)
        color_init = kwargs.get("color_init", "orangered")
        color_final = kwargs.get("color_final", "mediumorchid")
        title = kwargs.get("title", True)
        labels = kwargs.get("labels", True)

        if component_list is None:
            component_list = [_.stats.channel[-1] for _ in self._st_obs]

        # One row per component, one column for init and final each
        f, axes = self.setup_plot(nrows=len(component_list), ncols=2, **kwargs)

        # Plot each component in a different column
        for row, comp in enumerate(component_list):
            obs = self._st_obs.select(component=comp)[0]
            syn_init = self._st_init.select(component=comp)[0]
            syn_final = self._st_final.select(component=comp)[0]

            # Plot init synthetics to the left, final synthetics to the right
            axes[row][0].plot(syn_init.times(), syn_init.data, color_init,
                              zorder=10, label="INIT", linewidth=syn_linewidth
                              )
            axes[row][1].plot(syn_final.times(), syn_final.data, color_final,
                              zorder=10, label="FINAL", linewidth=syn_linewidth
                              )
            # Plot obs in both columns
            for col in [0, 1]:
                axes[row][col].plot(obs.times(), obs.data, "k", zorder=11,
                                    label="OBS", linewidth=linewidth)

            # Format the axes for a standardized look
            for col in range(len(axes[row])):
                format_axis(axes[row][col], percent_over=percent_over)

            # Component in the y-label
            if labels:
                axes[row][0].set_ylabel(comp.upper(), rotation="horizontal",
                                        ha="left", va="center")

            # Annotation in coriner of fig
            if len(component_list) == 1:
                plt.suptitle(f"{self.event_id} {self.station} {comp}",
                             fontsize=fontsize)
            f.text(0.13, 0.765, "M00", ha="left", fontsize=fontsize)
            f.text(0.525, 0.765, "M28", ha="left", fontsize=fontsize)

            # Set ylim based on ymax from either waveform
            max_yvals = []
            for data in [obs.data, syn_init.data, syn_final.data]:
                max_yvals.append(max(abs(data)))
            max_yval = max(max_yvals) * 1.2

            ylim = [-1 * max_yval, max_yval]
            axes[row][0].set_ylim(ylim)

            # Plot different windows for init and final model, if available
            if self.init_windows and self.final_windows:
                for i, windows in enumerate([self.init_windows[comp], 
                                             self.final_windows[comp]]):
                    ymin, ymax = axes[row][i].get_ylim()
                    for win in windows:
                        tleft = win.left * win.dt 
                        tright = win.right * win.dt
                        axes[row][i].axvline(x=tleft, c="k", lw=1.5)
                        axes[row][i].axvline(x=tright, c="k", lw=1.5)
                        axes[row][i].add_patch(
                            Rectangle(xy=(tleft, ymin), 
                                      width=tright - tleft,
                                      height=(ymax + np.abs(ymin)),
                                      fc="orange", ec="k",
                                      alpha=0.25,
                                      zorder=10
                                      ))

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
        if title:
            plt.suptitle(f"Waveform Improvement | "
                         f"{self.event_id} {self.station} "
                         f"[{self.min_period}, {self.max_period}]s",
                         fontsize=fontsize)

        if labels:
            # Common X and Y labels, manually decided values
            # f.text(0.375, 0.05, "Time [s]", ha="center", fontsize=fontsize)
            f.text(0.5, 0.025, "Time [s]", ha="center", fontsize=fontsize)
            f.text(0.08, 0.575, "Displacement [m]", ha="center", rotation=90,
                   fontsize=fontsize)


        # Title each of the models
        if title:
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
        gs = mpl.gridspec.GridSpec(1, 2, wspace=1.0, hspace=0.,
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


def main(event_id=None, station=None, component=None, xmin=None, xmax=None, 
         cfg="plot"):
    """
    Main call script to choose event and station based on what's available
    """
    if cfg == "calc":
        # VRL CALC PARAMETER CHOICE
        choice = "all"  # pick or all
        min_period = 6
        max_period = 30
        m_init = None
        m_final = "i11/s03"
        dsfid = f"aspen/{event_id}.h5"
        dsfid_final = f"birch/{event_id}.h5"
        show = False
        calc_vrl = True
        save_win = "final"
        plot = False
        plot_with_map = False
        xlim = None
        if component:
            component_list = [component]
        else:
            component_list = ["Z", "N", "E"]
    elif cfg == "plot":
        choice = "all"  # pick or all
        min_period = 6
        max_period = 30
        m_init = None
        m_final = "i11/s03"
        dsfid = f"aspen/{event_id}.h5"
        dsfid_final = f"birch/{event_id}.h5"
        show = False
        calc_vrl = False
        save_win = True
        plot = True
        plot_with_map = False   
        if xmin is not None:
            xlim = [float(xmin), float(xmax)]
        else:
            xlim = None
        if component:
            component_list = [component]
        else:
            component_list = ["Z", "N", "E"]

    kwargs = {"title": False, "labels": False, "xticks": True, "yticks": False,
              "minor_ticks": False}
    # =========================================================================

    # Get station information prior to plotting
    assert(os.path.exists(dsfid)), f"{dsfid} does not exist"
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
        # if os.path.exists(fid_out):
        #     print(f"{fid_out} exists")
        #     continue
        print(sta)
        try:
            cw = CompWave(dsfid=dsfid, dsfid_final=dsfid_final,
                          station=sta, min_period=min_period,
                          max_period=max_period)
            cw.gather(m_init, m_final, save_win)

            if calc_vrl:
                cw.calculate_vrl(save_win)
            if plot:
                if plot_with_map:
                    cw.plot_with_map(show=show, save=fid_out, xlim=xlim)
                else:
                    cw.plot(component_list=component_list, show=show, 
                            save=fid_out, xlim=xlim, **kwargs)
                plt.close("all")
        except Exception as e:
            traceback.print_exc()
            pass

def main_plot_specific():
    """
    """
    vals = [
["2013p617227", "NZ.TOZ", "Z", 47.5, 290],
["2014p952799", "NZ.NTVZ", "N", 10, 290],
["2016p105478", "NZ.PUZ", "Z", 90, 290],
["2016p881118", "NZ.MWZ", "E", 40, 275],
["2018p465580", "NZ.KHEZ", "E", 25, 275],
["2019p738432", "NZ.KHZ", "Z", 90, 290],
["2019p754447", "NZ.HIZ", "Z", 10, 290],
["2019p927023", "NZ.VRZ", "Z", 45, 275],]
    for val in vals:
        # event, station, comp = val
        # main(event, station, comp, cfg="plot")
        try:
            event, station, comp, xmin, xmax = val
            main(event, station, comp, xmin, xmax, "plot")
        except ValueError:
            continue


if __name__ == "__main__":
    main_plot_specific()
    a=1/0
    if len(sys.argv) > 1:
        main(*sys.argv[1:])
    else:
        for fid in sorted(glob("./aspen/*.h5")):
            event_id = os.path.splitext(os.path.basename(fid))[0]
            print(event_id)
            main(event_id, None)
