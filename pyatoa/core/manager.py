#!/usr/bin/env python3
"""
Main workflow components of Pyatoa.
Manager is the central workflow control object. It calls on mid and low level
classes to gather data, and then runs these through Pyflex for misfit window
identification, and then into Pyadjoint for misfit quantification. Config class
required to set the necessary parameters
"""
import warnings

import obspy
import pyflex
import pyadjoint
import numpy as np
from os.path import basename
from obspy.signal.filter import envelope

from pyatoa import logger
from pyatoa.core.gatherer import Gatherer

from pyatoa.utils.asdf.additions import write_adj_src_to_asdf
from pyatoa.utils.operations.source_receiver import gcd_and_baz,\
    seismogram_length
from pyatoa.utils.operations.formatting import create_window_dictionary, \
     channel_codes
from pyatoa.utils.processing.preprocess import preproc, trimstreams
from pyatoa.utils.processing.synpreprocess import stf_convolve_gaussian
from pyatoa.utils.configurations.external import set_pyflex_station_event
from pyatoa.utils.visuals.mapping import generate_map
from pyatoa.utils.visuals.waveforms import window_maker


class Manager:
    """
    Core object within Pyatoa.

    Workflow management function that internally calls on all other objects
    within the package in order to gather, process and analyze waveform data.
    """
    def __init__(self, config, ds=None, empty=True, station_code=None,
                 event=None, st_obs=None, st_syn=None, inv=None, windows=None,
                 staltas=None, adj_srcs=None):
        """
        If no pyasdf dataset is given in the initiation of the Manager, all
        data fetching will happen via given pathways in the config file,
        or through external getting via FDSN pathways

        :type config: pyatoa.core.config.Config
        :param config: configuration object that contains necessary parameters
            to run through the Pyatoa workflow
        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: ASDF data set from which to read and write data
        :type empty: bool
        :param empty: Do not instantiate Gatherer and look for event.
            Useful for when User provides own event object
        :type station_code: str
        :param station_code: station code for data gather, e.g. 'NZ.BFZ.10.HH?'
        :type event: obspy.core.event.Event
        :param event: An event object containing relevant earthquake information
        :type st_obs: obspy.core.stream.Stream
        :param st_obs: Stream object containing waveforms of observations
        :type st_syn: obspy.core.stream.Stream
        :param st_syn: Stream object containing waveforms of observations
        :type inv: obspy.core.inventory.Inventory
        :param inv: Inventory that should only contain the station of interest,
            it's relevant channels, and response information
        :type windows: dict of pyflex.Window objects
        :param windows: misfit windows calculated by Pyflex, stored in a
            dictionary based on component naming
        :type adj_srcs: dict of pyadjoint.AdjointSource objects
        :param adj_srcs: adjoint source waveforms stored in dictionaries

        """
        self.config = config
        self.ds = ds
        self.gatherer = None

        # Data required for each workflow
        self.station_code = station_code
        self.st_obs = st_obs
        self.st_syn = st_syn
        self.inv = inv

        # Data produced by the workflow
        self.windows = windows
        self.staltas = staltas
        self.adj_srcs = adj_srcs

        # Internally used statistics
        self.time_offset_sec = None
        self.num_windows = None
        self.total_misfit = None
        self.half_dur = None

        # Flags to show status of workflow
        self.dataset_id = None
        self.event_name = None
        self.st_obs_num = 0
        self.st_syn_num = 0
        self.inv_name = None
        self.half_dur_flag = None
        self.obs_process_flag = False
        self.syn_process_flag = False
        self.syn_shift_flag = False
        self.pyflex_windows = 0
        self.pyadjoint_misfit = 0

        # Launch the Manager class
        if not empty and not event:
            self.launch()
        else:
            self.event = event

    def __str__(self):
        """
        Print statement shows available information inside the workflow.
        """
        self._check()
        return ("DATA\n"
                "\tdataset (ds):                 {dataset}\n"
                "\tevent:                        {event}\n"
                "\tinventory (inv):              {inventory}\n"
                "\tobserved data (st_obs):       {obsstream}\n"
                "\tsynthetic data (st_syn):      {synstream}\n"
                "WORKFLOW\n"
                "\tobs data preprocessed:        {obsproc}\n"
                "\tsyn data preprocessed:        {synproc}\n"
                "\thalf duration (half_dur):     {synshift}\n"
                "\tpyflex (num_windows):         {pyflex}\n"
                "\tpyadjoint (total_misfit):     {pyadjoint}\n"
                ).format(dataset=self.dataset_id,
                         event=self.event_name,
                         obsstream=self.st_obs_num,
                         synstream=self.st_syn_num,
                         inventory=self.inv_name,
                         obsproc=self.obs_process_flag,
                         synproc=self.syn_process_flag,
                         synshift=self.half_dur_flag,
                         pyflex=self.pyflex_windows,
                         pyadjoint=self.pyadjoint_misfit
                         )

    def _check(self):
        """
        Update flag information for the User to know location in the workflow
        """
        # Give dataset filename if available
        if self.ds:
            self.dataset_id = basename(self.ds.filename)

        # Check to see if event is an Event object
        if isinstance(self.event, obspy.core.event.Event):
            self.event_name = self.event.resource_id
        else:
            self.event_name = None

        # Make sure observed waveforms are Stream objects,
        # Check if the observed waveforms have been preprocessed
        if isinstance(self.st_obs, obspy.Stream):
            self.st_obs_num = len(self.st_obs)
            self.obs_process_flag = (
                    hasattr(self.st_obs[0].stats, "processing") and
                    len(self.st_obs[0].stats.processing) >= 3
            )
        else:
            self.st_obs_num, self.obs_process_flag = 0, False

        # Make sure synthetic waveforms are Stream objects
        # Check if the synthetic waveforms have been preprocessed
        if isinstance(self.st_syn, obspy.Stream):
            self.st_syn_num = len(self.st_syn)
            self.syn_process_flag = (
                    hasattr(self.st_syn[0].stats, "processing") and
                    len(self.st_syn[0].stats.processing) >= 3
            )
        else:
            self.st_syn_num, self.syn_process_flag = 0, False

        # Check if preprocessing found a source time function
        if self.half_dur is not None:
            self.half_dur_flag = "{:.3f}".format(self.half_dur)

        # Check to see if inv is an Inventory object
        if isinstance(self.inv, obspy.Inventory):
            self.inv_name = "{net}.{sta}".format(net=self.inv[0].code,
                                                 sta=self.inv[0][0].code)
        else:
            self.inv_name = None

        # Check if Pyflex has been run, if so return the number of windows made
        if isinstance(self.windows, dict):
            self.pyflex_windows = self.num_windows

        # Check if Pyadjoint has been run, if so return the total misfit
        if isinstance(self.adj_srcs, dict):
            self.pyadjoint_misfit = "{:.3f}".format(self.total_misfit)

    def reset(self, hard_reset=False):
        """
        Convenience function to delete all data in the Manager, and hence start
        the workflow from the start without losing your event or gatherer.

        Soft reset retains event information so that another Station can be
        gathered for the same event

        :type hard_reset: bool
        :param hard_reset: hard or soft, soft reset does not re-instantiate
            gatherer class, and leaves the same event. Useful for short workflow
        """
        self.station_code = None
        self.st_obs = None
        self.st_syn = None
        self.inv = None
        self.windows = None
        self.staltas = None
        self.adj_srcs = None
        self.time_offset_sec = None
        self.num_windows = None
        self.total_misfit = None
        self.half_dur = None

        if hard_reset:
            self.launch(reset=True)

    @property
    def st(self):
        if isinstance(self.st_syn, obspy.Stream) and \
                isinstance(self.st_obs, obspy.Stream):
            return self.st_syn + self.st_obs
        elif isinstance(self.st_syn, obspy.Stream) and \
                not isinstance(self.st_obs, obspy.Stream):
            return self.st_syn
        elif isinstance(self.st_obs, obspy.Stream) and \
                not isinstance(self.st_syn, obspy.Stream):
            return self.st_obs
        else:
            return None

    def launch(self, reset=False, set_event=None):
        """
        Initiate the prerequisite parts of the Manager class. Populate with
        an obspy event object which is gathered from FDSN by default.
        Allow user to provide their own ObsPy event object so that the
        Gatherer does not need to query FDSN
        :type reset: bool
        :param reset: Reset the Gatherer class for a new run
        :type set_event: obspy.core.event.Event
        :param set_event: if given, will bypass gathering the event and manually
            set to the user given event
        """
        # First, launch the gatherer or reset the gatherer
        if (self.gatherer is None) or reset:
            logger.info("initiating/resetting gatherer")
            self.gatherer = Gatherer(config=self.config, ds=self.ds)

        # Second, get event information
        if set_event:
            self.gatherer.event = set_event
            self.event = self.gatherer.event
        # If no event given by User, start searching for the event
        else:
            # If the gatherer already has an event, have Manager copy it
            if self.gatherer.event:
                self.event = self.gatherer.event
            # If this is the first time launching the gatherer, gather event
            else:
                if self.config.event_id:
                    self.event = self.gatherer.gather_event()

    def gather_data(self, station_code):
        """
        Launch a gatherer object and gather event, station and waveform
        information given a station code. Fills the manager based on information
        most likely to be available (we expect an event to be available more
        often than waveform data).
        Catches general exceptions along the way, stops gathering if errors.

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        """
        if self.gatherer is None:
            self.launch()
        try:
            self.station_code = station_code
            logger.info("Gathering {station} for {event}".format(
                station=station_code, event=self.config.event_id)
            )
            logger.info("gathering station information")
            self.inv = self.gatherer.gather_station(station_code)
            logger.info("gathering observation waveforms")
            self.st_obs = self.gatherer.gather_observed(station_code)
            logger.info("gathering synthetic waveforms")
            self.st_syn = self.gatherer.gather_synthetic(station_code)
        except Exception as e:
            print(e)
            return

    def populate_dataset(self):
        """
        If gather_data is skipped altogether (e.g. User provides streams,
        inventory and event information in the Manager init), but data should
        still be saved into an ASDFDataSet for data storage, this function will
        fill that dataset in the same fashion as the data_gatherer class, and
        its internal subclasses.
        :return:
        """
        # Only populate if all requisite parts are given
        if self.event and self.inv and self.st_obs and self.st_syn:
            try:
                self.ds.add_quakeml(self.event)
            except ValueError:
                warnings.warn("Event already present, not added", UserWarning)
            try:
                self.ds.add_stationxml(self.inv)
            except TypeError:
                warnings.warn("Inv already present, not added", UserWarning)
            # PyASDF has its own warnings if observed data already present
            self.ds.add_waveforms(waveform=self.st_obs,
                                  tag=self.config.observed_tag)
            self.ds.add_waveforms(waveform=self.st_syn,
                                  tag=self.config.synthetic_tag)
        else:
            warnings.warn(
                "Manager does not contain all data, cannot fill dataset",
                UserWarning)

    def preprocess(self):
        """
        Preprocess observed and synthetic data in place on waveforms.
        External preprocess function called identical mannger for
        observation and synthetic waveforms.
        Waveforms are trimmed and sampled to the same values.
        Synthetic waveform is convolved with a source time function.
        """
        self._check()
        if not self.st_obs_num or not self.st_syn_num:
            warnings.warn("No waveform data", UserWarning)
            return
        if not isinstance(self.inv, obspy.core.inventory.Inventory):
            warnings.warn("No inventory", UserWarning)
            return

        # Process observation waveforms
        logger.info("preprocessing observation data")

        # If set in Config, rotate based on source receiver lat/lon values
        if self.config.rotate_to_rtz:
            _, baz = gcd_and_baz(event=self.event, sta=self.inv[0][0])
        else:
            baz = None

        # Adjoint sources require the same sampling_rate as the synthetics
        sampling_rate = self.st_syn[0].stats.sampling_rate

        # Run external preprocessing script
        self.st_obs = preproc(st=self.st_obs, inv=self.inv,
                              resample=sampling_rate, synthetic_unit=None,
                              back_azimuth=baz,
                              pad_length_in_seconds=self.config.zero_pad,
                              unit_output=self.config.unit_output,
                              corners=self.config.filter_corners,
                              filter_bounds=[self.config.min_period,
                                             self.config.max_period],
                              )
        
        # Run external synthetic waveform preprocesser
        logger.info("preprocessing synthetic data")
        self.st_syn = preproc(st=self.st_syn, inv=None, resample=None,
                              synthetic_unit=self.config.synthetic_unit,
                              back_azimuth=baz,
                              pad_length_in_seconds=self.config.zero_pad,
                              unit_output=self.config.unit_output,
                              corners=self.config.filter_corners,
                              filter_bounds=[self.config.min_period,
                                             self.config.max_period]
                              )
        
        # Mid-check to see if preprocessing failed
        self._check()
        if not self.obs_process_flag or not self.syn_process_flag:
            warnings.warn("Preprocessing failed", UserWarning)
            return

        # Trim observations and synthetics to the length of synthetics
        self.st_obs, self.st_syn = trimstreams(
            st_a=self.st_obs, st_b=self.st_syn, force="b")

        # Retrieve the first timestamp in the .sem? file from Specfem
        self.time_offset_sec = (self.st_syn[0].stats.starttime -
                                self.event.preferred_origin().time
                                )

        # Convolve synthetic data with a gaussian source-time-function
        try:
            moment_tensor = self.event.focal_mechanisms[0].moment_tensor
            self.half_dur = moment_tensor.source_time_function.duration / 2
            self.st_syn = stf_convolve_gaussian(
                st=self.st_syn, half_duration=self.half_dur, time_shift=False
            )
        except AttributeError:
            logger.info("moment tensor not found for event")

    def run_pyflex(self):
        """
        Call Pyflex to calculate best fitting misfit windows given observation
        and synthetic data. Return dictionaries of window objects,
        as well as STA/LTA traces. If a pyasdf dataset is present,
        save misfit windows in as auxiliary data.
        If no misfit windows are found for a given station, throw a warning
        because pyadjoint won't run.
        Pyflex configuration is given by the config as a list of values with
        the following descriptions:

        i  Standard Tuning Parameters:
        0: water level for STA/LTA (short term average/long term average)
        1: time lag acceptance level
        2: amplitude ratio acceptance level (dlna)
        3: normalized cross correlation acceptance level
        i  Fine Tuning Parameters
        4: c_0 = for rejection of internal minima
        5: c_1 = for rejection of short windows
        6: c_2 = for rejection of un-prominent windows
        7: c_3a = for rejection of multiple distinct arrivals
        8: c_3b = for rejection of multiple distinct arrivals
        9: c_4a = for curtailing windows w/ emergent starts and/or codas
        10:c_4b = for curtailing windows w/ emergent starts and/or codas
        """
        # Pre-check to see if data has already been gathered
        self._check()
        if not self.obs_process_flag or not self.syn_process_flag:
            warnings.warn("can't run Pyflex, no waveforms, or not processed")
            return

        # Create Pyflex Station and Event objects
        pf_sta, pf_ev = set_pyflex_station_event(inv=self.inv, event=self.event)

        # Empties to see if no windows were collected, windows and staltas
        # saved as dictionary objects by component name
        num_windows = 0
        windows, staltas = {}, {}
        for comp in self.config.component_list:
            try:
                # Run Pyflex to select misfit windows as Window objects
                window = pyflex.select_windows(
                    observed=self.st_obs.select(component=comp),
                    synthetic=self.st_syn.select(component=comp),
                    config=self.config.pyflex_config[1], event=pf_ev,
                    station=pf_sta,
                    )
                if window:
                    # Run Pyflex to collect STA/LTA information for plotting
                    # Only collect STA/LTA if windows are also collected
                    stalta = pyflex.stalta.sta_lta(
                        dt=self.st_syn.select(component=comp)[0].stats.delta,
                        min_period=self.config.min_period,
                        data=envelope(
                            self.st_syn.select(component=comp)[0].data)
                    )
                    staltas[comp] = stalta
                    windows[comp] = window
            except IndexError:
                window = []

            # Count windows and tell User
            num_windows += len(window)
            logger.info("{0} window(s) for comp {1}".format(len(window), comp))

        # Let the User know that no windows were found for this station
        if num_windows == 0:
            warnings.warn("Empty windows", UserWarning)

        # Store information for Pyadjoint and plotting
        self.windows = windows
        self.staltas = staltas
        self.num_windows = num_windows
    
        # If an ASDFDataSet is given, save the windows into auxiliary_data
        if self.ds is not None and (num_windows != 0):
            logger.info("saving misfit windows to PyASDF")
            for comp in windows.keys():
                for i, window in enumerate(windows[comp]):
                    tag = "{mod}/{net}_{sta}_{cmp}_{num}".format(
                        net=self.st_obs[0].stats.network,
                        sta=self.st_obs[0].stats.station,
                        cmp=comp, mod=self.config.model_number, num=i)

                    # ASDF auxiliary_data subgroups don't play nice with nested
                    # dictionaries, which the window parameters are. Format them
                    # a bit simpler for saving into the dataset
                    window_dict = create_window_dictionary(window)

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        self.ds.add_auxiliary_data(data=np.array([True]),
                                                   data_type="MisfitWindows",
                                                   parameters=window_dict,
                                                   path=tag
                                                   )
        # Let the User know the outcomes of Pyflex
        logger.info("{} window(s) total found".format(num_windows))

    def run_pyadjoint(self):
        """
        Run pyadjoint on observation and synthetic data given misfit windows
        calculated by pyflex. Method for caluculating misfit set in config,
        pyadjoint config set in external configurations. Returns a dictionary
        of adjoint sources based on component. Saves resultant dictionary
        to a pyasdf dataset if given.

        NOTE: This is not in the PyAdjoint docs, but in
        pyadjoint.calculate_adjoint_source, the window needs to be a list of
        lists, with each list containing the [left_window,right_window];
        each window argument should be given in units of time (seconds)

        NOTE: Pyatoa was developed with a dev branch of Pyadjoint, located here

        https://github.com/computational-seismology/pyadjoint/tree/dev

        Lion's version of Pyadjoint does not contain some of these functions
        """
        # Precheck to see if Pyflex has been run already
        self._check()
        if not self.pyflex_windows and self.num_windows != 0:
            warnings.warn("can't run Pyadjoint, no Pyflex windows", UserWarning)
            return

        logger.info("running Pyadjoint type {} ".format(
            self.config.pyadjoint_config[0]))

        # Iterate over given windows produced by Pyflex
        total_misfit = 0
        adjoint_sources = {}
        for key in self.windows:
            adjoint_windows = []

            # Prepare window indices to give to Pyadjoint
            for win in self.windows[key]:
                adj_win = [win.left * self.st_obs[0].stats.delta,
                           win.right * self.st_obs[0].stats.delta]
                adjoint_windows.append(adj_win)

            # Run Pyadjoint to retrieve adjoint source objects
            adj_src = pyadjoint.calculate_adjoint_source(
                adj_src_type=self.config.pyadjoint_config[0],
                config=self.config.pyadjoint_config[1],
                observed=self.st_obs.select(component=key)[0],
                synthetic=self.st_syn.select(component=key)[0],
                window=adjoint_windows, plot=False
                )
            
            # Save adjoint sources in dictionary object
            adjoint_sources[key] = adj_src
            logger.info("{misfit:.3f} misfit for comp {comp}".format(
                misfit=adj_src.misfit, comp=key)
                        )
            total_misfit += adj_src.misfit

            # If ASDFDataSet given, save Adjoint source into auxiliary data
            if self.ds:
                logger.info("saving adjoint sources {} to PyASDF".format(key))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    # The tag hardcodes an X as the second channel index
                    # to signify that these are synthetic waveforms
                    # This is required by Specfem3D
                    tag = "{mod}/{net}_{sta}_{ban}X{cmp}".format(
                        mod=self.config.model_number, net=adj_src.network,
                        sta=adj_src.station,
                        ban=channel_codes(self.st_syn[0].stats.delta),
                        cmp=adj_src.component[-1]
                    )
                    write_adj_src_to_asdf(adj_src=adj_src, ds=self.ds, tag=tag,
                                          time_offset=self.time_offset_sec)
        
        # Save adjoint source for plotting
        self.adj_srcs = adjoint_sources
        
        # Save total misfit, calucalated a la Tape (2010) Eq. 6
        self.total_misfit = 0.5 * total_misfit/self.num_windows

        # Let the User know the outcome of Pyadjoint
        logger.info("total misfit {:.3f}".format(self.total_misfit))

    def plot_wav(self, append_title='', figsize=(11.69, 8.27), dpi=100,
                 show=True, save=None, return_figure=False):
        """
        Waveform plots for all given components.
        If specific components are not given (e.g. adjoint source waveform),
        they are omitted from the final plot. Plotting should be dynamic, i.e.
        if only 2 components are present in the streams, only two subplots
        should be generated in the figure.
        :type append_title: str
        :param append_title: any string to append the title of the plot
        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :type figsize: tuple of floats
        :param figsize: length and width of the figure
        :type dpi: int
        :param dpi: dots per inch of the figure
        :type return_figure: bool
        :param return_figure: option to get the figure back
        """
        # Precheck for waveform data
        self._check()
        if not self.obs_process_flag and not self.syn_process_flag:
            warnings.warn("cannot plot, no/improper waveform data", UserWarning)
            return

        # Calculate the seismogram length
        length_s = seismogram_length(
            distance_km=gcd_and_baz(event=self.event, sta=self.inv[0][0])[0],
            slow_wavespeed_km_s=2, binsize=50, minimum_length=100
        )
        
        # Call on window making function to produce waveform plots
        fig_window = window_maker(
            st_obs=self.st_obs, st_syn=self.st_syn, windows=self.windows,
            staltas=self.staltas, adj_srcs=self.adj_srcs, length_s=length_s,
            time_offset=self.time_offset_sec,
            stalta_wl=self.config.pyflex_config[1].stalta_waterlevel,
            unit_output=self.config.unit_output, config=self.config,
            figsize=figsize, total_misfit=self.total_misfit,
            append_title=append_title, dpi=dpi, show=show, save=save
        )

        if return_figure:
            return fig_window

    def plot_map(self, map_corners, show_nz_faults=False,
                 annotate_names=False, color_by_network=False,
                 figsize=(8, 8.27), dpi=100, show=True, save=None):
        """
        Map plot showing a map of the given target region. All stations that
        show data availability (according to the station master list) are
        plotted as open markers. Event is plotted as a beachball if a moment
        tensor is given, station of interest highlighted, both are connected
        with a dashed line.
        Source receier information plotted in lower right hand corner of map.
        :type map_corners: dict
        :param map_corners: {lat_min, lat_max, lon_min, lon_max}
        :type annotate_names: bool
        :param annotate_names: annotate station names
        :type color_by_network: bool
        :param color_by_network: color based on network name
        :type show: bool
        :param show: show the plot once generated, defaults to False
        :type save: str
        :param save: absolute filepath and filename if figure should be saved
        :type show_nz_faults: bool
        :param show_nz_faults: plot active faults and hikurangi trench from
            internally saved coordinate files. takes extra time over simple plot
        :type figsize: tuple of floats
        :param figsize: length and width of the figure
        :type dpi: int
        :param dpi: dots per inch of the figure
        """
        self._check()

        # Warn user if no inventory is given
        if not isinstance(self.inv, obspy.Inventory):
            warnings.warn("no inventory given, plotting blank map", UserWarning)

        if not map_corners:
            map_corners = self.config.map_corners

        # Call external function to generate map
        generate_map(map_corners=map_corners, event=self.event, inv=self.inv,
                     show_nz_faults=show_nz_faults,
                     color_by_network=color_by_network,
                     annotate_names=annotate_names, show=show, figsize=figsize,
                     dpi=dpi, save=save
                     )
