"""
Pyasdf dataset information parser

Pyatoa's main workflow saves all relevant inforomation into pyasdf datasets.
Datasets are a good way to organize information by tags and nested structures,
but they can be difficult to navigate and to ascertain information from
The diver is a wrapper class which provies convenience functions for returning
statistics, information and counts of information inside a pyasdf datset
"""
import json


class Diver:
    """
    Information parsing class
    """
    def __init__(self):
        """
        Instantiate with a pyasdf dataset

        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset to deep dive into
        """
        self.ds = ds
        self.config = config
        self.misfit_values = None
        self.source_receiver_info = None
        self.cc_time_shifts = {}

    def __str__(self):
        """
        replace print() output
        """
        template = ("Deep dive for for {pat}\n"
                    "\t{nos:<4} stations\n"
                    "\t{now:<4} misfit_windows\n"
                    "\t{noa:<4} adj_srcs\n")
        return template.format(pat=self.path,
                               nos=len(self.stations),
                               now=len(self.misfit_windows),
                               noa=len(self.adj_srcs)
                               )

    def suggest(json_fid):
        """
        Takes the misfit json fid that should be outputted by the Pyatoa
        workflow, analyze the number of misfit windows and total misfits
        and provide some suggestions to the User about what the events may look
        like to avoid having to visualize kernels or waveforms
        :param fid:
        :return:
        """
        with open(json_fid, "r") as f:
            mft = json.load(f)
        for model in mft.keys():
            print("MODEL {}".format(model))
            for step in mft[model].keys():
                print("\tSTEP {}".format(step))
                extra_keys = ["misfit", "adjsrcs", "windows"]
                number_events = len(mft[model][step].keys()) - len(extra_keys)
                total_misfit = mft[model][step]["misfit"]
                for key in mft[model][step].keys():
                    if key not in extra_keys:
                        percent_misfit = (
                            (mft[model][step][key]["misfit"]/number_events) /
                            total_misfit
                        )
                        print(mft[model][step][key]["windows"] * percent_misfit)
                        # print(
                        #     "\t\t{ev}, {pc:.2f}% misfit w/ {wi} windows".format(
                        #         ev=key, pc=percent_misfit,
                        #         wi=mft[model][step][key]["windows"]
                        #     ))


    def count_windows(self):
        """
        Figure out which stations contain which windows, return a dictionary
        which lists available components.
        """
        stations = []
        for model_number in self.ds.auxiliary_data.MisfitWindows.list():
            for window in self.ds.auxiliary_data.MisfitWindows[model_number]:
                stations.append(
                    self.ds.auxiliary_data.MisfitWindows[model_number]
                    [window].parameters['channel_id']
                )
        uniqueid = set(stations)
        counted_windows = {}
        for id_ in uniqueid:
            counted_windows[id_] = stations.count(id_)
        return counted_windows

    def get_srcrcv_information(self):
        """calculate source receiver information for each pair
        """
        from obspy.geodetics import gps2dist_azimuth

        event_lat, event_lon = (self.ds.events[0].origins[0].latitude,
                                self.ds.events[0].origins[0].longitude)
        srcrcvdict = {}
        for sta in self.stations:
            coordict = self.ds.waveforms[sta].coordinates
            gcdist, az, baz = gps2dist_azimuth(lat1=event_lat, lon1=event_lon,
                                               lat2=coordict["latitude"],
                                               lon2=coordict["longitude"])
            srcrcvdict[sta] = {"great_circle_distance": gcdist, "azimuth": az,
                               "backazimuth": baz}
        self.source_receiver_info = srcrcvdict

    def collect_misfits(self, model="m00"):
        """
        for each station, collect the misfit value
        :param model: model number
        """
        misfit_values = {}
        for AS in self.ds.auxiliary_data.AdjointSources[model].list():
            parm = self.ds.auxiliary_data.AdjointSources[model][AS].parameters
            channel_id = '{}.{}'.format(parm["station_id"],parm["component"])
            misfit_values[channel_id] = parm["misfit_value"]
        self.misfit_values = misfit_values

    def collect_cc_time_shifts(self, model="m00"):
        """gather the values of time shifts from cross correlations provided
        in the misfit windows of pyflex
        """
        cc_time_shifts = {}
        for path in self.ds.auxiliary_data.MisfitWindows[model].list():
            miswin = self.ds.auxiliary_data.MisfitWindows[model][path]
            cc_time_shifts[miswin.path] = \
                                        miswin.parameters['cc_shift_in_seconds']
        self.cc_time_shifts[model] = cc_time_shifts
        
    def cc_time_shift_histogram(self, binsize=0.1):
        """plot histograms with cross correlation time shifts
        """
        if not self.cc_time_shifts:
            self.collect_cc_time_shifts()
        from pyatoa.visuals.plot_statistics import plot_cc_time_shift_histogram
        plot_cc_time_shift_histogram(self.cc_time_shifts, self.config,
                                     binsize=binsize)

    def misfit_histogram(self, binsize=0.1):
        """plot a histogram of misfits for models m0 and m_a (if available)
        default will be initial model compared to final model
        """
        if not self.misfit_values:
            self.collect_misfits()
        from pyatoa.visuals.plot_statistics import plot_misfit_histogram
        plot_misfit_histogram(self.misfit_values, self.config, binsize=binsize)


