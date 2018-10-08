"""
Pyasdf dataset information parser

Pyatoa's main workflow saves all relevant inforomation into pyasdf datasets.
Datasets are a good way to organize information by tags and nested structures,
but they can be difficult to navigate and to ascertain information from
The diver is a wrapper class which provies convenience functions for returning
statistics, information and counts of information inside a pyasdf datset
"""
import os
import copy
import sys
import pyasdf
import numpy as np


class Diver:
    """
    Information parsing class
    """
    def __init__(self, ds):
        """
        Instantiate with a pyasdf dataset

        :type ds: pyasdf.asdf_data_set.ASDFDataSet
        :param ds: dataset to deep dive into
        """
        self.ds = ds

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

    def write_adj_src_to_specfem(self, buf="./", model="m00"):
        """
        write adjoint source from asdf data format to specfem input. when
        writing to asdf data format, the adjoint source was time reversed for
        convenient plotting, so this function needs to re-reverse, i.e. put
        the data in the same format that pyadjoint outputted it
        :param buffer:
        :return:
        """
        shortcut = self.ds.auxiliary_data.AdjointSources[model]
        for adj_src in shortcut:
            file_template = "{id}.adj".format(
                shortcut.adj_src.path.replace('_', '.')
            )
            outfile = os.path.join(buf, file_template)
            np.savetxt(outfile, adj_src.data.value)

    def collect_misfits(self,model="m00"):
        """for each station, collect the misfit value
        """
        misfit_values = {}
        for AS in self.adj_srcs:
            parm = self.aux.AdjointSource[AS].parameters
            channel_id = '{}.{}'.format(parm["station_id"],parm["component"])
            misfit_values[channel_id] = parm["misfit_value"]

        self.misfit_values = misfit_values

    def get_srcrcv_information(self):
        """calculate source receiver information for each pair
        """
        from obspy.geodetics import gps2dist_azimuth

        event_lat,event_lon = (self.ds.events[0].origins[0].latitude,
                               self.ds.events[0].origins[0].longitude)

        srcrcvdict = {}
        for sta in self.stations:
            coordict = self.ds.waveforms[sta].coordinates
            GCDist,Az,BAz = gps2dist_azimuth(lat1=event_lat,lon1=event_lon,
                                             lat2=coordict["latitude"],
                                             lon2=coordict["longitude"])
            srcrcvdict[sta] = {"great_circle_distance":GCDist,
                               "azimuth":Az,
                               "backazimuth":BAz}

        self.source_receiver_info = srcrcvdict
    
    def collect_cc_time_shifts(self):
        """gather the values of time shifts from cross correlations provided
        in the misfit windows of pyflex
        """
        if not self.aux:
            self.populate()
        cc_time_shifts = {}
        for net in self.misfit_windows:
            for miswin in self.aux.MisfitWindows[net]:
                cc_time_shifts[miswin.path] = \
                                        miswin.parameters['cc_shift_in_seconds']
        
        self.cc_time_shifts = cc_time_shifts
        
    def cc_time_shift_histogram(self,binsize=0.1):
        """plot histograms with cross correlation time shifts
        """
        if not self.cc_time_shifts:
            self.collect_cc_time_shifts()
            
        from dataDepicter import Depicter
        depicter = Depicter(tack=self)
        depicter.plot_cc_time_shift_histogram(binsize=binsize)

    def misfit_histogram(self,m0=None,m_a=None,binsize=0.1):
        """plot a histogram of misfits for models m0 and m_a (if available)
        default will be initial model compared to final model
        """
        if not self.misfit_values:
            self.collect_misfits()
            
        from dataDepicter import Depicter
        depicter = Depicter(tack=self)
        depicter.plot_misfit_histogram(m0=m0, m_a=m_a, binsize=binsize)

