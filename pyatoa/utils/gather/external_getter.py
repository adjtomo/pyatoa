#!/usr/bin/env python3
"""
Wrapper class used to quickly interact with Obspy FDSN Client, based on Config

Attributes inhereted by the Gatherer class.
"""
from pyatoa import logger

from obspy.clients.fdsn.header import FDSNNoDataException


class ExternalGetter:
    def event_get(self):
        """
        return event information parameters pertaining to a given event id
        if an event id is given

        :rtype event: obspy.core.event.Event
        :return event: event object
        """
        event = None
        logger.debug(f"fetching event from {self.client}")
        if self.event_id is not None:
            try:
                event = self.Client.get_events(eventid=self.config.event_id)[0]
                self.origintime = event.origins[0].time
            except FDSNNoDataException:
                logger.warning(f"no event found for {self.config.event_id} "
                               f"from {self.config.client}")
                event = None

        if self.origintime and event is None:
            try:
                event = self.Client.get_events(starttime=self.origintime,
                                               endtime=self.origintime)
            except FDSNNoDataException:
                logger.warning(
                    f"no event found for origin time {self.origintime}"
                    f"from {self.config.client}"
                )
            if event is not None and len(event) > 1:
                logger.warning(f"{len(event)} events found, expected only 1,"
                               f"manual event gathering may be required."
                               )

        return event

    def station_get(self, station_code, level='response'):
        """
        return station information with level dependent on call, defaults to
        retrieving response information

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :type level: str
        :param level: level argument to be passed to obspy
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        """
        logger.debug(f"fetching station from {self.client}")
        net, sta, loc, cha = station_code.split('.')
        return self.Client.get_stations(
            network=net, station=sta, location=loc, channel=cha,
            starttime=self.origintime-self.config.start_pad,
            endtime=self.origintime+self.config.end_pad, level=level
        )
 
    def waveform_get(self, station_code):
        """
        Call for obspy to download data. For some reason obspy can return traces
        with varying sample lengths, so 10 second cushion and then trim after
        retrieval to make sure traces are the same length

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stream
        """
        logger.debug(f"fetching observations from {self.client}")

        net, sta, loc, cha = station_code.split('.')
        st = self.Client.get_waveforms(
            network=net, station=sta, location=loc, channel=cha,
            starttime=self.origintime - (self.config.start_pad + 10),
            endtime=self.origintime + (self.config.end_pad + 10)
        )
        # Sometimes FDSN queries return improperly cut start and end times, so
        # we retrieve +/-10 seconds and then cut down
        st.trim(starttime=self.origintime - self.config.start_pad,
                endtime=self.origintime + self.config.end_pad)

        logger.debug(f"stream got external {station_code}")
        return st

    def get_all(self, station_code):
        """
        Convenience function for retrieving station, inventory and receiver

        :type station_code: str
        :param station_code: Station code following SEED naming convention.
            This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
            L=location, C=channel). Allows for wildcard naming. By default
            the pyatoa workflow wants three orthogonal components in the N/E/Z
            coordinate system. Example station code: NZ.OPRZ.10.HH?
        :rtype event: obspy.core.event.Event
        :return event: event object
        :rtype: obspy.core.inventory.Inventory
        :return: inventory containing relevant network and stations
        :rtype stream: obspy.core.stream.Stream
        :return stream: waveform contained in a stream
        """
        event = self.event_get()
        inv = self.station_get(station_code)
        st = self.waveform_get(station_code)
        return event, inv, st


def get_gcmt_moment_tensor(origintime, magnitude, path=None,
                           time_wiggle_sec=120, magnitude_wiggle=0.5):
    """
    Query GCMT moment tensor catalog for moment tensor components

    :type origintime: UTCDateTime or str
    :param origintime: event origin time
    :type magnitude: float
    :param magnitude: centroid moment magnitude for event lookup
    :type path: str
    :param path: path to the gcmt moment tensor files, separated by year
    :type time_wiggle_sec: int
    :param time_wiggle_sec: padding on catalog filtering criteria realted to
        event origin time
    :type mag_wiggle: float
    :param mag_wiggle: padding on catalog filter for magnitude
    :rtype event: obspy.core.event.Event
    :return event: event object for given earthquake
    """
    from urllib.error import HTTPError
    from obspy import UTCDateTime, read_events

    if not isinstance(origintime, UTCDateTime):
        datetime = UTCDateTime(origintime)

    # Determine filename using datetime properties
    month = origintime.strftime('%b').lower()  # e.g. 'jul'
    year_short = origintime.strftime('%y')  # e.g. '19'
    year_long = origintime.strftime('%Y')  # e.g. '2019'

    fid = f"{month}{year_short}.ndk"
    logger.info("querying GCMT database for moment tensor")
    try:
        cat = read_events(
            "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
            f"catalog/NEW_MONTHLY/{year_long}/{fid}"
        )
    except HTTPError:
        cat = read_events(
            "http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
            "catalog/NEW_QUICK/qcmt.ndk"
        )

    # GCMT catalogs contain all events for a span of time
    # filter catalogs using ObsPy to find events with our specifications.
    # Magnitudes and origintimes are not always in agreement between agents
    # So allow fro some wiggle room
    cat_filt = cat.filter(f"time > {str(origintime - time_wiggle_sec)}",
                          f"time < {str(origintime + time_wiggle_sec)}",
                          f"magnitude >= {magnitude - magnitude_wiggle}",
                          f"magnitude <= {magnitude + mag_wiggle}",
                          )
    # Filtering may remove all events from catalog, return multiple events, or
    # may return the event of choice
    if not len(cat_filt):
        logger.info(f"no gcmt event found for {datetime} and M{magnitude}")
        raise FileNotFoundError("No events found")
    elif len(cat_filt) > 1:
        logger.info(f"multiple events found for {datetime} and M{magnitude}")
        print(f"{len(cat_filt)} events found, choosing first")
        return cat_filt[0]
    else:
        logger.info("gcmt event found matching criteria")
        return cat_filt[0]

