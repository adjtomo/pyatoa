"""
Various scripts to grab auxiliary data such as moment tensor information,
station information, and fault information. these are all specific to the
New Zealand tomography problem, and therefore paths are hard coded
"""
import os
import csv
import requests
import warnings

from obspy import UTCDateTime, read_events
from pyatoa import logger


def grab_geonet_moment_tensor(event_id, fid=None):
    """
    Get moment tensor information from a internal csv file,
    or from an external github repository query.
    Only relevant to the new zealand tomography problem.
    Geonet moment tensors stored with a specific column format.

    :type event_id: str
    :param event_id: unique event identifier
    :type fid: str
    :param fid: absolute path to the geonet moment tensor file, if None,
        search external github repository
    :rtype moment_tensor: dict
    :return moment_tensor: dictionary created from rows of csv file
    """
    # If a csv file is given, read that, else check external
    if fid:
        with open(fid) as f:
            reader = csv.reader(f)
            tag = "internal"
    else:
        # Request and open the CSV file. Assumed that GeoNet will keep their
        # moment-tensor information in their GitHub repository
        # OK, last accessed 23.6.19
        geonet_mt_csv = (
            "https://raw.githubusercontent.com/GeoNet/data/master/"
            "moment-tensor/GeoNet_CMT_solutions.csv"
        )
        response = requests.get(geonet_mt_csv)
        if not response.ok:
            warnings.warn("Github repo request failed", UserWarning)
            return None

        # Use CSV to parse through the returned repsonse
        reader = csv.reader(response.text.splitlines(), delimiter=',')
        tag = "external"

    # Parse the CSV file
    for i, row in enumerate(reader):
        # First row contains header information
        if i == 0:
            tags = row
        # First column gives event ids
        if row[0] == event_id:
            values = []
            # Grab the relevant information from the file
            for t, v in zip(tags, row):
                if t == "Date":
                    values.append(UTCDateTime(v))
                elif t == "PublicID":
                    values.append(v)
                else:
                    values.append(float(v))

            moment_tensor = dict(zip(tags, values))
            logger.info("geonet moment tensor {} for event: {}".format(
                tag, event_id)
            )
            return moment_tensor
    else:
        logger.info(
            "no geonet moment tensor found for event: {}".format(
                event_id))
        raise AttributeError("geonet moment tensor for event {}"
                             "doesn't exist".format(event_id))


def grab_gcmt_moment_tensor(datetime, magnitude, path=None,
                            time_wiggle_sec=120, mag_wiggle=0.5):
    """
    Fetch global centroid moment tensor information from internal ndk files,
    if nothing is found then raise some errors. If multiple events found (e.g.
    temporally close foreshock and mainshock), allow user choice.

    Expects GCMT in .ndk files in directory structure following:
    /path/to/directories/YYYY/mmmyy.ndk
    e.g. path/to/directories/2009/jun09.ndk

    :type datetime: UTCDateTime or str
    :param datetime: event origin time
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
    if not isinstance(datetime, UTCDateTime):
        datetime = UTCDateTime(datetime)

    # Determine filename using datetime properties
    month = datetime.strftime('%b').lower()  # e.g. 'jul'
    year_short = datetime.strftime('%y')  # e.g. '19'
    year_long = datetime.strftime('%Y')  # e.g. '2019'

    fid = "{m}{y}.ndk".format(m=month, y=year_short)

    # If a path is given to the GCMT catalogs (hardcoded), search
    if path:
        fpath = os.path.join(path, year_long, fid)
        cat = read_events(fpath)
    # If no path, query GCMT directly using Obspy read_events
    # Try looking at the new files first
    else:
        # import HTTPerror to be able to catch it
        from urllib.error import HTTPError

        if not isinstance(datetime, UTCDateTime):
            datetime = UTCDateTime(datetime)

        logger.info("querying GCMT database for moment tensor")
        try:
            cat = read_events(
                "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
                "catalog/NEW_MONTHLY/{y}/{fid}".format(y=year_long, fid=fid)
            )
        except HTTPError:
            cat = read_events(
                "http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
                "catalog/NEW_QUICK/qcmt.ndk"
            )

    # GCMT catalogs contain all events for a span of time
    # filter catalogs using Obspy to find events with our specifications.
    # Magnitudes and origintimes are not always in agreement between agents
    # So allow fro some wiggle room
    cat_filt = cat.filter("time > {}".format(str(datetime - time_wiggle_sec)),
                          "time < {}".format(str(datetime + time_wiggle_sec)),
                          "magnitude >= {}".format(magnitude - mag_wiggle),
                          "magnitude <= {}".format(magnitude + mag_wiggle)
                          )
    # Filtering may remove all events from catalog, return multiple events, or
    # may return the event of choice
    if not len(cat_filt):
        logger.info(
            "no gcmt event found for {0} and M{1}".format(datetime, magnitude)
        )
        raise FileNotFoundError("No events found")
    elif len(cat_filt) > 1:
        logger.info(
            "multiple events found for {0} and M{1}".format(datetime, magnitude)
        )
        print("{} events found, choosing first".format(len(cat_filt)))
        return cat_filt[0]
    else:
        logger.info("gcmt event found matching criteria")
        return cat_filt[0]




