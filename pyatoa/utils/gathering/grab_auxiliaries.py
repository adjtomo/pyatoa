"""
Various scripts to grab auxiliary data such as moment tensor information,
station information, and fault information. these are all specific to the
New Zealand tomography problem, and therefore paths are hard coded
"""
import os
import csv

from obspy import UTCDateTime, read_events
from pyatoa import logger


def hardcode_paths():
    """
    TO DO: try to remove this hardcoding, or avoid it if these aren't found
    personal development convenience function to hardcode in a path dictionary
    for fetching internally stored information

    :rtype paths: dict
    :return: dictionary containing hardcoded pathnames
    """
    split = os.getcwd().split('/')
    basecheck = os.path.join(split[0], split[1], split[2])
    import pdb;pdb.set_trace()
    if basecheck == "seis/prj":
        where = "GNS"
        datafolder = "/seis/prj/fwi/bchow/data"
    elif basecheck == "Users/chowbr":
        where = "VIC"
        datafolder = "/Users/chowbr/Documents/subduction/data"
    else:
        where = "MAUI"
    
    # Set hardcoded paths based on system
    if where != "MAUI":
        paths = {"faults": os.path.join(datafolder, "FAULTS", ''),
                 "stations": os.path.join(datafolder, "STATIONXML", "MASTER",
                                          "master_inventory.xml"),
                 "geonet_mt": os.path.join(datafolder, "GEONET", "data",
                                           "moment-tensor",
                                           "GeoNet_CMT_solutions.csv"),
                 "gcmt_mt": os.path.join(datafolder, "GCMT")
                 }
    else:
        datafolder = ("/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/"
                      "primer/auxiliary_data")

        paths = {"faults": os.path.join(datafolder, "faults", ''),
                 "stations": os.path.join(datafolder, "stationxml",
                                          "master_inventory.xml"),
                 "geonet_mt": os.path.join(datafolder, "geonet", "data",
                                           "moment-tensor",
                                           "GeoNet_CMT_solutions.csv"),
                 "gcmt_mt": os.path.join(datafolder, "gcmt")
                 }
    return paths


def grab_geonet_moment_tensor(event_id):
    """
    fetch moment tensor information from an internal csv file, only relevant
    to the new zealand tomography problem. geonet moment tensors stored with
    a specific column format, csv file can be found here:

    https://github.com/GeoNet/data/tree/master/moment-tensor

    :type event_id: str
    :param event_id: unique event identifier
    :rtype moment_tensor: dict
    :return moment_tensor: dictionary created from rows of csv file
    """
    with open(hardcode_paths()['geonet_mt']) as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if i == 0:
                tags = row
            if row[0] == event_id:
                values = []
                for t, v in zip(tags, row):
                    if t == "Date":
                        values.append(UTCDateTime(v))
                    elif t == "PublicID":
                        values.append(v)
                    else:
                        values.append(float(v))

                moment_tensor = dict(zip(tags, values))
                logger.info("geonet moment tensor found for event: {}".format(
                    event_id))
                return moment_tensor
        else:
            logger.info(
                "no geonet moment tensor found for event: {}".format(
                    event_id))
            raise AttributeError("geonet moment tensor for event {}"
                                 "doesn't exist".format(event_id))


def grab_gcmt_moment_tensor(datetime, magnitude):
    """
    Fetch global centroid moment tensor information from internal ndk files,
    if nothing is found then raise some errors. If multiple events found (e.g.
    temporally close foreshock and mainshock), allow user choice.

    :type datetime: UTCDateTime or str
    :param datetime: event origin time
    :type magnitude: float
    :param magnitude: centroid moment magnitude for event lookup
    :rtype event: obspy.core.event.Event
    :return event: event object for given earthquake
    """
    if not isinstance(datetime, UTCDateTime):
        datetime = UTCDateTime(datetime)

    year = str(datetime.year)
    month = {1: "jan", 2: "feb", 3: "mar", 4: "apr", 5: "may",
             6: "jun", 7: "jul", 8: "aug", 9: "sep", 10: "oct",
             11: "nov", 12: "dec"}[datetime.month]

    fid = "{m}{y}.ndk".format(m=month, y=year[2:])
    fpath = os.path.join(hardcode_paths()['gcmt_mt'], year, fid)
    try:
        cat = read_events(fpath)
    except FileNotFoundError:
        logger.info("no gcmt ndk file found internal, searching external")
        # import HTTPerror to be able to catch it
        from urllib.error import HTTPError
        try:
            cat = read_events(
                "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
                "catalog/NEW_MONTHLY/{y}/{fid}".format(y=year, fid=fid)
            )
        except HTTPError:
            cat = read_events(
                "http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
                "catalog/NEW_QUICK/qcmt.ndk"
            )

    cat_filt = cat.filter("time > {}".format(str(datetime - 60)),
                          "time < {}".format(str(datetime + 60)),
                          "magnitude >= {}".format(magnitude - .5),
                          "magnitude <= {}".format(magnitude + .5)
                          )
    if not len(cat_filt):
        logger.info(
            "no gcmt event found for datetime {0} and magnitude {1}".format(
                datetime, magnitude)
        )
        raise FileNotFoundError("No events found")
    elif len(cat_filt) > 1:
        logger.info(
            "multiple events found for datetime {0} and magnitude{1}".format(
                datetime,magnitude)
        )
        print("{} events found, choose from list:".format(len(cat_filt)))
        print("{0} {1}\n{2}".format(datetime, magnitude, cat_filt))
        choice = int(input("Event number (index from 0): "))
        return cat_filt[choice]
    else:
        logger.info("gcmt event found matching criteria")
        return cat_filt[0]


def timeshift_halfduration(gcmt_event, geonet_list):
    """
    Deprecated, and incorrect, will be deleted.
    calcuate the absolute time shift between centroid time and hypocenter time
    to shift the synthetic seismogram into absolute time using the equation:

    t_abs = t_pde + time shift + t_syn

    :type gcmt_event: obspy.core.event.Event
    :param gcmt_event:
    """
    import warnings
    warnings.warn("Incorrect function", DeprecationWarning)
    hypocenter_time = gcmt_event.origins[0].time
    centroid_time = [i.time for i in gcmt_event.origins
                     if i.origin_type == "centroid"][0]
    time_shift = abs(hypocenter_time - centroid_time)
    moment_tensor = gcmt_event.focal_mechanisms[0].moment_tensor
    half_duration = (moment_tensor.source_time_function['duration'])/2
    return time_shift, half_duration

