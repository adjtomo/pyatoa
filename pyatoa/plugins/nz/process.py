"""
Custom preprocessing tools to make accomodations for New Zealand temporary 
broadband stations which may require amplitude scaling or non-standard 
preprocessing techniques

Tools for processing obspy.Stream or obspy.Trace objects
Used for preprocessing data through filtering and tapering, zero padding etc.
Also contains tools for synthetic traces such as source time function
convolutions
"""
import warnings
import numpy as np

from pyatoa import logger


def _is_preprocessed(st):
    """
    Small check to make sure a stream object has not yet been run through
    preprocessing. Simple, as it assumes that a fresh stream will have no
    processing attribute in their stats, or if they do, will not have been
    filtered (getting cut waveforms from FDSN appends a 'trim' stat).
    :type st: obspy.stream.Stream
    :param st: stream to check processing on
    :rtype: bool
    :return: if preprocessing has occurred
    """
    for tr in st:
        if hasattr(tr.stats, 'processing'):
            for processing in tr.stats.processing:
                # A little hacky, but processing flag will have the str
                # ..': filter(options'... to signify that a filter is applied
                if 'filter(' in processing:
                    warnings.warn("stream already preprocessed", UserWarning)
                    return True

    # If nothing found, return False
    return False
    

def scale_beacon_amplitudes(st):
    """
    Return a waveform with scaled amplitudes based on station identifiers

    :type st: obspy.core.stream.Stream
    :param st: stream containing waveform data to be scaled
    """
    # Amplitude scaling
    scale_1 = 22.
    scale_2 = 0.4
    scale_3 = 57.5
    scales = [scale_1, scale_2, scale_3]

    # Stations that belong to each scaling group
    group_1 = [1, 2, 6, 7, 8, 9, 13, 14, 15]
    group_2 = [3, 10, 11, 12, 16, 17, 19, 20, 21]
    group_3 = [4, 5, 18]
    groups = [group_1, group_2, group_3]

    try: 
        idx = int(st.stats.station[:2])
    except ValueError:
        print("Station code does not match BEACON station formatting")
        return  
    for s, g in zip(scales, groups):
        if idx in g:
            for tr in st:
                tr.data *= s
            return st


def preproc(st_original, inv=None, unit_output="VEL", synthetic_unit=None,
            back_azimuth=None, filter_bounds=(10, 30), water_level=60,
            corners=4, taper_percentage=0.05):
    """
    Preprocess waveform data. Similar to the default preprocessing script
    except there are some checks on temporary data to differ processing slightly

    :type st_original: obspy.stream.Stream
    :param st_original: stream object to process
    :type inv: obspy.core.inventory.Inventory or None
    :param inv: inventory containing relevant network and stations
    :type unit_output: str
    :param unit_output: output of response removal, available:
        'DISP', 'VEL', 'ACC'
    :type synthetic_unit: str
    :param synthetic_unit: units of synthetic traces, same available as unit
    :type back_azimuth: float
    :param back_azimuth: back azimuth in degrees
    :type filter_bounds: list or tuple of float
    :param filter_bounds: (min period, max_period), units of Seconds
    :type water_level: int
    :param water_level: water level for response removal
    :type corners: int
    :param corners: value of the filter corners, i.e. steepness of filter edge
    :type taper_percentage: float
    :param taper_percentage: amount to taper ends of waveform
    :rtype st: obspy.stream.Stream
    :return st: preprocessed stream object
    """
    # Copy the stream to avoid editing in place
    st = st_original.copy()

    warnings.filterwarnings("ignore", category=FutureWarning)
    if _is_preprocessed(st):
        return st

    # Standard preprocessing
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=taper_percentage)

    # If inventory is given, assume working with observation data
    if inv:
        # Occasionally, inventory issues arise, as ValueErrors due to 
        # station availability, e.g. NZ.COVZ. Try/except to catch these.
        try:
            st.attach_response(inv)
            st.remove_response(output=unit_output,
                               water_level=water_level,
                               plot=False)
        except ValueError:
            logger.debug(f"Error removing response from {st[0].get_id()}")
            return None
    
        logger.debug("remove response, units of {}".format(unit_output))

        # Clean up streams after response removal
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=taper_percentage)

        # Rotate streams if they are not in the ZNE coordinate system
        st.rotate(method="->ZNE", inventory=inv)

    # No inventory means synthetic data
    elif not inv:
        if unit_output != synthetic_unit:
            logger.debug(
                "unit output and synthetic output do not match, adjusting")
            st.detrend("linear")
            st.detrend("demean")

        st.taper(max_percentage=taper_percentage)
    
    # Rotate the given stream from standard NEZ to RTZ
    if back_azimuth:
        st.rotate(method="NE->RT", back_azimuth=back_azimuth)
        logger.debug(f"rotating NE->RT by {back_azimuth} degrees")

    # Filter data using ObsPy Butterworth filters. Zerophase avoids phase shift
    if filter_bounds is not None:
        st.filter('bandpass', freqmin=1/filter_bounds[1],
                  freqmax=1/filter_bounds[0], corners=corners, zerophase=True
                  )
        logger.debug(f"filter {filter_bounds[0]}s to {filter_bounds[1]}s")

    return st

def _synthetic_units()

        # Determine the difference between synthetic unit and observed unit
        diff_dict = {"DISP": 1, "VEL": 2, "ACC": 3}
        difference = diff_dict[unit_output] - diff_dict[synthetic_unit]

        # Integrate or differentiate stream to retrieve correct units
        if difference == 1:
            logger.debug("integrating synthetic data")
            st.integrate(method="cumtrapz")
        elif difference == 2:
            logger.debug("double integrating synthetic data")
            st.integrate(method="cumtrapz").integrate(method="cumtrapz")
        elif difference == -1:
            logger.debug("differentiating synthetic data")
            st.differentiate(method="gradient")
        elif difference == -2:
            logger.debug("double differentiating synthetic data")
            st.differentiate(
                method="gradient").differentiate(method="gradient")


