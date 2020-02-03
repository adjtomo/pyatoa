"""
Custom preprocessing tools to make accomodations for New Zealand temporary 
broadband stations which may require amplitude scaling or non-standard 
preprocessing techniques

Tools for processing obspy.Stream or obspy.Trace objects
Used for preprocessing data through filtering and tapering, zero padding etc.
Also contains tools for synthetic traces such as source time function
convolutions
"""
from pyatoa import logger
from pyatoa.utils.tools.process import is_preprocessed, change_syn_units


def scale_beacon_amplitudes(st):
    """
    Return a waveform with scaled amplitudes based on station identifiers

    :type st: obspy.core.stream.Stream
    :param st: stream containing waveform data to be scaled
    :rtype: obspy.core.stream.Stream
    :return: stream with scaled waveforms
    """
    st_scale = st.copy()

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

    # Determine group by checking station name
    try: 
        idx = int(st_scale[0].stats.station[:2])
    except ValueError:
        logger.debug("Station code does not match BEACON station formatting")
        return  
    for s, g in zip(scales, groups):
        if idx in g:
            logger.debug(f"Scaling {st_scale[0].get_id()} by {s}")
            for tr in st_scale:
                tr.data *= s
            return st_scale


def preproc(st_original, inv=None, unit_output="VEL", synthetic_unit=None,
            back_azimuth=None, filter_bounds=(10, 30), water_level=60,
            corners=4, taper_percentage=0.05):
    """
    Preprocess waveform data. Almost identical to the default preprocessing
    script except there are some checks on temporary data to differ processing
    slightly

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
    if is_preprocessed(st):
        return st

    # Standard preprocessing
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=taper_percentage)

    # If inventory is given, assume working with observation data
    if inv:
        # Bannister data does not require instrument response removal
        if st[0].stats.network in ["Z8", "ZX"]:
            # Clean up streams after response removal
            st.detrend("linear")
            st.detrend("demean")
            st.taper(max_percentage=taper_percentage)

            # Rotate streams if they are not in the ZNE coordinate system
            st.rotate(method="->ZNE", inventory=inv)
        else:
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

            # Beacon data requires amplitude scaling
            if st[0].stats.network == "XX":
                st = scale_beacon_amplitudes(st)

    # No inventory means synthetic data
    elif not inv:
        if unit_output != synthetic_unit:
            logger.debug("unit output and synthetic output do not match, "
                         "adjusting")
            st = change_syn_units(st, current=unit_output,
                                  desired=synthetic_unit)
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
