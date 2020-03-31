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
from pyatoa.utils.process import is_preprocessed, change_syn_units


def scale_beacon_amplitudes(st):
    """
    Return a waveform with scaled amplitudes based on station identifiers.
    This scaling is empirical and is required to boost amplitudes of our
    CMG-40T60s instruments which for some reason have order of magnitude lower
    amplitudes. The 30s instruments are okay but amplitudes are larger than
    GeoNet stations (maybe site response?), so we put in a small scaling factor.

    :type st: obspy.core.stream.Stream
    :param st: stream containing waveform data to be scaled
    :rtype: obspy.core.stream.Stream
    :return: stream with scaled waveforms
    """
    st_scale = st.copy()

    # Amplitude scaling
    # scale_1 = 22.
    # scale_2 = 0.4
    # scale_3 = 57.5
    scale_1 = 100.
    scale_2 = 2.
    scale_3 = 250.
    scales = [scale_1, scale_2, scale_3]

    # Stations that belong to each scaling group
    group_1 = [1, 2, 6, 7, 8, 9, 13, 14, 15]  # 60s instruments
    group_2 = [3, 10, 11, 12, 16, 17, 19, 20, 21]  # 30s instruments
    group_3 = [4, 5, 18]  # somehow worse 60s instruments
    groups = [group_1, group_2, group_3]

    # Determine group by checking station name
    try: 
        idx = int(st_scale[0].stats.station[2:])
    except ValueError:
        logger.debug("station code does not match BEACON station formatting")
        return st 
    for s, g in zip(scales, groups):
        if idx in g:
            logger.debug(f"scaling {st_scale[0].get_id()} by {s}")
            for tr in st_scale:
                tr.data *= s
            return st_scale


def preproc(mgmt, choice, water_level=60, corners=4, taper_percentage=0.05):
    """
    Preprocess waveform data with differences in preprocessing for
    certain temporary network data

    :type mgmt: pyatoa.core.manager.Manager
    :param mgmt: Manager class that should contain a Config object as well as
        waveform data and inventory
    :type choice: str
    :param choice: option to preprocess observed, synthetic or both
        available: 'obs', 'syn'
    :type water_level: int
    :param water_level: water level for response removal
    :type corners: int
    :param corners: value of the filter corners, i.e. steepness of filter edge
    :type taper_percentage: float
    :param taper_percentage: amount to taper ends of waveform
    """
    # Copy the stream to avoid editing in place
    if choice == "syn":
        st = mgmt.st_syn.copy()
    elif choice == "obs":
        st = mgmt.st_obs.copy()
    if is_preprocessed(st):
        return st

    # Standard preprocessing before specific preprocessing
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=taper_percentage)

    # Observed specific data preprocessing includes response and rotating to ZNE
    # Bannister data does not require instrument response removal
    if choice == "obs" and not mgmt.config.synthetics_only and \
            st[0].stats.network not in ["Z8", "ZX"]:
        # Occasionally, inventory issues arise, as ValueErrors due to
        # station availability, e.g. NZ.COVZ. Try/except to catch these.
        try:
            st.attach_response(mgmt.inv)
            st.remove_response(output=mgmt.config.unit_output,
                               water_level=water_level,
                               plot=False)
        except ValueError:
            logger.debug(f"Error removing response from {st[0].get_id()}")
            return None
        logger.debug("remove response, units of {}".format(
            mgmt.config.unit_output)
        )

        # Clean up streams after response removal
        st.detrend("linear")
        st.detrend("demean")
        st.taper(max_percentage=taper_percentage)

        # Rotate streams if they are not in the ZNE coordinate system, e.g. Z12
        st.rotate(method="->ZNE", inventory=mgmt.inv)
    # Synthetic specific data processing includes changing units
    else:
        if mgmt.config.unit_output != mgmt.config.synthetic_unit:
            logger.debug("unit output and synthetic output do not match, "
                         "adjusting")
            st = change_syn_units(st, current=mgmt.config.unit_output,
                                  desired=mgmt.config.synthetic_unit)
            st.detrend("linear")
            st.detrend("demean")
            st.taper(max_percentage=taper_percentage)

    # Rotate the given stream from standard NEZ to RTZ
    if mgmt.baz:
        st.rotate(method="NE->RT", back_azimuth=baz)
        logger.debug(f"rotating NE->RT by {baz} degrees")

    # Filter data using ObsPy Butterworth filters. Zerophase avoids phase shift
    # Bandpass filter
    if mgmt.config.min_period and mgmt.config.max_period:
        st.filter("bandpass",
                  freqmin=1/mgmt.config.max_period,
                  freqmax=1/mgmt.config.min_period, corners=corners,
                  zerophase=True
                  )
        logger.debug(
            f"bandpass {mgmt.config.min_period}-{mgmt.config.max_period}s")
    # Highpass if only minimum period given
    elif mgmt.config.min_period:
        st.filter("highpass", freq=mgmt.config.min_period, corners=corners,
                  zerophase=True)
        logger.debug(f"highpass {mgmt.config.min_period}s")
    # Highpass if only minimum period given
    elif mgmt.config.max_period:
        st.filter("lowpass", freq=mgmt.config.max_period, corners=corners,
                  zerophase=True)
        logger.debug(f"lowpass {mgmt.config.max_period}s")

    # Beacon data requires amplitude scaling, perform after filtering to
    # avoid boosting the amplitude of noise
    if st[0].stats.network == "XX" and not mgmt.config.synthetics_only:
        st = scale_beacon_amplitudes(st)

    return st
