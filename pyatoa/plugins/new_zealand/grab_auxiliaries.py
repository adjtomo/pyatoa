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
        reader = csv.reader(open(fid))
        tag = "internal"
    else:
        # Request and open the CSV file. Assumed that GeoNet will keep their
        # moment-tensor information in their GitHub repository
        # Last accessed 23.6.19
        geonet_mt_csv = (
            "https://raw.githubusercontent.com/GeoNet/data/master/"
            "moment-tensor/GeoNet_CMT_solutions.csv"
        )
        response = requests.get(geonet_mt_csv)
        if not response.ok:
            warnings.warn("Github repo request failed", UserWarning)
            return None

        # Use CSV to parse through the returned response
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
            logger.info(f"geonet moment tensor {tag} for event: {event_id}")
            return moment_tensor
    else:
        logger.info(f"no geonet moment tensor found for event: {event_id}")
        raise AttributeError(f"geonet moment tensor for event {event_id}"
                             "doesn't exist")


def generate_focal_mechanism_from_geonet(mtlist, event=None):
    """
    For the New Zealand Tomography Problem

    Focal mechanisms created by John Ristau are written to a .csv file
    located on Github. This function will append information from the .csv file
    onto the Obspy event object so that all the information can be located in a
    single object

    :type mtlist: dict
    :param mtlist; row values from the GeoNet moment tensor csv file
    :type event: obspy.core.event.Event
    :param event: event to append focal mechanism to
    :rtype focal_mechanism: obspy.core.event.FocalMechanism
    :return focal_mechanism: generated focal mechanism
    """
    from obspy.core.event import source
    from obspy.core.event.base import Comment

    # Match the identifier with Goenet
    id_template = f"smi:local/geonetcsv/{mtlist['PublicID']}/{{}}"

    # Check that the input list is properly formatted
    if len(mtlist) != 32:
        print("geonet moment tensor list does not have the correct number"
              "of requisite components, should have 32")
        return

    # Generate the Nodal Plane objects containing strike-dip-rake
    nodal_plane_1 = source.NodalPlane(
        strike=mtlist['strike1'], dip=mtlist['dip1'], rake=mtlist['rake1']
    )
    nodal_plane_2 = source.NodalPlane(
        strike=mtlist['strike2'], dip=mtlist['dip2'], rake=mtlist['rake2']
    )
    nodal_planes = source.NodalPlanes(
        nodal_plane_1, nodal_plane_2, preferred_plane=1
    )

    # Create the Principal Axes as Axis objects
    tension_axis = source.Axis(
        azimuth=mtlist['Taz'], plunge=mtlist['Tpl'], length=mtlist['Tva']
    )
    null_axis = source.Axis(
        azimuth=mtlist['Naz'], plunge=mtlist['Npl'], length=mtlist['Nva']
    )
    pressure_axis = source.Axis(
        azimuth=mtlist['Paz'], plunge=mtlist['Ppl'], length=mtlist['Pva']
    )
    principal_axes = source.PrincipalAxes(
        t_axis=tension_axis, p_axis=pressure_axis, n_axis=null_axis
    )

    # Create the Moment Tensor object with correct units and scaling
    cv = 1E20 * 1E-7  # convert non-units, to dyne*cm, to N*m
    seismic_moment_in_nm = mtlist['Mo'] * 1E-7

    # Convert the XYZ coordinate system of GeoNet to an RTP coordinate system
    # expected in the CMTSOLUTION file of Specfem
    rtp = mt_transform(mt={"m_xx": mtlist['Mxx']*cv, "m_yy": mtlist['Myy']*cv,
                           "m_zz": mtlist['Mzz']*cv, "m_xy": mtlist['Mxy']*cv,
                           "m_xz": mtlist['Mxz']*cv, "m_yz": mtlist['Myz']*cv
                           },
                       method="xyz2rtp"
                       )
    tensor = source.Tensor(m_rr=rtp['m_rr'], m_tt=rtp['m_tt'],
                           m_pp=rtp['m_pp'], m_rt=rtp['m_rt'],
                           m_rp=rtp['m_rp'], m_tp=rtp['m_tp']
                           )
    # Create the source time function
    source_time_function = source.SourceTimeFunction(
        duration=2 * half_duration_from_m0(seismic_moment_in_nm)
    )

    # Generate a comment for provenance
    comment = Comment(force_resource_id=False,
                      text="Automatically generated by Pyatoa via GeoNet MT CSV"
                      )

    # Fill the moment tensor object
    moment_tensor = source.MomentTensor(
        force_resource_id=False, tensor=tensor,
        source_time_function=source_time_function,
        derived_origin_id=id_template.format('origin#ristau'),
        scalar_moment=seismic_moment_in_nm, double_couple=mtlist['DC']/100,
        variance_reduction=mtlist['VR'], comment=comment
        )

    # Finally, assemble the Focal Mechanism. Force a resource id so that
    # the event can identify its preferred focal mechanism
    focal_mechanism = source.FocalMechanism(
        force_resource_id=True, nodal_planes=nodal_planes,
        moment_tensor=moment_tensor, principal_axes=principal_axes,
        comments=[comment]
        )

    # Append the focal mechanisms to the event object. Set the preferred
    # focal mechanism so that this attribute can be used in the future
    if event:
        event.focal_mechanisms = [focal_mechanism]
        event.preferred_focal_mechanism_id = focal_mechanism.resource_id
        return event, focal_mechanism
    # If no event is given, just return the focal mechanism
    else:
        return None, focal_mechanism


