"""
Gather auxiliary data specifically relevant for New Zealand seismology.
"""
import csv
import requests
from obspy import UTCDateTime
from obspy.core.event import source
from obspy.core.event.base import Comment
from pyatoa import logger
from pyatoa.utils.srcrcv import mt_transform, half_duration_from_m0


def get_geonet_mt(event_id, csv_fid=None):
    """
    Get moment tensor information from a internal csv file,
    or from an external github repository query.
    Only relevant to the new zealand tomography problem.
    Geonet moment tensors stored with a specific column format.

    :type event_id: str
    :param event_id: unique event identifier
    :type csv_fid: str
    :param csv_fid: optional path to GeoNet CMT solution file that is stored 
        locally on disk, will be accessed before querying web service
    :rtype moment_tensor: dict
    :return moment_tensor: dictionary created from rows of csv file
    """
    reader = None
    if csv_fid is not None:
        try:
            reader = csv.reader(open(csv_fid, 'r'), delimiter=',')
        except FileNotFoundError:
            pass

    if reader is None:
        # Request and open the CSV file. Assumed that GeoNet will keep their
        # moment-tensor information in their GitHub repository
        # Last accessed 23.6.19
        geonet_mt_csv = (
            "https://raw.githubusercontent.com/GeoNet/data/master/"
            "moment-tensor/GeoNet_CMT_solutions.csv"
        )
        response = requests.get(geonet_mt_csv)
        if not response.ok:
            raise FileNotFoundError(f"Response from {geonet_mt_csv} not ok")

        reader = csv.reader(response.text.splitlines(), delimiter=',')

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
            logger.info(f"geonet moment tensor found for: {event_id}")
            return moment_tensor
    else:
        raise AttributeError(f"no geonet moment tensor found for: {event_id}")


def geonet_mt(event_id, units, event=None, csv_fid=None):
    """
    Focal mechanisms created by John Ristau are written to a .csv file
    located on Github. This function will append information from the .csv file
    onto the Obspy event object so that all the information can be located in a
    single object

    :type event_id: str
    :param event_id: unique event identifier
    :type units: str
    :param units: output units of the focal mechanism, either: 
        'dynecm': for dyne*cm  or 
        'nm': for Newton*meter
    :type event: obspy.core.event.Event
    :param event: event to append focal mechanism to
    :rtype focal_mechanism: obspy.core.event.FocalMechanism
    :return focal_mechanism: generated focal mechanism
    """
    assert(units in ["dynecm", "nm"]), "units must be 'dynecm' or 'nm'"

    mtlist = get_geonet_mt(event_id, csv_fid=csv_fid)

    # Match the identifier with Goenet
    id_template = f"smi:local/geonetcsv/{mtlist['PublicID']}/{{}}"

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
    if units == "nm":
        c = 1E-7  # conversion from dyne*cm to N*m
        logger.debug(f"GeoNet moment tensor is in units of Newton*meters")
    elif units == "dynecm":
        c = 1
        logger.debug(f"GeoNet moment tensor is in units of dyne*cm")

    # CV is the conversion from non-units to the desired output units
    cv = 1E20 * c
    seismic_moment_in_nm = mtlist['Mo'] * c

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
    comment = Comment(force_resource_id=True,
                      text="Automatically generated by Pyatoa via GeoNet MT CSV"
                      )

    # Fill the moment tensor object
    moment_tensor = source.MomentTensor(
        force_resource_id=True, tensor=tensor,
        source_time_function=source_time_function,
        # !!!
        # This doesn't play nice with obspy.Catalog.write(format='CMTSOLUTION')
        # so ignore the origin id
        # derived_origin_id=id_template.format('origin#ristau'),
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


