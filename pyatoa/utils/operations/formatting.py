"""
pyatoa relies on data structure being ordered and consistent throughout all the
various bits of data required. functions here will aid in reshaping data
into the correct formats
"""
import os


def distribute_dataless(path_to_response,inventory):
    """
    Response files written through obspy come out as a single object, but pyatoa
    will look for response information from individual components and individual
    stations. Distrubute this dataless information into the necessary components
    """
    inner_folder = '{STA}.{NET}'
    fid_template = 'RESP.{NET}.{STA}.{LOC}.{CHA}'
    full_template = os.path.join(path_to_response,inner_folder,fid_template)
    for net in inv:
        for sta in net:
            try:
                os.mkdir(os.path.join(path_to_response,inner_folder.format(
                    STA=sta.code, NET=net.code))
                    )
            except FileExistsError:
                pass
            for cha in sta:
                inv_temp = inv.select(network=net.code, station=sta.code,
                                      location=cha.location_code,
                                      channel=cha.code)
                inv_temp.write(full_template.format(
                    STA=sta.code, NET=net.code, LOC=cha.location_code,
                    CHA=cha.code),
                    format='STATIONXML'
                    )


def create_window_dictionary(window):
    """
    HDF5 doesnt play nice with nonstandard objects in dictionaries, e.g.
    nested dictionaries, UTCDateTime objects. So remake the pyflex window
    json dictionary into something that will sit well in a pyasdf object
    """
    win_dict = window._get_json_content()

    # change UTCDateTime objects into strings
    win_dict['absolute_endtime'] = str(win_dict['absolute_endtime'])
    win_dict['absolute_starttime'] = str(win_dict['absolute_starttime'])
    win_dict['time_of_first_sample'] = str(win_dict['time_of_first_sample'])

    phase_arrivals = win_dict['phase_arrivals']
    for phase in phase_arrivals:
        win_dict['phase_arrival_{}'.format(phase['name'])] = phase['time']
    win_dict.pop('phase_arrivals')

    return win_dict


def moment_tensor_list_to_objects(mtlist):
    """
    UNFINISHED
    event objects fetched by obspy do not natively come with any moment tensor
    or nodal plane information, that is stored separately in a .csv file
    located on github. For the tomography problem we need this information, so
    this functino will append information from the .csv file onto the obspy
    event object so that all the information can be located in a single object
    :param event:
    :param geonet_moment_tensor_list:
    :return:
    """
    from obspy.core.event import ResourceIdentifier
    import obspy.core.event.source as eventcore

    id_template = "smi:local/geonetcsv/{0}/{1}".format(mtlist['PublicID'],'{}')
    if len(mtlist) != 32:
        print("geonet moment tensor list does not have the correct number"
              "of requisite components, should have 32")
        return
    nodal_plane_1 = eventcore.NodalPlane(strike=mtlist['strike1'],
                                         dip=mtlist['dip1'],
                                         rake=mtlist['rake1']
                                         )
    nodal_plane_2 = eventcore.NodalPlane(strike=mtlist['strike2'],
                                         dip=mtlist['dip2'],
                                         rake=mtlist['rake2']
                                         )
    nodal_planes = eventcore.NodalPlanes(nodal_plane_1, nodal_plane_2,
                                         preferred_plane=1)
    moment_tensor = eventcore.MomentTensor(
        resource_id=id_template.format('momenttensor'),
        derived_origin_id=id_template.format('origin#ristau'),
        scalar_moment=mtlist['Mo']*1E-7
        )

    focal_mechanism = eventcore.FocalMechanism(
        resource_id=ResourceIdentifier(
            "smi:local/geonetcsv/{}/focal_mechanism".format(mtlist['PublicID'])
            ),
        nodal_planes=nodal_planes, moment_tensor=moment_tensor
        )

    def write_to_asdf(adjoint_source, ds, time_offset, coordinates=None):
        """
        !!! Stolen and modified from Pyadjoint source code

        Writes the adjoint source to an ASDF file.
        Note: For now it is assumed SPECFEM will be using the adjoint source
        :param ds: The ASDF data structure read in using pyasdf.
        :type ds: str
        :param time_offset: The temporal offset of the first sample in seconds.
            This is required if using the adjoint source as input to SPECFEM.
        :type time_offset: float
        :param coordinates: If given, the coordinates of the adjoint source.
            The 'latitude', 'longitude', and 'elevation_in_m' of the adjoint
            source must be defined.
        :type coordinates: list
        .. rubric:: SPECFEM
        SPECFEM requires one additional parameter: the temporal offset of the
        first sample in seconds. The following example sets the time of the
        first sample in the adjoint source to ``-10``.
        >>> adj_src.write_to_asdf(ds, time_offset=-10,
        ...               coordinates={'latitude':19.2,
        ...                            'longitude':13.4,
        ...                            'elevation_in_m':2.0})
        """
        # Import here to not have a global dependency on pyasdf
        from pyasdf.exceptions import NoStationXMLForStation

        # Convert the adjoint source to SPECFEM format
        l = len(self.adjoint_source)
        specfem_adj_source = np.empty((l, 2))
        specfem_adj_source[:, 0] = np.linspace(0, (l - 1) * self.dt, l)
        specfem_adj_source[:, 1] = time_offset
        specfem_adj_source[:, 1] = self.adjoint_source[::-1]

        tag = "%s_%s_%s" % (self.network, self.station, self.component)
        min_period = self.min_period
        max_period = self.max_period
        component = self.component
        station_id = "%s.%s" % (self.network, self.station)

        if coordinates:
            # If given, all three coordinates must be present
            if {"latitude", "longitude", "elevation_in_m"}.difference(
                    set(coordinates.keys())):
                raise ValueError(
                    "'latitude', 'longitude', and 'elevation_in_m'"
                    " must be given")
        else:
            try:
                coordinates = ds.waveforms[
                    "%s.%s" % (self.network, self.station)].coordinates
            except NoStationXMLForStation:
                raise ValueError("Coordinates must either be given "
                                 "directly or already be part of the "
                                 "ASDF file")

        # Safeguard against funny types in the coordinates dictionary
        latitude = float(coordinates["latitude"])
        longitude = float(coordinates["longitude"])
        elevation_in_m = float(coordinates["elevation_in_m"])

        parameters = {"dt": self.dt, "misfit_value": self.misfit,
                      "adjoint_source_type": self.adj_src_type,
                      "min_period": min_period, "max_period": max_period,
                      "latitude": latitude, "longitude": longitude,
                      "elevation_in_m": elevation_in_m,
                      "station_id": station_id, "component": component,
                      "units": "m"}

        # Use pyasdf to add auxiliary data to the ASDF file
        ds.add_auxiliary_data(data=specfem_adj_source,
                              data_type="AdjointSource", path=tag,
                                parameters=parameters)