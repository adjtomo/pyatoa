"""
Pyatoa relies on data structure being ordered and consistent throughout all the
various bits of data required. functions here will aid in reshaping data
into the correct formats
"""
import os
import numpy as np


def distribute_dataless(path_to_response, inv):
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

    :type window: pyflex.Window
    :param window: misfit window calcualted by pyflex
    :rtype win_dict: dict
    :return win_dict: formatted dictionary of values for pyasdf auxiliary data
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


def write_adj_src_to_ascii(ds, model_number, comp_list=["N", "E", "Z"]):
    """
    take AdjointSource auxiliary data from a pyasdf dataset and write out
    the adjoint sources into ascii files with proper formatting, for input
    into PyASDF
    :param ds:
    :param model_number:
    :param comp_list:
    :return:
    """
    def write_to_ascii(f, array):
        """
        function used to write the ascii in the correct format
        :param array:
        :return:
        """
        for dt, amp in array:
            if dt == 0. and not amp == 0.:
                dt = 0
                adj_formatter = "{dt:>14}{amp:18.6f}\n"
            if not dt == 0 and amp == 0.:
                amp = 0
                adj_formatter = "{dt:14.6f}{amp:>18}\n"
            else:
                adj_formatter = "{dt:14.6f}{amp:18.6f}\n"
            f.write(adj_formatter.format(dt=dt, amp=amp))

    shortcut = ds.auxiliary_data.AdjointSources[model_number]
    for adj_src in shortcut.list():
        station = shortcut[adj_src].path.replace('_', '.')
        fid = "./{}.adj".format(station)
        with open(fid, "w") as f:
            write_to_ascii(f, shortcut[adj_src].data.value)
        if comp_list:
            for comp in comp_list:
                station_check = station[:-1] + comp
                if station_check.replace('.', '_') not in shortcut.list():
                    blank_adj_src = shortcut[adj_src].data.value
                    blank_adj_src[:, 1] = np.zeros(len(blank_adj_src[:, 1]))
                    blank_template = "./{}.adj".format(station_check)
                    with open(blank_template, "w") as b:
                        write_to_ascii(b, blank_adj_src)


def write_adj_src_to_asdf(adj_src, ds, tag, time_offset):
    """
    NOTE: Stolen and modified from Pyadjoint source code:
          pyadjoint.adjoint_source.write_to_asdf()

    Writes the adjoint source to an ASDF file.
    Note: For now it is assumed SPECFEM will be using the adjoint source

    :type adj_src: pyadjoint.asdf_data_set.ASDFDataSet
    :param adj_src: adjoint source to save
    :type ds: pyasdf.asdf_data_set.ASDFDataSet
    :type tag: str
    :param tag: internal pathing for save location in the auxiliary data attr.
    :param ds: The ASDF data structure read in using pyasdf.
    :type time_offset: float
    :param time_offset: The temporal offset of the first sample in seconds.
        This is required if using the adjoint source as input to SPECFEM.
    .. rubric:: SPECFEM
    SPECFEM requires one additional parameter: the temporal offset of the
    first sample in seconds. The following example sets the time of the
    first sample in the adjoint source to ``-10``.
    >>> adj_src.write_to_asdf(ds, time_offset=-10,
    ...               coordinates={'latitude':19.2,
    ...                            'longitude':13.4,
    ...                            'elevation_in_m':2.0})
    """
    # Convert the adjoint source to SPECFEM format
    l = len(adj_src.adjoint_source)
    specfem_adj_source = np.empty((l, 2))
    specfem_adj_source[:, 0] = np.linspace(0, (l - 1) * adj_src.dt, l)
    specfem_adj_source[:, 0] += time_offset
    specfem_adj_source[:, 1] = adj_src.adjoint_source[::-1]

    station_id = "{net}.{sta}".format(net=adj_src.network, sta=adj_src.station)
    coordinates = ds.waveforms["{net}.{sta}".format(
        net=adj_src.network, sta=adj_src.station)].coordinates

    # Safeguard against funny types in the coordinates dictionary
    latitude = float(coordinates["latitude"])
    longitude = float(coordinates["longitude"])
    elevation_in_m = float(coordinates["elevation_in_m"])

    parameters = {"dt": adj_src.dt, "misfit_value": adj_src.misfit,
                  "adjoint_source_type": adj_src.adj_src_type,
                  "min_period": adj_src.min_period,
                  "max_period": adj_src.max_period,
                  "latitude": latitude, "longitude": longitude,
                  "elevation_in_m": elevation_in_m,
                  "station_id": station_id, "component": adj_src.component,
                  "units": "m"}

    ds.add_auxiliary_data(data=specfem_adj_source, data_type="AdjointSources",
                          path=tag, parameters=parameters)