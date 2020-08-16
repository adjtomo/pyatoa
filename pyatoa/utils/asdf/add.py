"""
ASDF Datasets can be given auxiliary data to supplement the existing waveform,
event and station information contained. The functions contained in this script
add new auxiliary data structures to existing ASDF datasets
"""
import warnings
import numpy as np
from pyatoa import logger
from pyatoa.utils.form import channel_code


def add_misfit_windows(windows, ds, path):
    """
    Write Pyflex misfit windows into the auxiliary data of an ASDFDataSet

    :type windows: dict of list of pyflex.Window
    :param windows: dictionary of lists of window objects with keys
        corresponding to components related to each window
    :type ds: pyasdf.ASDFDataSet
    :param ds: ASDF data set to save windows to
    :type path: str
    :param path: internal pathing to save location of auxiliary data
    """
    # Save windows by component
    for comp in windows.keys():
        for i, win in enumerate(windows[comp]):
            # Figure out how to tag the data in the dataset
            net, sta, loc, cha = win.channel_id.split(".")

            # net_sta_comp_i: e.g. NZ_BFZ_Z_0
            window_tag = "_".join([net, sta, cha[-1], str(i)])

            # Turn UTCDateTime objects into strings
            wdict = win._get_json_content()
            wdict["absolute_endtime"] = str(wdict["absolute_endtime"])
            wdict["absolute_starttime"] = str(wdict["absolute_starttime"])
            wdict["time_of_first_sample"] = str(wdict["time_of_first_sample"])

            # Flatten nested dictionary
            phase_arrivals = wdict["phase_arrivals"]
            for phase in phase_arrivals:
                wdict[f"phase_arrival_{phase['name']}"] = phase["time"]
            wdict.pop("phase_arrivals")

            # Write windows into dataset
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ds.add_auxiliary_data(data=np.array([win.left, win.right]),
                                      data_type="MisfitWindows",
                                      parameters=wdict,
                                      path=f"{path}/{window_tag}"
                                      )


def add_adjoint_sources(adjsrcs, ds, path, time_offset):
    """
    NOTE: Borrowed and modified from Pyadjoint source code:
          pyadjoint.adjoint_source.write_to_asdf()

    Writes the adjoint source to an ASDF file.
    Note: For now it is assumed SPECFEM will be using the adjoint source

    :type adjsrcs: list of pyadjoint.asdf_data_set.ASDFDataSet
    :param adjsrcs: adjoint source to save
    :type ds: pyasdf.asdf_data_set.ASDFDataSet
    :type path: str
    :param path: internal pathing for save location in the auxiliary data attr.
    :type ds: pyasdf.ASDFDataSet
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
    # Save adjoint sources per component
    for key, adj_src in adjsrcs.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Create the standardized tag that identifies the adjoint source
            # Assumes the component is formatted properly
            adj_src_tag = "_".join([adj_src.network,
                                    adj_src.station,
                                    adj_src.component])

            # Convert the adjoint source to SPECFEM format
            srclen = len(adj_src.adjoint_source)
            specfem_adj_source = np.empty((srclen, 2))
            # Create the time axis in the 0th column
            specfem_adj_source[:, 0] = np.linspace(0, (srclen - 1) * adj_src.dt,
                                                   srclen)
            specfem_adj_source[:, 0] += time_offset
            # Time-reverse waveform
            specfem_adj_source[:, 1] = adj_src.adjoint_source[::-1]

            station_id = f"{adj_src.network}.{adj_src.station}"
            try:
                coordinates = ds.waveforms[station_id].coordinates
                # Safeguard against funny types in the coordinates dictionary
                latitude = float(coordinates["latitude"])
                longitude = float(coordinates["longitude"])
                elevation_in_m = float(coordinates["elevation_in_m"])
            except KeyError:
                latitude = longitude = elevation_in_m = 0
                logger.warning(
                    f"cannot find matching StationXML for {station_id} "
                    "coordinate information will be set to 0"
                    )
            # Parameters saved as a dictionary object
            parameters = {
                "dt": adj_src.dt, "misfit_value": adj_src.misfit,
                "adjoint_source_type": adj_src.adj_src_type,
                "min_period": adj_src.min_period,
                "max_period": adj_src.max_period,
                "latitude": latitude, "longitude": longitude,
                "elevation_in_m": elevation_in_m, "station_id": station_id,
                "component": adj_src.component, "units": "m"
            }

            ds.add_auxiliary_data(data=specfem_adj_source,
                                  data_type="AdjointSources",
                                  path=f"{path}/{adj_src_tag}",
                                  parameters=parameters
                                  )

