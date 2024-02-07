"""
ASDF Datasets can be given auxiliary data to supplement the existing waveform,
event and station information contained. The functions contained in this script
add new auxiliary data structures to existing ASDF datasets
"""
import numpy as np
from pyatoa import logger
from pyatoa.utils.asdf.clean import del_auxiliary_data_path


def add_misfit_windows(windows, ds, path, overwrite=True, 
                       _data_type="MisfitWindows"):
    """
    Write Pyflex misfit windows into the 'MisfitWindows' attribute of the 
    auxiliary data of an ASDFDataSet.  Misfit windows will be tagged by the 
    `path`, the receiver name, component, and window number.

    :type windows: dict of list of pyflex.Window
    :param windows: dictionary of lists of window objects with keys
        corresponding to components related to each window
    :type ds: pyasdf.ASDFDataSet
    :param ds: ASDF data set to save windows to
    :type path: str
    :param path: internal pathing to save location of auxiliary data, separated
        by delimeter '/' for nesting. E.g., i01/s00 will save misfit windows to 
        ASDFDataSet.auxiliary_data.MisfitWindows.i01.s00. Required by PyASDF
    :type overwrite: bool
    :param overwrite: if auxiliary data for the automatically generated path 
        already exists in the dataset, deletes the existing misfit window.
        If False, ASDFWarning will be thrown and new misfit windows will not 
        be written.
    :type _data_type: str
    :param _data_type: the main 'directory' under AuxiliaryData that misfit 
        windows are stored under. By default, Pyatoa expects this to be 
        'MisfitWindows', but it can be changed by advanced users, changing this
        will likely have unintended downstream consequences.
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

            # If overwrite and windows with this particular path already exist
            # then delete them so new windows can be written
            if overwrite:
                try:
                    fullpath = "/".join[_data_type, path, window_tag]
                    del_auxiliary_data_path(ds=ds, path=fullpath)
                    logger.debug(f"overwriting existing windows {fullpath}")
                except KeyError:
                    continue

            # Write windows into dataset
            ds.add_auxiliary_data(data=np.array([win.left, win.right]),
                                  data_type=_data_type, parameters=wdict,
                                  path=f"{path}/{window_tag}")


def add_adjoint_sources(adjsrcs, ds, path, time_offset, overwrite=True):
    """
    Writes the adjoint source to an ASDF file.

    .. warning::
        It is inherently assumed SPECFEM will be using the adjoint source

    .. note::
        Modified from Pyadjoint source code:
        pyadjoint.adjoint_source.write_to_asdf()

    .. note::
        SPECFEM requires one additional parameter: the temporal offset of the
        first sample in seconds. This will have been set by the parameter
        USER_T0 in the constants.h.in file of SPECFEM's setup directory

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
    :type overwrite: bool
    :param overwrite: if auxiliary data for the automatically generated path 
        already exists in the dataset, deletes the existing adjoint source.
        If False, ASDFWarning will be thrown and new adjoint source will not 
        be written.
    """
    # Save adjoint sources per component
    for key, adj_src in adjsrcs.items():
        # Create the standardized tag that identifies the adjoint source
        # Assumes the component is formatted properly by the Manager
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

        # Parameters saved as a dictionary object to match the variables of
        # the AdjointSource object, with additional identifiers
        parameters = {"adj_src_type": adj_src.adjsrc_type,
                        "misfit": adj_src.misfit,
                        "dt": adj_src.dt,
                        "component": adj_src.component,
                        "min_period": adj_src.min_period or "None",
                        "max_period": adj_src.max_period or "None",
                        "network": adj_src.network,
                        "station": adj_src.station,
                        "location": adj_src.location,
                        "starttime": str(adj_src.starttime)
                        }
        ds.add_auxiliary_data(data=specfem_adj_source.astype(np.float64),
                                data_type="AdjointSources",
                                path=f"{path}/{adj_src_tag}",
                                parameters=parameters
                                )

