"""
ASDF Datasets can be given auxiliary data to supplement the existing waveform,
event and station information contained. The functions contained in this script
add new auxiliary data structures to existing ASDF datasets
"""
import numpy as np


def write_stats_to_asdf(ds, model, step):
    """
    Misfit stats are useful for quickly investigating the information contained
    in the ASDF Datasets. Calls on functions from pyatoa.utils.asdf.extraction
    to extract misfit information and then writes it into auxiliary_data
    The data of auxiliary_data.Statistics is the total misfit for a given
     
    :type ds: pyasdf.asdf_data_set.ASDFDataSet
    :param ds: asdf dataset containing Statistics auxiliary data
    :type model: str
    :param model: model tag, e.g. "m00"
    :type step: int
    :param step: trial step number (from Seisflows)
    """
    import pyatoa.utils.asdf.extractions as extract
    
    stats = extract.misfit_stats(ds, model)
    misfit = extract.sum_misfits(ds, model)
   
    ds.add_auxiliary_data(data=np.array([misfit]), 
                          data_type="Statistics",
                          path="{mod}/s{step:0>2}".format(mod=model, step=step),
                          parameters=stats
                          )
    

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
