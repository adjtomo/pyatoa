"""
Functions to convert ASDFDataSet into individual data files that can be 
stored in a directory structure. Not very smart, just dumps all data into a 
given path. 
"""
import os
import json
import numpy as np
from pyatoa.utils.form import format_event_name
from pyatoa.utils.write import write_adj_src_to_ascii


def write_all(ds, path="./"):
    """
    Convenience function to dump everything inside a dataset

    :type ds: pyasdf.ASDFDataSet
    :param ds: Dataset containing info the write
    :type path: str
    :param path: path to save data to
    """
    write_events(ds, path)
    write_stations(ds, path)
    write_waveforms(ds, path)
    write_windows(ds, path)
    write_adjoint_sources(ds, path)


def write_events(ds, path="./"):
    """
    Write Event object as a QuakeML file.

    :type ds: pyasdf.ASDFDataSet
    :param ds: Dataset containing info the write
    :type path: str
    :param path: path to save data to
    """
    for event in ds.events:
        event.write(os.path.join(path, f"{format_event_name(event)}.xml"),
                    format="QUAKEML"
                    )


def write_stations(ds, path="./"):
    """
    Write Stations dataless files as STATIONXML files

    :type ds: pyasdf.ASDFDataSet
    :param ds: Dataset containing info the write
    :type path: str
    :param path: path to save data to
    """
    for sta_name in ds.waveforms.list():
        ds.waveforms[sta_name].StationXML.write(
            os.path.join(path, f"{sta_name.replace('.','_')}.xml"), 
            format="STATIONXML")


def write_waveforms(ds, path="./", station=None, tag=None):
    """
    Write waveforms as MSEED files

    :type ds: pyasdf.ASDFDataSet
    :param ds: Dataset containing info the write
    :type path: str
    :param path: path to save data to
    """
    # Set up the directory structure
    for sta in ds.waveforms.list():
        if station is not None and sta != station:
            continue
        for tag_ in ds.waveforms[sta].get_waveform_tags():
            if tag is not None and tag_ != tag:
                continue
            st = ds.waveforms[sta][tag_]
            st.write(os.path.join(path, f"{sta.replace('.','_')}_{tag_}.ms"), 
                     format="MSEED")


def write_windows(ds, path="./"):
    """
    Write MisfitWindows as .JSON files

    :type ds: pyasdf.ASDFDataSet
    :param ds: Dataset containing info the write
    :type path: str
    :param path: path to save data to
    """
    class WindowEncoder(json.JSONEncoder):
        """
        So that JSON plays nice with numpy objects, UTCDateTimes are already str
        taken from Pyflex.window_selector.write()
        """
        def default(self, obj):
            # Numpy objects also require explicit handling.
            if isinstance(obj, np.int64):
                return int(obj)
            elif isinstance(obj, np.int32):
                return int(obj)
            elif isinstance(obj, np.float64):
                return float(obj)
            elif isinstance(obj, np.float32):
                return float(obj)
            # Let the base class default method raise the TypeError
            return json.JSONEncoder.default(self, obj)

    windows = ds.auxiliary_data.MisfitWindows
    for model in windows.list():
        for step in windows[model].list():
            window_dict = {}
            for win in windows[model][step].list():
                window_dict[win] = windows[model][step][win].parameters
            with open(os.path.join(path, 
                      f"windows_{model}{step}.json"), "w") as f:
                json.dump(window_dict, f, cls=WindowEncoder, indent=4, 
                          separators=(',', ':')
                          )   


def write_adjoint_sources(ds, path="./"):
    """
    Write AdjointSources as ASCII files in directories corresponding to model
    number and step count.

    :type ds: pyasdf.ASDFDataSet
    :param ds: Dataset containing info the write
    :type path: str
    :param path: path to save data to
    """
    adjsrcs = ds.auxiliary_data.AdjointSources
    for model in adjsrcs.list():
        for step in adjsrcs[model].list():
            pathout = os.path.join(path, f"{model}{step}")
            if not os.path.exists(pathout):
                os.makedirs(pathout)
            write_adj_src_to_ascii(ds, model, step, pathout)






