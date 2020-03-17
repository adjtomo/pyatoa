"""
For reading and writing files. Some functions to read the outputs of Specfem,
others to write output files for Specfem, or for any external files required
for codes that interact with Pyatoa.
"""
import os
import glob
import json
import time
import random
import numpy as np
from obspy import Stream, Trace, UTCDateTime


def read_fortran_binary(path):
    """
    Convert a Specfem3D fortran .bin file into a NumPy array,
    Copied verbatim from Seisflows/plugins/solver_io/fortran_binary.py/_read()

    :type path: str
    :param path: path to fortran .bin file
    :rtype: np.array
    :return: fortran binary data as a numpy array
    """
    nbytes = os.path.getsize(path)
    with open(path, "rb") as f:
        f.seek(0)
        n = np.fromfile(f, dtype="int32", count=1)[0]
        if n == nbytes - 8:
            f.seek(4)
            data = np.fromfile(f, dtype="float32")
            return data[:-1]
        else:
            f.seek(0)
            data = np.fromfile(f, dtype="float32")
            return data


def read_ascii(path, origintime=None, location=''):
    """
    Specfem3D outputs seismograms to ASCII (.sem?) files

    Pyatoa expects seismograms as Obspy Stream objects.
    This convenience function converts the .sem? files into Stream objects
    with the correct header information.

    Works with Specfem3D git version 6895e2f7

    Note: if origintime is None, the default start time is 1970-01-01T00:00:00,
          this is fine for quick-looking the data but it is recommended that an
          actual starttime is given.

    :type path: str
    :param path: path of the given ascii file
    :type origintime: obspy.UTCDateTime
    :param origintime: UTCDatetime object for the origintime of the event
    :type location: str
    :param location: location value for a given station/component
    :rtype st: obspy.Stream.stream
    :return st: stream containing header and data info taken from ascii file
    """
    # This was tested up to version 6895e2f7
    try:
        times = np.loadtxt(fname=path, usecols=0)
        data = np.loadtxt(fname=path, usecols=1)

    # At some point in 2018, the Specfem developers changed how the ascii files
    # were formatted from two columns to comma separated values, and repeat
    # values represented as 2*value_float where value_float represents the data
    # value as a float
    except ValueError:
        times, data = [], []
        with open(path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            try:
                time_, data_ = line.strip().split(',')
            except ValueError:
                if "*" in line:
                    time_ = data_ = line.split('*')[-1]
                else:
                    raise ValueError
            times.append(float(time_))
            data.append(float(data_))

        times = np.array(times)
        data = np.array(data)

    if origintime is None:
        print("No origintime given, setting to default 1970-01-01T00:00:00")
        origintime = UTCDateTime("1970-01-01T00:00:00")

    # We assume that dt is constant after 3 decimal points
    delta = round(times[1] - times[0], 3)

    # Honor that Specfem doesn't start exactly on 0
    origintime += times[0]

    # Write out the header information
    net, sta, cha, fmt = os.path.basename(path).split('.')
    stats = {"network": net, "station": sta, "location": location,
             "channel": cha, "starttime": origintime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'},
             "time_offset": times[0], "format": fmt
             }
    st = Stream([Trace(data=data, header=stats)])

    return st


def parse_output_optim(path_to_optim):
    """
    Seisflows creates a file 'output.optim' which provides a log of the
    optimization procedures, including trial step lengths, associated misfits.
    This function parses that file into a few arrays for use in the workflow.

    :type path_to_optim: str
    :param path_to_optim: path to 'output.optim' created by Seisflows
    :rtype: np.arrays
    :return: numpy arrays which define the iterations, steplengths and misfits
    """
    with open(path_to_optim, 'r') as f:
        lines = f.readlines()

    # Parse the file, skip the header, ignore any tail
    iterations, steplens, misfits = [], [], []
    for line in lines[2:]:
        line = line.strip().split()
        # Each iteration will have an integer to represent the iter number
        if len(line) == 3:
            iteration = line[0]
            iterations.append(int(iteration))
            steplens.append(float(line[1]))
            misfits.append(float(line[2]))
        # Each trial step length will follow and will carry the same iteration
        elif len(line) == 2:
            iterations.append(int(iteration))
            steplens.append(float(line[0]))
            misfits.append(float(line[1]))
        elif len(line) == 0:
            continue
        else:
            print(line)
            print("invalid line length encountered in output.optim")
            return

    # Set the lists as numpy arrays for easier manipulation
    iterations = np.asarray(iterations)
    steplens = np.asarray(steplens)
    misfits = np.asarray(misfits)

    return iterations, steplens, misfits


def write_misfit_json(ds, model, step, fidout="./misfits.json"):
    """
    Write a JSON file containing misfit information for a given dataset,
    model and step count. Sums misfit, number of windows, and number of
    adjoint sources for a given event, model and step.
    ---
    NOTE:
    Operates on a crude file lock system, which renames the .json file to keep
    other compute nodes from accessing a file being written.

    If this function crashes while writing, Pyatoa will get stuck in a loop,
    and Seisflows will need to be stopped and resumed.
    
    Before running this function, make sure /path/to/misfits.json_lock
    is removed or renamed.

    TO DO: introduce a better file lock system
    ---

    As per Tape (2010) Eq. 7, the total misfit function F^T is given as:
            F^T(m) = (1/S) * sum[s=1:S] (F^T_s(m))
    
    where S is the number of sources

    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset, assumed to contain auxiliary_data.Statistics
    :type model: str
    :param model: model number, e.g. "m00"
    :type step: str
    :param step: step count, e.g. "s00"
    :type fidout: str
    :param fidout: output file to write the misfit
    """
    # organize information to be written
    stats = ds.auxiliary_data.Statistics[model][step].parameters
    misfit = ds.auxiliary_data.Statistics[model][step].data[()][0]
    windows = stats["number_misfit_windows"]
    adjsrcs = stats["number_adjoint_sources"]
    event_id = os.path.basename(ds.filename).split(".")[0]
   
    # build nested dictioanaries 
    event_dict = {"misfit": float(f"{misfit:.3f}"),
                  "windows": int(windows), 
                  "adjsrcs": int(adjsrcs),
                  }
    step_dict = {event_id: event_dict}
    model_dict = {step: step_dict}
    misfit_dict = {model: model_dict}
    
    # To allow multiple CPUs clean, simultaenous write capability, create a
    # lock file that can only be accessed by one cpu at a time
    fidout_lock = fidout + "_lock"
    while True:
        # another process has control of the misfit file, wait random time
        if os.path.exists(fidout_lock):
            # print("file is locked, waiting")
            time.sleep(random.randint(5, 10))
        # misfit file is available for writing
        elif os.path.exists(fidout) and not os.path.exists(fidout_lock):
            # Lock the file so that other processes don't try to access it
            os.rename(fidout, fidout_lock)
            with open(fidout_lock, "r") as f:
                misfit_dict = json.load(f)
                # Figure out if this model/step combination is already present
                if model in misfit_dict.keys():
                    if step in misfit_dict[model].keys():
                        misfit_dict[model][step][event_id] = event_dict
                    else:
                        misfit_dict[model][step] = step_dict
                else:
                    misfit_dict[model] = model_dict
       
                # Parse dictionary to give the total misfit by summing all parts
                misfits, windows_all, adjsrcs_all = [], [], []
                for key in misfit_dict[model][step].keys():
                    if key in ["misfit", "windows", "adjsrcs"]:
                        continue
                    misfits.append(misfit_dict[model][step][key]["misfit"])
                    windows_all.append(misfit_dict[model][step][key]["windows"])
                    adjsrcs_all.append(misfit_dict[model][step][key]["adjsrcs"])

                # Save summed values into the 'step' dictionary
                misfit_dict[model][step]["misfit"] = sum(misfits)/len(misfits)
                misfit_dict[model][step]["windows"] = sum(windows_all)
                misfit_dict[model][step]["adjsrcs"] = sum(adjsrcs_all)
                f.close()
    
            # rewrite new misfit into lock file, then rename when finished
            with open(fidout_lock, "w") as f:
                json.dump(misfit_dict, f, indent=4, separators=(',', ':'), 
                          sort_keys=True)
                f.close()
            os.rename(fidout_lock, fidout)   
            return
        # misfit file has not been written yet
        else:
            # print("file has not been written")
            with open(fidout, "w") as f:
                json.dump(misfit_dict, f, indent=4, separators=(',', ':'),
                          sort_keys=True)
            return

         
def write_misfit_stats(ds, model, pathout="./", fidout=None):
    """
    A simpler alternative to write_misfit_json()

    This function simply writes a new text file for each event, which contains 
    the total misfit for that event.

    e.g. path/to/misfits/{model_number}/{event_id}
    
    These files will then need to be read by: seisflows.workflow.write_misfit()

    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset, assumed to contain auxiliary_data.Statistics
    :type model: str
    :param model: model number, e.g. "m00"
    :type pathout: str
    :param pathout: output path to write the misfit. fid will be the event name
    :type fidout: str
    :param fidout: allow user defined filename, otherwise default to name of ds
        note: if given, var 'pathout' is not used, this must be a full path
    """
    from pyatoa.utils.asdf.extractions import sum_misfits

    # By default, name the file after the name of the asdf dataset
    if fidout is None:
        event_id = os.path.basename(ds.filename).split(".")[0]
        fidout = os.path.join(pathout, event_id)
    
    # calculate misfit 
    misfit = sum_misfits(ds, model)

    # save in the same format as seisflows 
    np.savetxt(fidout, [misfit], '%11.6e')


def rcv_vtk_from_specfem(path_to_data, path_out, utm_zone=-60, z=3E3):
    """
    Creates source and receiver VTK files based on the STATIONS and
    CMTSOLUTIONS from a Specfem3D DATA directory.

    :type path_to_data: str
    :param path_to_data: path to specfem3D/DATA directory
    :type path_out: str
    :param path_out: path to save the fiels to
    :type utm_zone: int
    :param utm_zone: utm zone for converting lat lon coordinates
    :type z: float
    :param z: elevation to put stations at
    """
    from pyatoa.utils.srcrcv import lonlat_utm

    # Templates for filling
    vtk_line = "{x:18.6E}{y:18.6E}{z:18.6E}\n"
    vtk_header = """
    # vtk DataFile Version 2.0
    Source and Receiver VTK file from Pyatoa
    ASCII
    DATASET POLYDATA
    POINTS\t{} float\n
    """

    stations = np.loadtxt(os.path.join(path_to_data, "STATIONS"),
                          usecols=[2, 3], dtype=str)
    lats = stations[:, 0]
    lons = stations[:, 1]

    with open(os.path.join(path_out, "rcvs.vtk"), "w") as f:
        f.write(vtk_header.format(len(stations)))
        for lat, lon in zip(lats, lons):
            rx, ry = lonlat_utm(lon_or_x=lon, lat_or_y=lat, utm_zone=utm_zone,
                                inverse=False)
            f.write(vtk_line.format(x=rx, y=ry, z=z))


def src_vtk_from_specfem(path_to_data, path_out, utm_zone=-60, cx=None,
                         cy=None, cz=False):
    """
    Creates source and receiver VTK files based on the STATIONS and
    CMTSOLUTIONS from a Specfem3D DATA directory.

    :type path_to_data: str
    :param path_to_data: path to specfem3D/DATA directory
    :type path_out: str
    :param path_out: path to save the fiels to
    :type utm_zone: int
    :param utm_zone: utm zone for converting lat lon coordinates
    :type cx: float
    :param cx: Constant X-value for creating an Y-slice, should be in units of
        meters, in the UTM coordinate system
    :type cy: float
    :param cy: Constant Y-value for creating an X-slice, should be in units of
        meters, in the UTM coordinate system
    :type cz: float
    :param cz: Constant Z-value for creating a Z-slice, should be in units of
        meters with positve z-axis so negative Z-values are deeper
    """
    from pyatoa.utils.srcrcv import lonlat_utm

    def read_cmtsolution(path):
        """utility function to read cmtsolution file into dictionary object"""
        dict_out = {}
        cmt = np.loadtxt(path, skiprows=1, dtype="str", delimiter=":")
        for arr in cmt:
            # Replace spaces with underscores in the key
            key, value = arr[0].replace(" ", "_"), arr[1].strip()
            # Most values will be float except event name
            try:
                value = float(value)
            except ValueError:
                pass
            dict_out[key] = value
        return dict_out

    # Templates for filling
    vtk_line = "{x:18.6E}{y:18.6E}{z:18.6E}\n"
    vtk_header = """
    # vtk DataFile Version 2.0
    Source and Receiver VTK file from Pyatoa
    ASCII
    DATASET POLYDATA
    POINTS\t{} float\n
    """

    # Gather all the sources
    sources = glob.glob(os.path.join(path_to_data, "CMTSOLUTION*"))

    # Open files that need to be written
    f_xslice = f_yslice = f_zslice = None
    f_std = open(os.path.join(path_out, "srcs.vtk"))
    # Constant X-value means a slice parallel to the Y-Axis, a bit confusing
    if cx:
        f_yslice = open(os.path.join(path_out, "srcs_yslice_{cx}.vtk"))
    if cy:
        f_xslice = open(os.path.join(path_out, "srcs_xslice_{cy}.vtk"))
    if cz:
        f_zslice = open(os.path.join(path_out, "srcs_zslice_{cz}.vtk"))

    # Write in the headers, use a try-except to write all even if None
    for f in [f_std, f_xslice, f_yslice, f_zslice]:
        try:
            f.write(vtk_header.format(len(sources)))
        except AttributeError:
            continue

    # Create VTK file for Sources, assuming CMTSOLUTION format
    for source in sources:
        src = read_cmtsolution(source)
        sx, sy = lonlat_utm(lon_or_x=src["longitude"], lat_or_y=src["latitude"],
                            utm_zone=utm_zone, inverse=False)
        sz = src["depth"] * -1E3
        # Write data to all files using try-except
        f.std.write(vtk_line.format(x=sx, y=sy, z=sz))
        if cy:
            f_xslice.write(vtk_line.format(x=sx, y=cy, z=sz))
        if cx:
            f_yslice.write(vtk_line.format(x=cx, y=sy, z=sz))
        if cz:
            f_zslice.write(vtk_line.format(x=sx, y=sy, z=cz))

    # Close all the files
    for f in [f_std, f_xslice, f_yslice, f_zslice]:
        try:
            f.close()
        except AttributeError:
            continue


def create_stations_adjoint(ds, model, specfem_station_file, pathout=None):
    """
    Generate the STATIONS_ADJOINT file for Specfem input by reading in the
    STATIONS file and cross-checking which adjoint sources are available in the
    Pyasdf dataset.
    
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing AdjointSources auxiliary data
    :type model: str
    :param model: model number, e.g. "m00"
    :type specfem_station_file: str
    :param specfem_station_file: path/to/specfem/DATA/STATIONS
    :type pathout: str
    :param pathout: path to save file 'STATIONS_ADJOINT'
    """
    event_id = os.path.basename(ds.filename).split('.')[0]

    # Check which stations have adjoint sources
    stas_with_adjsrcs = []
    for code in ds.auxiliary_data.AdjointSources[model].list():
        stas_with_adjsrcs.append(code.split('_')[1])
    stas_with_adjsrcs = set(stas_with_adjsrcs)

    # Figure out which stations were simulated
    with open(specfem_station_file, "r") as f:
        lines = f.readlines()

    # If no output path is specified, save into current working directory with
    # an event_id tag to avoid confusion with other files, else normal naming
    if pathout is None:
        write_out = f"./STATIONS_ADJOINT_{event_id}"
    else:
        write_out = os.path.join(pathout, "STATIONS_ADJOINT")

    # Rewrite the Station file but only with stations that contain adjoint srcs
    with open(write_out, "w") as f:
        for line in lines:
            if line.split()[0] in stas_with_adjsrcs:
                    f.write(line)


def write_adj_src_to_ascii(ds, model, pathout=None, comp_list=["N", "E", "Z"]):
    """
    Take AdjointSource auxiliary data from a Pyasdf dataset and write out
    the adjoint sources into ascii files with proper formatting, for input
    into PyASDF.

    Note: Specfem dictates that if a station is given as an adjoint source,
        all components must be present, even if some components don't have
        any misfit windows. This function writes blank adjoint sources
        (an array of 0's) to satisfy this requirement.

    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing adjoint sources
    :type model: str
    :param model: model number, e.g. "m00"
    :type pathout: str
    :param pathout: path to write the adjoint sources to
    :type comp_list: list of str
    :param comp_list: component list to check when writing blank adjoint sources
        defaults to N, E, Z, but can also be e.g. R, T, Z
    """
    def write_to_ascii(f, array):
        """
        Function used to write the ascii in the correct format.
        Columns are formatted like the ASCII outputs of Specfem, two columns
        times written as float, amplitudes written in E notation, 6 spaces
        between.

        :type f: _io.TextIO
        :param f: the open file to write to
        :type array: numpy.ndarray
        :param array: array of data from obspy stream
        """
        for dt, amp in array:
            if dt == 0. and amp != 0.:
                dt = 0
                adj_formatter = "{dt:>13d}      {amp:13.6E}\n"
            elif dt != 0. and amp == 0.:
                amp = 0
                adj_formatter = "{dt:13.6f}      {amp:>13d}\n"
            else:
                adj_formatter = "{dt:13.6f}      {amp:13.6E}\n"

            f.write(adj_formatter.format(dt=dt, amp=amp))

    # Shortcuts
    adjsrcs = ds.auxiliary_data.AdjointSources[model]
    event_id = ds.filename.split('/')[-1].split('.')[0]

    # Set the path to write the data to.
    # If no path is given, default to current working directory
    if pathout is None:
        pathout = os.path.join("./", event_id)
    if not os.path.exists(pathout):
        os.makedirs(pathout)

    # Loop through adjoint sources and write out ascii files
    # ASDF datasets use '_' as separators but Specfem wants '.' as separators
    already_written = []
    for adj_src in adjsrcs.list():
        station = adj_src.replace('_', '.')
        fid = os.path.join(pathout, f"{station}.adj")
        with open(fid, "w") as f:
            write_to_ascii(f, adjsrcs[adj_src].data[()])

        # Write blank adjoint sources for components with no misfit windows
        for comp in comp_list:
            station_blank = (adj_src[:-1] + comp).replace('_', '.')
            if station_blank.replace('.', '_') not in adjsrcs.list() and \
                    station_blank not in already_written:
                # Use the same adjoint source, but set the data to zeros
                blank_adj_src = adjsrcs[adj_src].data[()]
                blank_adj_src[:, 1] = np.zeros(len(blank_adj_src[:, 1]))

                # Write out the blank adjoint source
                fid_blank = os.path.join(pathout, f"{station_blank}.adj")
                with open(fid_blank, "w") as b:
                    write_to_ascii(b, blank_adj_src)

                # Append to a list to make sure we don't write doubles
                already_written.append(station_blank)


def tile_combine_imgs(ds, wavs_path, maps_path, save_pdf_to,
                      sort_by="baz", purge_wavs=False, purge_tiles=False):
    """
    Utility function to tile and combine the output figures from the workflow.
    Tiles maps and waveform plots together, and then combines them into a
    composite pdf. Options to delete intermediate files. This function is useful
    for cutting down on the number of files generated by Pyatoa.

    Maps should not be purged because they can be used by future workflows

    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset used for sorting
    :type wavs_path: str
    :param wavs_path: path to the waveform figures
    :type maps_path: str
    :param maps_path: path to the map figures
    :type save_pdf_to: str
    :param save_pdf_to: path and filename to save final PDF file
    :type sort_by: str
    :param sort_by: method to sort stations by when combining into a pdf,
        available: 'baz', 'misfit' (misfit not implemented)
    :type purge_wavs: bool
    :param purge_wavs: delete the waveform files after tiling them
    :type purge_tiles: bool
    :param purge_tiles: delete the tile files after combining into pdf
    :return:
    """
    # Intra-function imports because this is usually only called once in a while
    from pyatoa.visuals.plot_tools import tile_imgs, imgs_to_pdf
    from pyatoa.utils.srcrcv import sort_by_backazimuth

    # Set the template filenames to look for/ use
    wav_fid = "wav_{sta}.png"
    map_fid = "map_{sta}.png"
    tile_fid = "tile_{sta}.png"

    # Get station names from waveforms. Maps will be named similarly
    files = glob.glob(os.path.join(wavs_path, wav_fid.format(sta="*")))
    stanames = []
    for f in files:
        sta = os.path.basename(f).split("_")[1].split(".")[0]
        stanames.append(sta)
    stanames = set(stanames)
    stanames = list(stanames)

    # combine map and waveform figures into tiles
    tile_names = []
    for sta in stanames:
        wav_name = os.path.join(wavs_path, wav_fid.format(sta=sta))
        map_name = os.path.join(maps_path, map_fid.format(sta=sta))
        tile_name = os.path.join(wavs_path, tile_fid.format(sta=sta))

        if os.path.exists(map_name) and os.path.exists(wav_name):
            tile_imgs(fids=[map_name, wav_name], fid_out=tile_name)
            tile_names.append(tile_name)
        else:
            raise FileNotFoundError(
                f"Either {wav_name} or {map_name} doesn't exist when it should")

    # remove old waveform images
    if purge_wavs:
        for f in files:
            os.remove(f)

    # combine the tiles into a single .pdf
    # sort stations by backazimuth for easier visualization
    if sort_by:
        if sort_by == "baz":
            sorted_station_names = sort_by_backazimuth(ds)
        # Sort by largest to smallest misfit
        elif sort_by == "misfit":
            raise NotImplementedError

        # Sort the tile names based on the sort method
        tile_names_sorted = []
        for name in sorted_station_names:
            net, sta = name.split('.')
            for tile in tile_names:
                if sta in tile:
                    tile_names_sorted.append(tile)
        tile_names = tile_names_sorted

    # Create the pdf
    imgs_to_pdf(fids=tile_names, fid_out=save_pdf_to)

    # Remove tiles
    if purge_tiles:
        for tile in tile_names:
            os.remove(tile)


def srcrcv_vtk_from_dataset(pathin, pathout, model, utm_zone=-60,
                            event_depth=False):
    """
    !!! Deprecated in favor of srcrcv_vtk_from_specfem !!!

    Same as create_srcrcv_vtk_single, except instead of taking an asdf
    dataset input, takes a path, reads in datasets and creates one
    large vtk file containing all stations and all events.

    -Useful for visualizations of misfit kernels and gradients.
    -Automatically creates a separate event vtk file.

    :type pathin: str
    :param pathin: path containing .h5 files, will loop through all
        available h5 files in the folder
    :type model: str
    :param model: model number, e.g. 'm00'
    :type pathout: str
    :param pathout: output path to save vtk file
    :type utm_zone: int
    :param utm_zone: the utm zone of the mesh, 60 for NZ
    :type event_depth: bool
    :param event_depth: if True, uses the real event depth, if False, places
        event at 5km above the surface
    """
    import pyasdf
    from pyatoa.utils.srcrcv import lonlat_utm

    vtk_header = ("# vtk DataFile Version 2.0\n"
                  "Source and Receiver VTK file from Pyatoa\n"
                  "ASCII\n"
                  "DATASET POLYDATA\n"
                  "POINTS\t{} float\n"
                  )
    vtk_line = "{X:18.6E}{Y:18.6E}{E:18.6E}\n"

    # Loop through available datasets
    datasets = glob.glob(os.path.join(pathin, '*.h5'))
    if not datasets:
        return

    event_ids, sta_ids = [], []
    ev_x, ev_y, ev_z, sta_x, sta_y, sta_elv = [], [], [], [], [], []
    for fid in datasets:
        with pyasdf.ASDFDataSet(fid) as ds:
            # Check if dataset contains adjoint sources
            if not bool(ds.auxiliary_data.AdjointSources):
                continue

            # Loop through stations with adjoint sources
            if hasattr(ds.auxiliary_data.AdjointSources, model):
                for adjsrc in ds.auxiliary_data.AdjointSources[model].list():
                    sta = ds.auxiliary_data.AdjointSources[model][adjsrc]

                    # make sure no repeat stations
                    if sta.parameters["station_id"] in sta_ids:
                        continue

                    # Convert lat lon to UTM
                    sx, sy = lonlat_utm(lon_or_x=sta.parameters["longitude"],
                                        lat_or_y=sta.parameters["latitude"],
                                        utm_zone=utm_zone, inverse=False
                                        )
                    sta_x.append(sx)
                    sta_y.append(sy)
                    sta_elv.append(sta.parameters["elevation_in_m"])
                    sta_ids.append(sta.parameters["station_id"])
            else:
                continue

            # Get event location information in UTM
            event_id = os.path.basename(ds.filename).split(".")[0]
            ex, ey = lonlat_utm(
                lon_or_x=ds.events[0].preferred_origin().longitude,
                lat_or_y=ds.events[0].preferred_origin().latitude,
                utm_zone=utm_zone, inverse=False
            )
            # Depth in units of meters
            if event_depth:
                ez = ds.events[0].preferred_origin().depth
            else:
                # set event epicentral depth to 5km to keep it above topography
                ez = 5E3

            event_ids.append(event_id)
            ev_x.append(ex)
            ev_y.append(ey)
            ev_z.append(ez)

    # Write header for VTK file and then print values for source receivers
    fid_out = os.path.join(pathout, f"rcvs_{model}.vtk")
    with open(fid_out, "w") as f:
        f.write(vtk_header.format(len(sta_x)))
        # Loop through stations and write them to vtk file
        for sx, sy, se in zip(sta_x, sta_y, sta_elv):
            f.write(vtk_line.format(X=sx, Y=sy, E=se))

    # Make a separate VTK file for all events so they can be formatted different
    event_fid_out = os.path.join(pathout, "srcs.vtk")
    if not os.path.exists(event_fid_out):
        with open(event_fid_out, "w") as f:
            f.write(vtk_header.format(len(ev_x)))
            for ex, ey, ez in zip(ev_x, ev_y, ev_z):
                f.write(vtk_line.format(X=ex, Y=ey, E=ez))

    # Make a separate VTK file for each event.
    # This only needs to be run once so just skip over if the files exist
    for event_id, ex, ey, ez in zip(event_ids, ev_x, ev_y, ev_z):
        event_fid_out = os.path.join(pathout, f"{event_id}.vtk")
        if os.path.exists(event_fid_out):
            continue
        with open(event_fid_out, "w") as f:
            f.write(vtk_header.format(1))
            f.write(vtk_line.format(X=ex, Y=ey, E=ez))


def create_srcrcv_vtk_single(ds, model, pathout, event_separate=False,
                             utm_zone=-60, event_depth=False):
    """
    !!! Deprectated, I don't really use this !!!

    It's useful to visualize source receiver locations in Paraview, alongside
    sensitivity kernels. VTK files are produced by Specfem, however this sr.vtk
    file contains all receivers, and a source at depth which is sometimes
    confusing for top down visualization.

    NOTE:
    -This function will create source_receiver vtk files using the .h5 files
    with only those receiver that were used in the misfit analysis, and only
    an epicentral source location.
    -Gives the option to create an event vtk file separate to receivers, for
    more flexibility in the visualization.

    :type ds: pyasdf.ASDFDataSet
    :param ds: pyasdf dataset outputted by pyatoa
    :type model: str
    :param model: model number, e.g. 'm00'
    :type pathout: str
    :param pathout: output path to save vtk file
    :type event_separate: str
    :param event_separate: if event vtk file to be made separately
    :type utm_zone: int
    :param utm_zone: the utm zone of the mesh, 60 for NZ
    :type event_depth: bool
    :param event_depth: if True, uses the real event depth, if False, places
        event at 5km above the surface
    """
    from warnings import warn
    warn("Deprecated", DeprecationWarning)

    from pyatoa.utils.srcrcv import lonlat_utm

    # Check that this can be run, if dataset contains adjoint sources
    if not bool(ds.auxiliary_data.AdjointSources):
        return

    # Some information that is used a few times
    vtk_header = ("# vtk DataFile Version 2.0\n"
                  "Source and Receiver VTK file from Pyatoa\n"
                  "ASCII\n"
                  "DATASET POLYDATA\n"
                  )

    # Get receiver location information in lat-lon,
    event_id = os.path.basename(ds.filename).split(".")[0]
    sta_x, sta_y, sta_elv, sta_ids = [], [], [], []

    for adjsrc in ds.auxiliary_data.AdjointSources[model].list():
        sta = ds.auxiliary_data.AdjointSources[model][adjsrc]

        # make sure no repeat stations
        if sta.parameters["station_id"] in sta_ids:
            continue

        # Convert lat lon to UTM
        x, y = lonlat_utm(lon_or_x=sta.parameters["longitude"],
                          lat_or_y=sta.parameters["latitude"],
                          utm_zone=utm_zone, inverse=False)
        sta_x.append(x)
        sta_y.append(y)
        sta_elv.append(sta.parameters["elevation_in_m"])
        sta_ids.append(sta.parameters["station_id"])

    # Get event location information in UTM
    ev_x, ev_y = lonlat_utm(lon_or_x=ds.events[0].preferred_origin().longitude,
                            lat_or_y=ds.events[0].preferred_origin().latitude,
                            utm_zone=utm_zone, inverse=False
                            )
    # Depth in units of meters
    if event_depth:
        ev_z = ds.events[0].preferred_origin().depth
    else:
        # set event epicentral depth to 5km to keep it above topography
        ev_z = 5E3

    # Write header for VTK file and then print values for source receivers
    fid_out = os.path.join(pathout, f"{event_id}_{model}.vtk")
    with open(fid_out, "w") as f:
        f.write(vtk_header)
        # num points equal to number of stations plus 1 event
        f.write(f"POINTS\t{len(sta_x)+1} float\n")
        f.write(f"{ev_x:18.6E}{ev_y:18.6E}{ev_z:18.6E}\n")
        for x, y, e in zip(sta_x, sta_y, sta_elv):
            f.write(f"{x:18.6E}{y:18.6E}{e:18.6E}\n")

    # Make a separate VTK file for the source
    if event_separate:
        event_fid_out = os.path.join(pathout, f"{event_id}_{model}_event.vtk")
        with open(event_fid_out, "w") as f:
            f.write(vtk_header)
            f.write("POINTS\t1 float\n")
            f.write(f"{ev_x:18.6E}{ev_y:18.6E}{ev_z:18.6E}\n")
