"""
For writing various output files used by Pyatoa, Specfem and Seisflows
"""
import os
import glob
import numpy as np
from pyatoa.utils.form import format_event_name


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


def write_misfit(ds, model, step, path="./", fidout=None):
    """
    This function writes a text file containing event misfit.
    This misfit value corresponds to F_S^T of Eq 6. Tape et al. (2010)

    e.g. path/to/misfits/{model_number}/{event_id}
    
    These files will then need to be read by: seisflows.workflow.write_misfit()

    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset, assumed to contain auxiliary_data.Statistics
    :type model: str
    :param model: model number, e.g. "m00"
    :type step: str
    :param step: step count, e.g. "s00"
    :type path: str
    :param path: output path to write the misfit. fid will be the event name
    :type fidout: str
    :param fidout: allow user defined filename, otherwise default to name of ds
        note: if given, var 'pathout' is not used, this must be a full path
    """
    from pyatoa.utils.asdf.fetch import sum_misfits

    # By default, name the file after the event id
    if fidout is None:
        fidout = os.path.join(path, event_name(ds=ds))
    
    # calculate misfit 
    misfit = sum_misfits(ds, model, step)

    # save in the same format as seisflows 
    np.savetxt(fidout, [misfit], '%11.6e')


def rcv_vtk_from_specfem(path_to_data, path_out="./", utm_zone=-60, z=3E3):
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
    vtk_header = ("# vtk DataFile Version 2.0\n" 
                  "Source and Receiver VTK file from Pyatoa\n"
                  "ASCII\n"
                  "DATASET POLYDATA\n"
                  "POINTS\t{} float\n")

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
        f.write("\n")


def src_vtk_from_specfem(path_to_data, path_out="./", utm_zone=-60, cx=None,
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
    vtk_header = ("# vtk DataFile Version 2.0\n" 
                  "Source and Receiver VTK file from Pyatoa\n"
                  "ASCII\n"
                  "DATASET POLYDATA\n"
                  "POINTS\t{} float\n")

    # Gather all the sources
    sources = glob.glob(os.path.join(path_to_data, "CMTSOLUTION*"))

    # Open files that need to be written
    f_std = open(os.path.join(path_out, "srcs.vtk"), "w")
    f_xslice = f_yslice = f_zslice = None
    # Constant X-value means a slice parallel to the Y-Axis, a bit confusing
    if cx:
        f_yslice = open(os.path.join(path_out, f"srcs_yslice_{cx}.vtk"), "w")
    if cy:
        f_xslice = open(os.path.join(path_out, f"srcs_xslice_{cy}.vtk"), "w")
    if cz:
        f_zslice = open(os.path.join(path_out, f"srcs_zslice_{abs(cz)}.vtk"),
                        "w")

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
        f_std.write(vtk_line.format(x=sx, y=sy, z=sz))
        if cy:
            f_xslice.write(vtk_line.format(x=sx, y=cy, z=sz))
        if cx:
            f_yslice.write(vtk_line.format(x=cx, y=sy, z=sz))
        if cz:
            f_zslice.write(vtk_line.format(x=sx, y=sy, z=cz))

    # Close all the files
    for f in [f_std, f_xslice, f_yslice, f_zslice]:
        try:
            f.write("\n")
            f.close()
        except AttributeError:
            continue


def write_stations_adjoint(ds, model, specfem_station_file, step=None,
                           pathout=None):
    """
    Generate the STATIONS_ADJOINT file for Specfem input by reading in the
    STATIONS file and cross-checking which adjoint sources are available in the
    Pyasdf dataset.
    
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing AdjointSources auxiliary data
    :type model: str
    :param model: model number, e.g. "m00"
    :type step: str
    :param step: step count, e.g. "s00"
    :type specfem_station_file: str
    :param specfem_station_file: path/to/specfem/DATA/STATIONS
    :type pathout: str
    :param pathout: path to save file 'STATIONS_ADJOINT'
    """
    eid = event_name(ds=ds)

    # Check which stations have adjoint sources
    stas_with_adjsrcs = []
    adj_srcs = ds.auxiliary_data.AdjointSources[model]
    if step:
        adj_srcs = adj_srcs[step]

    for code in adj_srcs.list():
        stas_with_adjsrcs.append(code.split('_')[1])
    stas_with_adjsrcs = set(stas_with_adjsrcs)

    # Figure out which stations were simulated
    with open(specfem_station_file, "r") as f:
        lines = f.readlines()

    # If no output path is specified, save into current working directory with
    # an event_id tag to avoid confusion with other files, else normal naming
    if pathout is None:
        write_out = f"./STATIONS_ADJOINT_{eid}"
    else:
        write_out = os.path.join(pathout, "STATIONS_ADJOINT")

    # Rewrite the Station file but only with stations that contain adjoint srcs
    with open(write_out, "w") as f:
        for line in lines:
            if line.split()[0] in stas_with_adjsrcs:
                    f.write(line)


def write_adj_src_to_ascii(ds, model, step=None, pathout=None, 
                           comp_list=["N", "E", "Z"]):
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
    :type step: str
    :param step: step count e.g. "s00"
    :type pathout: str
    :param pathout: path to write the adjoint sources to
    :type comp_list: list of str
    :param comp_list: component list to check when writing blank adjoint sources
        defaults to N, E, Z, but can also be e.g. R, T, Z
    """
    def write_to_ascii(f_, array):
        """
        Function used to write the ascii in the correct format.
        Columns are formatted like the ASCII outputs of Specfem, two columns
        times written as float, amplitudes written in E notation, 6 spaces
        between.

        :type f_: _io.TextIO
        :param f_: the open file to write to
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

            f_.write(adj_formatter.format(dt=dt, amp=amp))

    # Shortcuts
    adjsrcs = ds.auxiliary_data.AdjointSources[model]
    if step:
        adjsrcs = adjsrcs[step]

    eid = event_name(ds=ds)

    # Set the path to write the data to.
    # If no path is given, default to current working directory
    if pathout is None:
        pathout = os.path.join("./", eid)
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

