"""
For writing various output files used by Pyatoa, Specfem and Seisflows
"""
import os
import glob
import numpy as np
from obspy.core.inventory.channel import Channel
from pyatoa import logger
from pyatoa.utils.form import (format_event_name, format_iter, format_step, 
                               channel_code)


def write_stations(inv, fid="./STATIONS", elevation=False, burial=0.):
    """
    Write a SPECFEM3D STATIONS file given an ObsPy inventory object

    .. note::
        If topography is implemented in your mesh, elevation values should be
        set to 0 which means 'surface' in SPECFEM.

    .. note::
        burial information is not contained in the ObsPy inventory so it is
        always set to a constant value, which can be given as input. default 0

    :type inv: obspy.core.inventory.Inventory
    :param inv: Inventory object with station locations to write
    :type fid: str
    :param fid: path and file id to save the STATIONS file to.
    :type elevation: bool
    :param elevation: if True, sets the actual elevation value from the
        inventory. if False, sets elevation to 0.
    :type burial: float
    :param burial: the constant value to set burial values to. defaults to 0.
    """
    with open(fid, "w") as f:
        for net in inv:
            for sta in net:
                lat = sta.latitude
                lon = sta.longitude
                if elevation:
                    elv = sta.elevation
                else:
                    elv = 0.

                f.write(f"{sta.code:>6}{net.code:>6}{lat:12.4f}{lon:12.4f}"
                        f"{elv:7.1f}{burial:7.1f}\n")


def write_inv_seed(inv, path="./", dir_structure="{sta}.{net}",
                   file_template="RESP.{net}.{sta}.{loc}.{cha}",
                   components="ZNE", channel_code="HX{comp}", **kwargs):
    """
    Pyatoa requires stations to be discoverable in SEED format, i.e., in a data
    center repository. This structure dictates that each component of each 
    station has its own individual StationXML file, saved in a specific 
    directory structure with a unique file naming schema. 

    This utility is useful for creating the necessary StationXML files for 
    temporary or synthetic stations which are not discoverable via
    FDSN or through datacenters.

    .. note::
        kwargs are passed to obspy.core.inventory.channel.Channel

    :type path: str
    :param path: location to save StationXML files to
    :type dir_structure: str
    :param dir_structure: template for directory structure, likely should not 
        need to change from the default 
    :type file_template: str
    :param file_template: template for file naming, likely should not change 
        from default template
    :type channels: str
    :param channels: OPTIONAL, if inventory does not contain channels (e.g.,
        if read from a SPECFEM STATIONS file), channels will be generated here.
    :type channel_code: str
    :param channel_code: Explicitely defined default channel values for
    generating channels on the fly when none are provided by the inventory
    """
    # Default values for other values required by the inventory
    location_code = kwargs.get("location_code", "")
    elevation = kwargs.get("elevation", 0.)
    depth = kwargs.get("depth", 0.)

    assert(os.path.exists(path)), f"output path does not exist: {path}"
    
    for net in inv:
        for sta in net:
            # Create the directory to store individual channels
            sta_dir = os.path.join(path, dir_structure.format(sta=sta.code, 
                                                              net=net.code)
                                   )
            if not os.path.exists(sta_dir):
                os.makedirs(sta_dir)
          
            # If the station has no channels inherently, generate them on the 
            # fly based on default and user-defined information 
            if len(sta) == 0:
                channel_list = []
                for comp in components:
                    cha_ = Channel(code=channel_code.format(comp=comp),
                                   location_code=location_code,
                                   latitude=sta.latitude, 
                                   longitude=sta.longitude, elevation=elevation,
                                   depth=depth, **kwargs
                                   )
                    channel_list.append(cha_)
                sta.channels = channel_list

            # Cycle through the channels and generate individual StationXMLs
            for cha in sta:
                # Select the channel out of the inventory
                channel_inv = inv.select(network=net.code, station=sta.code,
                                         channel=cha.code
                                         )
                # If select returns properly, write out channel as StationXML
                if channel_inv:
                    fid_out = os.path.join(
                        path, 
                        dir_structure.format(sta=sta.code, net=net.code),
                        file_template.format(net=net.code, sta=sta.code,
                                             loc=cha.location_code,
                                             cha=cha.code)
                    )
                    channel_inv.write(fid_out, format="STATIONXML")
                else:
                    print("{}.{} could not be selected".format(sta.code,
                                                               cha.code)
                          )


def write_sem(st, unit, path="./", time_offset=0):
    """
    Write an ObsPy Stream in the two-column ASCII format that Specfem uses

    :type st: obspy.core.stream.Stream
    :param st: stream containing synthetic data to be written
    :type unit: str
    :param unit: the units of the synthetic data, used for file extension, must 
        be 'd', 'v', 'a' for displacement, velocity, acceleration, resp.
    :type path: str
    :param path: path to save data to, defaults to cwd
    :type time_offset: float
    :param time_offset: optional argument to offset the time array. Sign matters
        e.g. time_offset=-20 means t0=-20
    """
    assert(unit.lower() in ["d", "v", "a"]), "'unit' must be 'd', 'v' or 'a'"
    for tr in st:
        s = tr.stats
        fid = f"{s.network}.{s.station}.{channel_code(s.delta)}X{s.channel[-1]}"
        fid = os.path.join(path, f"{fid}.sem{unit.lower()}")
        data = np.vstack((tr.times() + time_offset, tr.data)).T
        np.savetxt(fid, data, ["%13.7f", "%17.7f"])


def write_misfit(ds, iteration, step_count=None, path="./", fidout=None):
    """
    This function writes a text file containing event misfit.
    This misfit value corresponds to F_S^T of Eq 6. Tape et al. (2010)

    e.g. path/to/misfits/{iteration}/{event_id}
    
    These files will then need to be read by: seisflows.workflow.write_misfit()

    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset, assumed to contain auxiliary_data.Statistics
    :type iteration: str or int
    :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
    :type step_count: str or int
    :param step_count: step count e.g. "s00". Will be formatted so int ok.
    :type path: str
    :param path: output path to write the misfit. fid will be the event name
    :type fidout: str
    :param fidout: allow user defined filename, otherwise default to name of ds
        note: if given, var 'pathout' is not used, this must be a full path
    """
    iter_tag = format_iter(iteration)
    step_tag = format_step(step_count)

    # By default, name the file after the event id
    if fidout is None:
        fidout = os.path.join(path, format_event_name(ds))
    
    # Collect the total misfit calculated by Pyadjoint
    total_misfit = 0
    adjoint_sources = ds.auxiliary_data.AdjointSources[iter_tag]
    if step_tag:
        adjoint_sources = adjoint_sources[step_tag]

    for adjsrc in adjoint_sources.list():
        total_misfit += adjoint_sources[adjsrc].parameters["misfit"]

    # Count up the number of misfit windows
    win = ds.auxiliary_data.MisfitWindows[iter_tag]
    if step_tag:
        win = win[step_tag]
    number_windows = len(win)

    scaled_misfit = 0.5 * total_misfit / number_windows

    # save in the same format as seisflows 
    np.savetxt(fidout, [scaled_misfit], '%11.6e')

    return scaled_misfit


def write_stations_adjoint(ds, iteration, specfem_station_file, step_count=None,
                           pathout=None):
    """
    Generate the STATIONS_ADJOINT file for Specfem input by reading in the
    STATIONS file and cross-checking which adjoint sources are available in the
    Pyasdf dataset.
    
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing AdjointSources auxiliary data
    :type iteration: str or int
    :param iteration: iteration number, e.g. "i01". Will be formatted so int ok.
    :type step_count: str or int
    :param step_count: step count e.g. "s00". Will be formatted so int ok.
        If NoneType, final step of the iteration will be chosen automatically.
    :type specfem_station_file: str
    :param specfem_station_file: path/to/specfem/DATA/STATIONS
    :type pathout: str
    :param pathout: path to save file 'STATIONS_ADJOINT'
    """
    # Check which stations have adjoint sources
    stas_with_adjsrcs = []
    adj_srcs = ds.auxiliary_data.AdjointSources[format_iter(iteration)]
    # Dynamically determine final step count in the iteration
    if step_count is None:
        step_count = adj_srcs.list()[-1]
    logger.debug(f"writing stations adjoint for "
                f"{format_iter(iteration)}{format_step(step_count)}"
                )
    adj_srcs = adj_srcs[format_step(step_count)]

    for code in adj_srcs.list():
        stas_with_adjsrcs.append(code.split('_')[1])
    stas_with_adjsrcs = set(stas_with_adjsrcs)

    # Figure out which stations were simulated
    with open(specfem_station_file, "r") as f:
        lines = f.readlines()

    # If no output path is specified, save into current working directory with
    # an event_id tag to avoid confusion with other files, else normal naming
    if pathout is None:
        write_out = f"./STATIONS_ADJOINT_{format_event_name(ds)}"
    else:
        write_out = os.path.join(pathout, "STATIONS_ADJOINT")

    # Rewrite the Station file but only with stations that contain adjoint srcs
    with open(write_out, "w") as f:
        for line in lines:
            if line.split()[0] in stas_with_adjsrcs:
                    f.write(line)


def write_adj_src_to_ascii(ds, iteration, step_count=None, pathout=None, 
                           comp_list="ZNE"):
    """
    Take AdjointSource auxiliary data from a Pyasdf dataset and write out
    the adjoint sources into ascii files with proper formatting, for input
    into PyASDF.

    .. note::
        Specfem dictates that if a station is given as an adjoint source,
        all components must be present, even if some components don't have
        any misfit windows. This function writes blank adjoint sources
        (an array of 0's) to satisfy this requirement.

    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing adjoint sources
    :type iteration: str or int
    :param iteration: iteration number, e.g. "i00". Will be formatted so int ok.
    :type step_count: str or int
    :param step_count: step count e.g. "s00". Will be formatted so int ok.
            If NoneType, final step of the iteration will be chosen automatically.
    :type pathout: str
    :param pathout: path to write the adjoint sources to
    :type comp_list: str
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
    adjsrcs = ds.auxiliary_data.AdjointSources[format_iter(iteration)]
    if step_count is None:
        step_count = adjsrcs.list()[-1]
    adjsrcs = adjsrcs[format_step(step_count)]
    logger.debug(f"writing adjoint sources to ascii for "
                f"{format_iter(iteration)}{format_step(step_count)}"
                )

    # Set the path to write the data to.
    # If no path is given, default to current working directory
    if pathout is None:
        pathout = os.path.join("./", format_event_name(ds))
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
        for comp in list(comp_list):
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

    def quick_read_cmtsolution(path):
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
        src = quick_read_cmtsolution(source)
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


