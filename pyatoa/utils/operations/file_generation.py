"""
For generation of input files for Specfem, or for any external files required
for codes that interact with Pyatoa
"""
import os
import json
import time

def write_misfit_json(ds, model, step_count=0, fidout="./misfits.json"):
    """
    Misfit text file useful for quickly determining misfit information garnered
    from a swath of pyatoa runs, used by Seisflows

    As per Tape (2010) Eq. 7, the total misfit function F^T is given as:
            F^T(m) = (1/S) * sum[s=1:S] (F^T_s(m))
    
    where S is the number of sources

    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset, assumed to contain auxiliary_data.Statistics
    :type model: str
    :param model: model number, e.g. "m00"
    :type step: int
    :param step: line search step count
    :type fidout: str
    :param fidout: output file to write the misfit
    """
    # organize information to be written
    step = "s{:0>2}".format(step_count)
    stats = ds.auxiliary_data.Statistics[model][step].parameters
    misfit = ds.auxiliary_data.Statistics[model][step].data.value[0]
    windows = stats["number_misfit_windows"]
    adjsrcs = stats["number_adjoint_sources"]
    event_id = os.path.basename(ds.filename).split(".")[0]
   
    # build nested dictioanaries 
    event_dict = {"misfit": float(misfit), 
                  "windows": int(windows), 
                  "adjsrcs": int(adjsrcs)
                  }
    step_dict = {"{}".format(event_id): event_dict}
    model_dict = {"{}".format(step): step_dict}
    misfit_dict = {model: model_dict}
    
    # To allow multiple cpu's clean, simultaenous write capability, create a 
    # lock file that can only be accessed by one cpu at a time
    fidout_lock = fidout + "_lock"
    while True:
        # another process has control of the misfit file
        if os.path.exists(fidout_lock):
            print("file is locked, waiting")
            time.sleep(5)
        # misfit file is available for writing
        elif os.path.exists(fidout) and not os.path.exists(fidout_lock):
            print("file is available at {}".format(time.asctime()))
            os.rename(fidout, fidout_lock)
            with open(fidout_lock, "r") as f:
                misfit_dict = json.load(f)
                if model in misfit_dict.keys():
                    if step in misfit_dict[model].keys():
                        misfit_dict[model][step][event_id] = event_dict
                    else:
                        misfit_dict[model][step] = step_dict
                else:
                    misfit_dict[model] = model_dict
       
                # Parse misfit dict to give the total misfit
                misfits, windows_all, adjsrcs_all = [], [], []
                for key in misfit_dict[model][step].keys():
                    if key in ["misfit", "windows", "adjsrcs"]:
                        continue
                    misfits.append(misfit_dict[model][step][key]["misfit"])
                    windows_all.append(misfit_dict[model][step][key]["windows"])
                    adjsrcs_all.append(misfit_dict[model][step][key]["adjsrcs"])

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
            print("file has not been written")
            with open(fidout, "w") as f:
                json.dump(misfit_dict, f, indent=4, separators=(',', ':'),
                          sort_keys=True)
            return

         
def write_misfit_stats(ds, model, step_count=0, fidout="./pyatoa.misfits"):
    """
    Deprecated, and not working. See write_misfit_json

    Misfit text file useful for quickly determining misfit information garnered
    from a swatch of pyatoa runs, used by Seisflows

    Caution, Seisflows uses the misfit written from this function to determine
    if misfit is reduced per trial step. The formatter is 8.6f which follows
    the Seisflows convention, but may be too limited for other purposes.

    :type ds: pyasdf.ASDFDataSet
    :param ds: processed dataset, assumed to contain auxiliary_data.Statistics
    :type model: str
    :param model: model number, e.g. "m00"
    :type step: int
    :param step: line search step count
    :type fidout: str
    :param fidout: output file to write the misfit
    """
    import warnings
    warnings.warn("This function is deprecated", DeprecationWarning)

    step = "s{:0>2}".format(step_count)
    stats = ds.auxiliary_data.Statistics[model][step].parameters
    misfit = ds.auxiliary_data.Statistics[model][step].data.value[0]
    windows = stats["number_misfit_windows"]
    adjsrcs = stats["number_adjoint_sources"]
   
    # reformat to write. set format string
    model = int(model[1:])
    event_id = os.path.basename(ds.filename).split(".")[0]
    header = "MODEL\tSTEP\t   EVENT_ID\t\t\t\tMSFT\tWNDW\tADJS\n"
    format_string = "{:>5d}\t{:>4d}\t{:>10s}\t\t{:8.6e}\t{:>4d}\t{:>4d}\n"
   
    # First iteration, first step, write header and first line of data
    if not os.path.exists(fidout):
        with open(fidout, "a") as f:
            f.write(header)
            f.write(format_string.format(
                          model, step_count, event_id, misfit, windows, adjsrcs)
                    )
        return
        
     
    # File has already been written, check for duplicates, overwrite dupes 
    with open(fidout, "rw") as f:
        lines = f.readlines()
        # look at misfits already written, skip header
        for i, line in enumerate(lines[1:]):
            i += 1  # hacky way to keep up with the fact that we start from 1
            _model, _step, _event, _misfit, _window, _adjsrc = line.split()
            # duplicate found, replace misfit, window and adjsrc
            if model == _model and step == _step and  event_id == _event:
                if misfit != _misfit:
                    lines[i][3] = misfit 
                    lines[i][4] = windows
                    lines[i][5] = adjsrcs 
                break
        # LEFT UNFINISHED HERE, MOVED TO JSON WRITING        
 
        # rewrite file from the beginning
        f.seek(0,0)
        f.write(header)
        for line in lines[1:]:
            f.wwrite(format_string.format(
                     line[0], line[1], line[2], line[3], line[4], line[5])
    
                      )

           
def generate_srcrcv_vtk_file(h5_fid, fid_out, model="m00", utm_zone=60,
                             event_fid_out=None):
    """
    It's useful to visualize source receiver locations in Paraview, alongside
    sensitivity kernels. VTK files are produced by Specfem, however they are for
    all receivers, and a source at depth which is sometimes confusing. 
    This function will create source_receiver vtk files using the asdf h5 files,
    with only those receiver that were used in the misfit analysis, and only
    an epicentral source location, such that the source is visible on a top
    down view from Paraview.
    
    Gives the option to create an event vtk file separate to receivers, for
    more flexibility in the visualization.
    
    :type h5_fid: str
    :param h5_fid: path to pyasdf h5 file outputted by pyatoa
    :type fid_out: str
    :param fid_out: output path and filename to save vtk file e.g. 'test.vtk'
    :type model: str
    :param model: h5 is split up by model iteration, e.g. 'm00'
    :type utm_zone: int
    :param utm_zone: the utm zone of the mesh, 60 for NZ
    :type event_fid: str
    :param event_fid: if event vtk file to be made separately
    """
    import pyasdf
    from pyatoa.utils.operations.source_receiver import lonlat_utm

    # get receiver location information in UTM_60 coordinate system from
    # pyasdf auxiliary_data. make sure no repeat stations
    ds = pyasdf.ASDFDataSet(h5_fid)
    sta_x, sta_y, sta_elv, sta_ids = [], [], [], []
    if bool(ds.auxiliary_data):
        for adjsrc in ds.auxiliary_data.AdjointSources[model].list():
            sta = ds.auxiliary_data.AdjointSources[model][adjsrc]
            station_id = sta.parameters["station_id"]
            if station_id in sta_ids:
                continue
            latitude = sta.parameters["latitude"]
            longitude = sta.parameters["longitude"]
            elevation_in_m = sta.parameters["elevation_in_m"]

            x, y = lonlat_utm(lon_or_x=longitude, lat_or_y=latitude,
                              utm_zone=utm_zone, inverse=False)
            sta_x.append(x)
            sta_y.append(y)
            sta_elv.append(elevation_in_m)
            sta_ids.append(station_id)

    # get event location information in UTM_60. set depth at 100 for epicenter
    ev_x, ev_y = lonlat_utm(lon_or_x=ds.events[0].preferred_origin().longitude,
                            lat_or_y=ds.events[0].preferred_origin().latitude,
                            utm_zone=utm_zone, inverse=False
                            )
    ev_elv = 100.0

    # write header for vtk file and then print values for source receivers
    with open(fid_out, "w") as f:
        f.write("# vtk DataFile Version 2.0\n"
                "Source and Receiver VTK file from Pyatoa\n"
                "ASCII\n"
                "DATASET POLYDATA\n"
                )
        # num points equal to number of stations plus 1 event
        f.write("POINTS\t{} float\n".format(len(sta_x)+1))
        f.write("{X:18.6E}{Y:18.6E}{E:18.6E}\n".format(
            X=ev_x, Y=ev_y, E=ev_elv)
        )
        for x, y, e in zip(sta_x, sta_y, sta_elv):
            f.write("{X:18.6E}{Y:18.6E}{E:18.6E}\n".format(X=x, Y=y, E=e))

    # make a separate vtk file for the source
    if event_fid_out:
        with open(event_fid_out, "w") as f:
            f.write("# vtk DataFile Version 2.0\n"
                    "Source and Receiver VTK file from Pyatoa\n"
                    "ASCII\n"
                    "DATASET POLYDATA\n"
                    )
            f.write("POINTS\t1 float\n".format(len(sta_x) + 1))
            f.write("{X:18.6E}{Y:18.6E}{E:18.6E}\n".format(
                X=ev_x, Y=ev_y, E=ev_elv)
            )


def create_stations_adjoint(ds, model, filepath=None):
    """
    TO DO: remove the hardcoded paths for station list
    
    Generate an adjoint stations file for Specfem input by reading in the master
    station list and checking which adjoint sources are available in the
    pyasdf dataset
    
    :type ds: pyasdf.ASDFDataSet
    :param ds: dataset containing AdjointSources auxiliary data
    :type model: str
    :param model: model number, e.g. "m00"
    :type filepath: str
    :param filepath: path/to/specfem/DATA/STATIONS
    :return:
    """
    for f in ['/seis/prj/fwi/bchow/data/STATIONXML/MASTER/'
              'master_station_list.txt',
              '/Users/chowbr/Documents/subduction/data/'
              'STATIONXML/MASTER/master_station_list.txt',
              '/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/primer/'
              'auxiliary_data/stationxml/master_station_list.txt',
              ]:
        if os.path.exists(f):
            master_station_list = f

    stas_with_adjsrcs = []
    for code in ds.auxiliary_data.AdjointSources[model].list():
        stas_with_adjsrcs.append(code.split('_')[1])
    stas_with_adjsrcs = set(stas_with_adjsrcs)

    with open(master_station_list, "r") as f:
        lines = f.readlines()

    event_id = ds.filename.split('/')[-1].split('.')[0]

    # if no output path is specified, save into current working directory with
    # an event_id tag to avoid confusion with other files, else normal naming
    if not filepath:
        pathout = "./STATIONS_ADJOINT_{}".format(event_id)
    else:
        pathout = os.path.join(filepath, "STATIONS_ADJOINT")
    with open(pathout, "w") as f:
        for line in lines:
            if line.split()[0] in stas_with_adjsrcs:
                    f.write(line)


def write_adj_src_to_ascii(ds, model, filepath=None,
                           comp_list=["N", "E", "Z"]):
    """
    take AdjointSource auxiliary data from a pyasdf dataset and write out
    the adjoint sources into ascii files with proper formatting, for input
    into PyASDF
    :param ds:
    :param model_number:
    :param comp_list:
    :return:
    """
    import numpy as np
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
    
    adjsrcs = ds.auxiliary_data.AdjointSources[model]
    event_id = ds.filename.split('/')[-1].split('.')[0]
    
    if not filepath:
        pathcheck = os.path.join("./", event_id)
    else:
        pathcheck = filepath

    if not os.path.exists(pathcheck):
        os.makedirs(pathcheck)
    for adj_src in adjsrcs.list():
        station = adjsrcs[adj_src].path.replace('_', '.')
        fid = "{path}/{sta}.adj".format(path=pathcheck, sta=station)
        with open(fid, "w") as f:
            write_to_ascii(f, adjsrcs[adj_src].data.value)
        # write blank adjoint sources for components with no misfit windows
        if comp_list:
            for comp in comp_list:
                station_check = station[:-1] + comp
                if station_check.replace('.', '_') not in adjsrcs.list():
                    blank_adj_src = adjsrcs[adj_src].data.value
                    blank_adj_src[:, 1] = np.zeros(len(blank_adj_src[:, 1]))
                    blank_template = "{path}/{sta}.adj".format(
                                           path=pathcheck, sta=station_check
                                                               )
                    with open(blank_template, "w") as b:
                        write_to_ascii(b, blank_adj_src)

