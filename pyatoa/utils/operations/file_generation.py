"""
For generation of input files for Specfem
"""
import os


def sum_residuals(ds, model_number, suffix=None, filepath="./"):
    """
    Residuals are stored in adjoint trace dictionaries, and are needed for
    Seisflows. This will sum the residuals and place them into a specific
    filepath.
    This replaces seisflows.workflow.inversion.write_misfit()

    :param ds:
    :param model_number:
    :param filepath:
    :return:
    """
    adjoint_sources = ds.auxiliary_data.AdjointSources[model_number]
    residuals = []
    for adjsrc in adjoint_sources.list():
        residuals.append(adjoint_sources[adjsrc].parameters["misfit_value"])
    summed_residuals = sum(residuals)

    if suffix:
        import numpy as np
        np.savetxt(fname=os.path.join(filepath, "f_{}".format(suffix)),
                   X=[summed_residuals], fmt="%11.6e"
                   )

    return summed_residuals


def create_stations_adjoint(ds, model_number, filepath=None):
    """
    TO DO: remove the hardcoded paths for station list
    generate an adjoint stations file for specfem input by reading in the master
    station list and checking which adjoint sources are available in the
    pyasdf dataset
    :param ds:
    :param model_number:
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
    for code in ds.auxiliary_data.AdjointSources[model_number].list():
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

