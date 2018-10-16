"""
For generation of input files for Specfem
"""


def create_stations_adjoint(ds, model_number):
    """
    generate an adjoint stations file for specfem input by reading in the master
    station list and checking which adjoint sources are available in the
    pyasdf dataset
    :param ds:
    :param model_number:
    :return:
    """
    master_station_list = ('/Users/chowbr/Documents/subduction/data/'
                           'STATIONXML/MASTER/master_station_list.txt'
                           )
    stas_with_adjsrcs = []
    for code in ds.auxiliary_data.AdjointSources[model_number].list():
        stas_with_adjsrcs.append(code.split('_')[1])
    stas_with_adjsrcs = set(stas_with_adjsrcs)

    with open(master_station_list, "r") as f:
        lines = f.readlines()

    with open("STATIONS_ADJOINT", "w") as f:
        for line in lines:
            if line.split()[0] in stas_with_adjsrcs:
                    f.write(line)