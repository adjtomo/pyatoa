"""
As with all packages, a configuration file provides a simple, single-entry
point for the User with the package. Here we provide a python script that will
create a .json file that will be parsed by process.py, containing
parameters that control the expectations and outputs of Pyatoa
"""


def sfconfig():
    """
    user set dictionary parameters
    :type min_period: int
    :param min_period: minimum filter period in seconds
    :type max_period: int
    :param max_period: maximum filter period in seconds
    :type rotate_to_rtz: bool
    :param rotate_to_rtz: rotate compoments from north east to radial transverse
    :type unit_output: str
    :param unit_output: ground motion units, can be 'DISP', 'VEL' or 'ACC'
    :type pyflex_config: str
    :param pyflex_config: user-set parameters to be given to pyflex, for
        avilable config types, see:
        ~pyatoa.utils.configurations.external.pyflex_configs()~
    :type adj_src_type: str
    :param adj_src_type: adjoint source calculation method for Pyadjoint,
        available: 'waveform', 'cc_traveltime_misfit', 'multitaper_misfit'
    :type paths_to_*: list of str
    :param paths_to_synthetics: paths to folder containing synthetic data
        Pyatoa expects .sem files inside, e.g. paths_to_synthetics/*semd
    :param paths_to_waveforms: paths to folder containg observation data
        Pyatoa expects the resulting folder structure to follow SEED convention
    :param paths_to_responses: paths to STATIONXML files for observations
        Pyatoa expects the resulting folder structure to follow SEED convention
    :type set_logging: bool
    :param set_logging: turn on logging for detailed output information
    :type fix_windows: bool
    :param fix_windows: if misfit windows already exist in the ASDFDataSet
        do not calculate new windows, instead recalculate adjoint sources inside
        the already given windows
    :type plot_waveforms: bool
    :param plot_waveforms: plot 3 component waveform data with misfit windows
        and adjoint sources, save them into 'figure_directory'
    :type plot_maps: bool
    :param plot_maps: plot source-receiver map and save into 'figure_directory'
    :type write_misfits_json: bool
    :param write_misfits_json: write misfit information, number of windows,
        number of adjoint sources to a json file with name 'misfit_json_file'
    :type tile_and_combine: bool
    :param tile_and_combine: combine the waveforms and maps into a single pdf
        file, saved into 'figures'/composites
    :type purge_*: bool
    :param purge_originals: purge the wav_*.png and map_*.png from Pyatoa, only
        if 'tile_and_combine'==True
    :param purge_tiles: purge intermediate tile files which combine wav and map
        only if 'tile_and_combine'==True
    :type *_diry: str
    :param figure_dir: name of the figure directory, default='figures'
    :param data_dir: name of the data directory, default='data'
    :param misfit_dir: name of scratch misfit directory, default='misfits'
    :type misfit_json: str
    :param misfit_json: name of the misfits json file, only if
        write_misfits_json==True, default='misfits.json'

    """
    pd = {
        # Pyatoa Config object
        "min_period": 10,
        "max_period": 30,
        "filter_corners": 4,
        "rotate_to_rtz": False,
        "unit_output": "DISP",
        "pyflex_config": "default",
        "adj_src_type": "multitaper_misfit",
        "paths_to_waveforms": [],
        "paths_to_responses": [],

        # Pyatoa Considerations
        "set_logging": True,
        "fix_windows": False,
        "plot_waveforms": True,
        "plot_maps": True,

        # Pyatoa Outputs
        "write_misfit_json": True,
        "tile_and_combine": True,
        "purge_originals": True,
        "purge_tiles": True,
        "create_srcrcv_vtk": True,
        "create_src_vtk": False,

        # Pyatoa Paths (output directory set in Seisflows)
        "figure_dir": "figures",
        "data_dir": "data",
        "misfit_dir": "misfits",
        "vtk_dir": "vtks",
        "misfits_json": "misfits.json",
    }
    return pd


