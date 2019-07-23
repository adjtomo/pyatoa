"""
Create a config file in .json format

As with all packages, a configuration file provides a simple, single-entry
point for the User with the package. Here we provide a python script that will
create a .json file that will be parsed by process.py, containing
parameters that control the expectations and outputs of Pyatoa
"""
import json


def sfconfig(min_period=10, max_period=30, filter_corners=4,
             rotate_to_rtz=False, unit_output="DISP", pyflex_config="default",
             adj_src_type="multitaper_misfit", paths_to_waveforms=[],
             paths_to_responses=[], set_logging="info", fix_windows=False,
             plot_waveforms=True, plot_maps=True, synthetics_only=False,
             window_amplitude_ratio=.0, write_misfit_json=True, 
             tile_and_combine=False,  purge_originals=False, purge_tiles=True,
             create_srcrcv_vtk=True, snapshot=True, misfits_json="misfits.json",
             fidout='./sfconfig.json'):
    """
    user set dictionary parameters
    :type min_period: int
    :param min_period: minimum filter period in seconds
    :type max_period: int
    :param max_period: maximum filter period in seconds
    :type filter_corners: int
    :param filter_corners: number of corners to use on the waveform filter
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
    :type synthetics_only: bool
    :param synthetics_only: synthetic synthetic example, will search for
        observation data based on synthetic search routines
    :type window_amplitude_ratio: float
    :param window_amplitude_ratio: global amplitude ratio to exclude misfit
        windows. Evaluates as peak_window/peak_waveform > window_amplitude_ratio
    :type write_misfit_json: bool
    :param write_misfit_json: write misfit information, number of windows,
        number of adjoint sources to a json file with name 'misfit_json_file'
    :type tile_and_combine: bool
    :param tile_and_combine: combine the waveforms and maps into a single pdf
        file, saved into 'figures'/composites
    :type purge_*: bool
    :param purge_originals: purge the wav_*.png and map_*.png from Pyatoa, only
        if 'tile_and_combine'==True
    :param purge_tiles: purge intermediate tile files which combine wav and map
        only if 'tile_and_combine'==True
    :type create_srcrcv_vtk: bool
    :param create_srcrcv_vtk: If finalize should create vtk files for src n rcvs
    :type snapshot: bool
    :param snapshot: save a copy of the HDF5 files from each iteration, incase
        e.g. file corruption, allows for temporary backups
    :type misfits_json: str
    :param misfits_json: name of the misfits json file, only if
        write_misfits_json==True, default='misfits.json'
    :type fidout: str
    :param fidout: filename to save the .json config file
    """
    pd = {
        # Pyatoa Config object
        "min_period": min_period,
        "max_period": max_period,
        "filter_corners": filter_corners,
        "rotate_to_rtz": rotate_to_rtz,
        "unit_output": unit_output,
        "pyflex_config": pyflex_config,
        "adj_src_type": adj_src_type,
        "paths_to_waveforms": paths_to_waveforms,
        "paths_to_responses": paths_to_responses,

        # Pyatoa Considerations
        "set_logging": set_logging,
        "fix_windows": fix_windows,
        "plot_waveforms": plot_waveforms,
        "plot_maps": plot_maps,
        "synthetics_only": synthetics_only,
        "window_amplitude_ratio": window_amplitude_ratio,

        # Pyatoa Outputs
        "write_misfit_json": write_misfit_json,
        "tile_and_combine": tile_and_combine,
        "purge_originals": purge_originals,
        "purge_tiles": purge_tiles,
        "create_srcrcv_vtk": create_srcrcv_vtk,
        "snapshot": snapshot,

        # Pyatoa Paths (output directory set in Seisflows)
        "misfits_json": misfits_json,
    }

    with open(fidout, "w") as f:
        json.dump(pd, f, indent=4, separators=(',', ':'))     

    return pd


if __name__ == "__main__":
    sfconfig()
