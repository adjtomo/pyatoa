"""
Hacky functions to combine all the outputs figures of pyatoa into a single pdf
"""
import os
import glob
import subprocess


def tile_images(wavs_path, maps_path, purge=False):
    """
    Tile the outputted waveform and map images into single pngs placed
    horizontally. Expects the output of the images to be in the form of
    wav_*.png and map_*.png where * represents the station name, i.e. HAZ
    :type image_path: str
    :param image_path: where the .png files are located, outputted by pyatoa
    :type purge: bool
    :param purge: delete images after tiling to save on space, defeaults to no
    """
    # get set of station names
    files = glob.glob(os.path.join(wavs_path, "*_*.png"))
    stanames = []
    for f in files:
        sta = os.path.basename(f).split("_")[1].split(".")[0]
        stanames.append(sta)
    stanames = set(stanames)
    stanames = list(stanames)

    # combine map and waveform figures
    for name in stanames:
        map_name = os.path.join(maps_path, f"map_{name}.png")
        wav_name = os.path.join(wavs_path, f"wav_{name}.png")
        tile_name = os.path.join(wavs_path, f"tile_{name}.png")
        if os.path.exists(map_name) and os.path.exists(wav_name):
            subprocess.run(
                ["montage", map_name, wav_name, "-tile", "2x1", 
                 "-geometry", "+0+0", tile_name]
            )
            # remove old images to save space
            if purge:
                os.remove(map_name)
                os.remove(wav_name)

        # purge the maps which have no wav counterparts
        # this was necessary when we were making maps each step, but it's been
        # since changed to only make maps on the first go
        # elif os.path.exists(map_name) and purge:
        #     os.remove(map_name)


def combine_images(ds, tiles_path, save_to="./composite.pdf", purge=False):
    """
    Combine the tiled images into a single composite pdf.
    Needs to be run after tile_images
    :type ds: pyasdf.ASDFDataSet()
    :param ds: dataset for the event to get station, event information
    :type image_path: str
    :param image_path: where the .png files are located, outputted by pyatoa
    :type purge: bool
    :param purge: delete images after tiling to save on space, defeaults to no
    """
    from pyatoa.utils.operations.source_receiver import sort_by_backazimuth
   
    # sort stations by backazimuth for easier visualization 
    sorted_station_names = sort_by_backazimuth(ds)

    tiles_available = glob.glob(os.path.join(tiles_path, "tile_*.png"))
    tiles_available = [os.path.basename(_) for _ in tiles_available]    
 
    # Loop through sorted station names, if that name is available as a tile
    # add to the list that will be sent to subprocess
    sorted_tile_list = []
    for name in sorted_station_names:
        for i, avail in eumerate(tiles_available):
            if name in avail:
                sorted_tile_list.append(tiles_available[i])            
                break
    
    # run imagemagick convert using the subprocess module
    subprocess.run(["convert"] + sub_list  + [save_to])
    
    # remove tile_*.png files 
    if purge:
        for tile in tiles_available:
            os.remove(tile)


def tile_and_combine(ds, model, step, figure_path, maps_path,
                                      purge_originals=False, purge_tiles=False):
    """
    Tile waveform and map images, combine all images into a single .pdf file
    To be called from within a process script
    """
    event_id = os.path.basename(ds.filename).split(".")[0]
    
    # set up directories and pathnames
    composite_path = os.path.join(figure_path, "composites")
    wavs_path = os.path.join(figure_path, model, event_id)
    if not maps_path:
        maps_path = wavs_path
    if not os.path.exists(composite_path):
        os.makedirs(composite_path)

    # Create the name of the pdf to save to
    save_to = os.path.join(composite_path,
                           "{e}_{m}_{s}.pdf".format(e=event_id, m=model, s=step) 
                           )
     
    # tile and combine images 
    tile_images(wavs_path=wavs_path, maps_path=maps_path, purge=purge_originals)
    combine_images(ds=ds, tiles_path=wavs_path, save_to=save_to, 
                   purge=purge_tiles
                   )

