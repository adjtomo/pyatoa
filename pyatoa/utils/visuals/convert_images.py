"""
Hacky functions to combine all the outputs figures of pyatoa into a single pdf
Makes use of subprocess and the external program Imagemagick
"""
import os
import glob
import subprocess


def tile_images(wavs_path, maps_path, purge=False):
    """
    Tile the outputted waveform and map images into single pngs placed
    horizontally. Expects the output of the images to be in the form of
    wav_*.png and map_*.png where * represents the station name, i.e. HAZ

    Requires imagemagick

    :type wavs_path: str
    :param wavs_path: where the wav_*.png files are located, outputted by pyatoa
    :type maps_path: str
    :param maps_path: where the map_*.png files are located, outputted by pyatoa
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


def combine_tiles(ds, tiles_path, save_to="./composite.pdf", purge=False):
    """
    Combine the tiled images into a single composite pdf.
    Needs to be run after tile_images. Requires imagemagick

    :type ds: pyasdf.ASDFDataSet()
    :param ds: dataset for the event to get station, event information
    :type tiles_path: str
    :param tiles_path: where the .png files are located, outputted by pyatoa
    :type save_to: str
    :param save_to: where the pdf should be saved and the fid
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
        for i, avail in enumerate(tiles_available):
            if name in avail:
                sorted_tile_list.append(tiles_available[i])            
                break
    
    # run imagemagick convert using the subprocess module
    subprocess.run(["convert"] + sorted_tile_list + [save_to])
    
    # remove tile_*.png files 
    if purge:
        for tile in tiles_available:
            os.remove(tile)


def combine_images(path, globfid="*.png", save_to="./composite.pdf",
                   purge=False):
    """
    Combine a collection of images into a single composite pdf.
    Needs to be run after tile_images. Requires imagemagick

    :type path: str
    :type path: path to search for the images
    :param globfid: wildcard search term to look for images, will be sorted
    :type save_to: str
    :param save_to: path and fid to save the pdf
    :type purge: bool
    :param purge: delete images after tiling to save on space, defeaults to no
    """
    images = glob.glob(os.path.join(path, globfid))
    images.sort()

    # run imagemagick convert using the subprocess module
    subprocess.run(["convert", images, save_to])

    # remove tile_*.png files
    if purge:
        for img in images:
            os.remove(img)


def tile_and_combine(ds, model, figure_path, maps_path, save_to,
                     purge_originals=False, purge_tiles=False):
    """
    Tile waveform and map images, combine all images into a single .pdf file
    To be called from within a process script

    TO DO: take the path naming out of here and move it to the
    """
    event_id = os.path.basename(ds.filename).split(".")[0]
    
    # set up directories and pathnames
    wavs_path = os.path.join(figure_path, model, event_id)
    if not maps_path:
        maps_path = wavs_path
     
    # tile and combine images 
    tile_images(wavs_path=wavs_path, maps_path=maps_path, purge=purge_originals)
    combine_tiles(ds=ds, tiles_path=wavs_path, save_to=save_to,
                  purge=purge_tiles
                  )

