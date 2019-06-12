"""
Hacky functions to combine all the outputs figures of pyatoa into a single pdf
"""
import os
import glob
import shutil
import subprocess


def tile_images(image_path, purge=False):
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
    files = glob.glob(os.path.join(image_path, "*_*.png"))
    stanames = []
    for f in files:
        sta = os.path.basename(f).split("_")[1].split(".")[0]
        stanames.append(sta)
    stanames = set(stanames)
    stanames = list(stanames)

    # combine map and waveform figures
    os.chdir(image_path)
    for name in stanames:
        map_name = "map_{}.png".format(name)
        wav_name = "wav_{}.png".format(name)
        if os.path.exists(map_name) and os.path.exists(wav_name):
            subprocess.run(
                ["montage", map_name, wav_name, "-tile", "2x1", 
                 "-geometry", "+0+0", "tile_{}.png".format(name)]
            )
            # remove old images to save space
            if purge:
                os.remove(map_name)
                os.remove(wav_name)

def combine_images(ds, image_path, composite_name="composite.pdf", purge=False):
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
    stations_available = glob.glob(os.path.join(image_path, "tile_*.png"))
    sub_list = []
    for name in sorted_station_names:
        tilename = "tile_{}.png".format(name.split(".")[1])
        # hacky way to check that the tile image exists
        for check in stations_available:
            if tilename in check:
                sub_list.append(tilename)
                break
    
    # run imagemagick convert using the subprocess module
    os.chdir(image_path)
    subprocess.run(["convert"]
                   + sub_list 
                   + [composite_name]
                   )
    
    # remove tile_*.png files 
    if purge:
        for tile in stations_available:
            os.remove(tile)


def tile_and_combine(ds, model, step, figure_path, 
                                      purge_originals=False, purge_tiles=False):
    """
    Tile waveform and map images, combine all images into a single .pdf file
    To be called from within a process script
    """
    event_id = os.path.basename(ds.filename).split(".")[0]
    
    # set up directories and pathnames
    composite_path = os.path.join(figure_path, "composites")
    composite_name = "{e}_{m}_{s}.pdf".format(e=event_id, m=model, s=step) 
    if not os.path.exists(composite_path):
        os.makedirs(composite_path)
     
    # tile and combine images 
    image_path = os.path.join(figure_path, model, event_id)
    tile_images(image_path, purge=purge_originals)
    combine_images(ds, image_path, composite_name, purge=purge_tiles)
    
    # move composite figure to new path
    src = os.path.join(image_path, composite_name)
    dst = os.path.join(composite_path, composite_name)
    os.rename(src, dst)         

