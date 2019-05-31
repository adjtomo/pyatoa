"""
hacky way to combine all the outputs figures of pyatoa into a single pdf
"""
import os
import glob
import shutil
import pyasdf
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
    for n in stanames:
        map_name = "map_{}.png".format(n)
        wav_name = "wav_{}.png".format(n)
        if os.path.exists(map_name) and os.path.exists(wav_name):
            subprocess.run(
                ["montage", map_name, wav_name, "-tile", "2x1", 
                 "-geometry", "+0+0", "tile_{}.png".format(n)]
            )
        # remove old images to save space
        if purge:
            for tag in ["map", "wav"]:
                os.remove("{}_{}.png".format(tag, n))


def combine_images(ds, image_path, purge=False):
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

    event_id = os.path.basename(ds.filename).split(".")[0]
    os.chdir(image_path)
    subprocess.run(["convert"] + sub_list + ["{}_composite.pdf".format(event_id)])
    if purge:
        for tile in stations_available:
            os.remove(tile)



if __name__ == "__main__":
    model_number = "m00"
    pyatoa_output = ("/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/"
                     "tests/seisflows_test/hikurangi_trial/pyatoa.output")
    image_path = os.path.join(pyatoa_output, "figures", model_number, "*")
    data_path = os.path.join(pyatoa_output, "data", "{}.h5")
    
    for images in glob.glob(image_path):
        event_id = os.path.basename(images)
        print(event_id)
        ds_fid = data_path.format(event_id)
        ds_temp = shutil.copy(src=ds_fid, dst=ds_fid + "_temp")
        ds = pyasdf.ASDFDataSet(ds_temp)
        tile_images(image_path=images)
        combine_images(ds, images, purge=True)
        os.remove(ds_temp)
