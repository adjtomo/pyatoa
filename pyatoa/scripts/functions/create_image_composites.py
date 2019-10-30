"""
Short script to call tile_and_combine()
"""
import os
import glob
import shutil
import pyasdf

from pyatoa.utils.tools.io import tile_combine_imgs

# set parameters
model_number = "m00"
pyatoa_output = ("/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/"
                 "tests/seisflows_test/hikurangi_trial/pyatoa.output")
purge_originals = False
purge_tiles = True


data_path_template = os.path.join(pyatoa_output, "data", "{}.h5")
figure_path = os.path.join(pyatoa_output, "figures")
image_paths = os.path.join(pyatoa_output, "figures", model, "*")

for image_path in glob.glob(image_paths):
    event_id = os.path.basename(image_path)
    print(event_id)
    ds_fid = data_path_template.format(event_id)
    composite_name = "{}_composite.pdf".format(event_id)
   
    # call tile and combine 
    with pyasdf.ASDFDataSet(ds_fid) as ds:
        tile_and_combine(ds, model, figure_path, purge_originals, purge_tiles)
