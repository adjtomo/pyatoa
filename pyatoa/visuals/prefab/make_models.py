"""
A script to call the pyatoa.visuals.model.Model class and plot standard looking
slices of a model for easy visualization. Contains a set of configurations
to control the look of the figures.

Note:
    One can easily make a .gif to visualize model slices in psuedo-3d by using
    the imagemagick package, with the command. This could be wrapped in
    subprocess but since it's a one liner, it's pretty easy to do quickly.

    $ convert -delay 20 -loop 1 *png model_animated.gif
"""
import os
from glob import glob
from pyatoa.visuals.model import Model
from pyatoa.visuals.plot_tools import imgs_to_pdf


def call_models(fid, path_out="./figures"):
    """
    Set the parameters in this function, which will then be passed to Model

    Slice attributes:
    :z_slices (list): list of slices to make for top-down figures (normal Z).
        use "surface" for surface projection, otherwise use floats corresponding
        to depth. e.g. slices at surface, 5 and 10km depth would be
        ['surface', 5, 10]
    :x_slices (list): list of slices for side on figures (normal Y-axis)
        these are given in percentage of the axis from 0 to 1. e.g. if you want
        a slice at the middle of the axes, [0.5]
    :y_slices (list): list of slices for side on figures (normal X-axis)

    Parameter attributes:
    :font_factor (float): size of font on the axes objects, default 1.1
    :convert (float): convert the axis by some constant value, e.g. if axes are
                      in units of meters, use 1E-3 to get into units of km
    :zero_origin (bool): set the origin of the X and Y axes to 0
    :cmap (str): colormap to use, some options available,
    :reverse (bool): reverse the colormap
    :title (str): title of the colorbar
    :min_max (list of float): minimum and maximum bounds of the colorbar
    :default_range (bool): overrides min_max, if True, use the default visible
        range for the colorbar
    :num_colors (int): number of discrete colors in colorbar
    :num_labels (int): number of labels on the colorbar
    :round_to (int): only if default_range == True, rounds the bounds of the
        default range to get cleaner min, max and step values in colorbar
    """
    figsize = (1000, 1000)
    x_slices, y_slices, z_slices = None, None, None  # empty initialization

    # File ID's for VTK objects that are to be plotted
    coast_fid = "coast.npy"
    srcs_fid = "srcs.vtk"
    rcvs_fid = "rcvs.vtk"

    # File ID's for output files
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    fid_out = os.path.join(path_out, os.path.basename(fid).split(".")[0])
    save_fid_z = fid_out + "_z_slice_{tag}.png"
    save_fid_x = fid_out + "_x_slice_{tag}.png"
    save_fid_y = fid_out + "_y_slice_{tag}.png"

    # Pick specific slices to create
    z_slices = ["surface", 5, 10, 15, 25, 30, 35, 40, 45, 50]
    x_slices = [0.25, 0.5, 0.75]
    y_slices = [0.25, 0.5, 0.75]

    # Set the parameters for all figures generated, passed as kwargs
    par = {"font_factor": 1.1,
           "convert": 1E-3,
           "zero_origin": True,
           "cmap": "RdYlBu",
           "reverse": False,
           "title": "Vs log(m/m00)",
           "min_max": None,
           "default_range": False,
           "num_colors": 51,
           "num_clabels": 5,
           "round_to": None,
           "show": False
           }

    # Determine title and parameters for different plot types
    if "log" in fid_out:
        tit = os.path.basename(fid).split(".")[0].split("_")[-1]
        num = os.path.basename(fid).split(".")[0].split("_")[-2]
        if num == "init":
            num = 0
        par["title"] = f"{tit.capitalize()} log(m{int(num):0>2}/m00)"
        par["cmap"] = "RdYlBu"
        if "vs" in fid_out:
            par["default_range"] = False
        elif "vs" in fid_out:
            par["default_range"] = True
    elif "poisson" in fid_out:
        par["cmap"] = "viridis" 
        par["title"] = "Poissons Ratio"
        par["default_range"] = False
        par["min_max"] = None

    # Call the model maker
    model = Model(fid=fid, srcs=srcs_fid, rcvs=rcvs_fid, coast=coast_fid,
                  figsize=figsize, zero_origin=par["zero_origin"],
                  convert=par["convert"], offscreen=True, logging=False)

    # Make the models
    if z_slices:
        par["save"] = save_fid_z
        z_fids = []
        for z in z_slices:
            model.startup()
            fid_out = model.plot_model_topdown(depth_km=z, **par)
            z_fids.append(fid_out)
            imgs_to_pdf(z_fids, 
                        save_fid_z.format(tag="all").replace("png", "pdf"))
    if x_slices:
        par["save"] = save_fid_x
        x_fids = []
        for x in x_slices:
            model.startup()
            fid_out = model.plot_depth_cross_section(choice="X",
                                                     slice_at=x, **par)
            x_fids.append(fid_out)
            imgs_to_pdf(x_fids,
                        save_fid_x.format(tag="all").replace("png", "pdf"))
    if y_slices:
        par["save"] = save_fid_y
        y_fids = []
        for y in y_slices:
            model.startup()
            fid_out = model.plot_depth_cross_section(choice="Y",
                                                     slice_at=y, **par)
            y_fids.append(fid_out)
            imgs_to_pdf(y_fids,
                        save_fid_y.format(tag="all").replace("png", "pdf"))

    # Remove all the pngs as theyre in the pdf now
    for pngfids in [x_fids, y_fids, z_fids]:
        for fid_ in pngfids:
            os.remove(fid_)


if __name__ == "__main__":
    fids = glob("*0008*vp*.vtk")
    for fid in fids:
        call_models(fid)

