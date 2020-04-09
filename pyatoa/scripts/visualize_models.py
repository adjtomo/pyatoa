"""
An automation script to take the output models from Seisflows, convert the .bin
files into .vtk files, collect into approriate directories, create differences
and plot standard visualizations using Mayavi
"""
import os
import sys
import time
import numpy as np
import subprocess
from glob import glob
from pyatoa.utils.read import read_specfem_vtk

# Avoids 'QXcbConnection: Could not connect to display' error
# os.environ['QT_QPA_PLATFORM']='offscreen'
from pyatoa.visuals.model import Model


# Globally accessible paths and parameters
cwd = os.getcwd()
path_to_specfem = os.path.join(cwd, "scratch", "solver", "mainsolver")
path_to_vtks = os.path.join(cwd, "pyatoa.io", "figures", "vtks")
path_to_xcomb = "/home/chowbr/primer/run_xcombine_vol_data_vtk.sh"

# Paths to auxiliary data
coast_fid = "/home/chowbr/primer/auxiliary/coastline/coast.npy"
srcs_fid = os.path.join(path_to_vtks, "srcs.vtk")
rcvs_fid = os.path.join(path_to_vtks, "rcvs.vtk")

# Parameter choices
parameters = ["vs", "vp"]
z_slices = ["surface", 5, 10, 15, 25, 30, 35, 40, 45, 50]
x_slices = [0.25, 0.5, 0.75]
y_slices = [0.25, 0.5, 0.75]

# Set the parameters for all figures generated, passed as kwargs
figsize = (1000, 1000)
model_par = {"font_factor": 1.1,
             "convert": 1E-3,
             "zero_origin": True,
             "cmap": "RdYlBu",
             "reverse": False,
             "title": "{} log(m/m00)",
             "min_max": [-0.25, 0.25],
             "default_range": False,
             "num_colors": 51,
             "num_clabels": 5,
             "round_to": None,
             "show": False
             }


def clean_dir():
    """
    Clean directory before progressing
    """
    fids = glob(os.path.join(path_to_specfem, "SUM", "*"))
    if fids:
        check = input("files found in SUM directory, clean? [y/(n)]: ")
        if check == "y":
            for f in fids:
                os.remove(f)
        else:
            sys.exit(0)


def pre_organize(file_tag="model"):
    """
    Bookkeeping function that symlinks all files with the unique tag names
    so that only one call to xcombine_vol_data_vtk is required

    :type file_tag: str
    :param file_tag: tag from Seisflows, ['model', 'gradient', 'kernel']
    :rtype: int
    :return: number of unique models that need to be turned into .vtk files
    """
    # Number of counts relates to time required for xcombine_vol_data_vtk
    count = 0

    directories = glob(os.path.join("output", f"{file_tag}_????"))
    for directory in directories:
        unique_id = os.path.basename(directory)
        for par in parameters:
            count += 1

            srcs = glob(os.path.join(directory, f"proc*_{par}.bin"))
            for src in srcs:
                # e.g. proc000000_vs.bin
                filename = os.path.basename(src)
                # Add a unique tag to the filename
                parts = filename.split("_")
                new_filename = "_".join([parts[0], unique_id, parts[1]])
                # Set the source and destination for symlink
                dst = os.path.join(path_to_specfem, "SUM", new_filename)

                os.symlink(os.path.abspath(src), dst)

    print(f"{count} unique ids symlinked")
    return count


def run_xcombine_vol_data_vtk(count):
    """
    Use subprocess to call sbatch of an existing script. A bit hacky but avoids
    having to generate an entire sbatch script in here.

    :type count: int
    :param count: number of unique models that should be generated
    :return:
    """
    os.chdir(path_to_specfem)
    subprocess.call(f"sbatch {path_to_xcomb}", shell=True)
    # Wait for all files to be created by sbatch
    while True:
        files = glob(os.path.join(path_to_specfem, "SUM", "*.vtk"))
        if len(files) == count:
            break
        time.sleep(10)


def post_organize():
    """
    Remove symlinks and organize files in tagged directories, ensure that the
    directory containing symlinks is empty after organization
    """
    # Move files
    for src in glob(os.path.join(path_to_specfem, "SUM", "*.vtk")):
        # e.g. model, kernel, gradient
        tag = os.path.basename(src).split("_")[0]
        dst_dir = os.path.join(path_to_vtks, tag)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        dst = os.path.join(dst_dir, os.path.basename(src))

        os.rename(src, dst)

    # Remove symlinks
    for fid in glob(os.path.join(path_to_specfem, "SUM", "*.bin")):
        if os.path.islink(fid):
            os.remove(fid)

    # Just make sure the directory is empty, otherwise something went wrong
    assert(len(glob(os.path.join(path_to_specfem, "SUM", "*"))) == 0)


def difference_models(log=True, poissons=True, outdir="model"):
    """
    Run diff_vtk() on all models, only if files don't already exist.
    Need to construct filenames before checking that they exist:
    Create either net model updates with log differences, or poissons ratios
    with poissons

    :type log: bool
    :param log: create net model update using log differences
    :type poissons: bool
    :param poissons: create poissons ratio vtk files
    :type outdir: str
    :param outdir: if given, name of the directory to save to the differenced 
        models to. If not given, saved into the model directory
    """
    outdir = os.path.join(path_to_vtks, outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create net model update, or log differences
    if log:
        for par in parameters:
            model_init = os.path.join(
                path_to_vtks, "model", f"model_init_{par}.vtk")
            models = glob(os.path.join(path_to_vtks, "model", "model*"))
            for model in models:
                if model == model_init:
                    continue
                fid_out = os.path.join(outdir, f"log_{os.path.basename(model)}")
                if not os.path.exists(fid_out):
                    diff_vtk(model_a=model, model_b=model_init, fidout=fid_out,
                             method="log")

    # Calculate Poissons ratio
    if poissons:
        models_vp = glob(os.path.join(path_to_vtks, "model", "model*vp.vtk"))
        for model_vp in models_vp:
            model_vs = model_vp.replace("vp", "vs")
            fid_out = model_vp.replace("vp", "poissons")
            fid_out = os.path.join(outdir, os.path.basename(fid_out))
            if not os.path.exists(fid_out):
                diff_vtk(model_a=model_vp, model_b=model_vs, fidout=fid_out,
                         method="poissons")


def make_models(tags=["models"]):
    """
    Call vis_models for all the models created
    """
    for t in tags:
        for fid in glob(os.path.join(path_to_vtks, t, "*.vtk")):
            if "init" in fid:
                continue
            vis_model(fid)


def vis_model(fid):
    """
    Call the Pyatoa Model class to generate standardized model plots for all
    available models. Parameters are set internally in this function
    :return:
    """
    # File ID's for output files
    figure_dir = os.path.join(path_to_vtks, "figures")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    fid_out = os.path.join(figure_dir, os.path.basename(fid).split(".")[0])

    save_fid_z = fid_out + "_z_{tag}.png"
    save_fid_x = fid_out + "_x_{tag}.png"
    save_fid_y = fid_out + "_y_{tag}.png"

    # Determine title for colorbar
    if "log" in fid_out:
        par = os.path.basename(fid).split(".")[0].split("_")[-1]
        num = os.path.basename(fid).split(".")[0].split("_")[-2]
        model_par["title"] = f"{par.capitalize()} log(m{int(num):0>2}/m00)"
    elif "poisson" in fid_out:
        model_par["title"] = "Poissons Ratio"

    # Call the model maker
    model = Model(fid=fid, srcs=srcs_fid, rcvs=rcvs_fid, coast=coast_fid,
                  figsize=figsize, zero_origin=model_par["zero_origin"],
                  convert=model_par["convert"])

    # Make the models
    if z_slices:
        model_par["save"] = save_fid_z
        for z in z_slices:
            model.startup()
            model.plot_model_topdown(depth_km=z, **model_par)
    if x_slices:
        model_par["save"] = save_fid_x
        for x in x_slices:
            model.startup()
            model.plot_depth_cross_section(choice="X", slice_at=x, **model_par)
    if y_slices:
        model_par["save"] = save_fid_y
        for y in y_slices:
            model.startup()
            model.plot_depth_cross_section(choice="Y", slice_at=y, **model_par)


def diff_vtk(model_a, model_b, fidout, method="subtract"):
    """
    Read each model and scan line by line, difference all necessary values,
    write out a new vtk file that is the combination of A and B

    Different methods are:
    1) subtract: c = a - b
    2) log: c = log(a / b) ~ (a-b)/b, where a=m_i, b=m_0
    3) poisson: c = 0.5 * (a**2 - 2 * b**2) / (a**2 - b**2)
        where a=Vp, b=Vs

    :type model_a: str
    :param model_a: path to model file A
    :type model_b: str
    :param model_b: path to model file B
    :type method: str
    :param method: method for differencing,
        choice = ['subtract', 'log', 'poissons']
    :type fidout: str
    :param fidout: name of the output file to be written
    :rtype differences: list
    :return differences: list of the differences in values between a and b
    """
    # read files
    model_a, header_dict_a = read_specfem_vtk(model_a)
    model_b, header_dict_b = read_specfem_vtk(model_b)

    # check that the files have the same characteristics before parsing
    for key in header_dict_a.keys():
        if key == "scalars":
            continue
        elif header_dict_a[key] != header_dict_b[key]:
            sys.exit("{} not equal".format(key))

    # parse through models together and separate by len of line, skip header
    differences = []
    start = header_dict_a["data_line"]
    for a, b in zip(model_a[start:-1], model_b[start:-1]):
        try:
            a = float(a.strip())
            b = float(b.strip())
            if method == "subtract":
                difference = a - b
            # Take the natural log of the the quotient of a and b, this gives
            # to first order approximation, the percent difference. Yoshi said
            # Albert Tarantola said, "always view models in log space"
            elif method == "log":
                difference = np.log(a / b)
            elif method == "poissons":
                difference = 0.5 * (a**2 - 2 * b**2) / (a**2 - b**2)

            differences.append(difference)
        except ValueError:
            print("value error")

    # Write out the differences to a new file, newline at the end to play nice
    with open(fidout, "w") as f:
        f.writelines(model_a[:header_dict_a["data_line"]])
        for diff in differences:
            if diff == 0:
                f.write("{:13.5f}    \n".format(float(diff)))
            elif abs(diff) > 1:
                f.write("{:13.5f}    \n".format(float(diff)))
            else:
                f.write("{:13.10f}    \n".format(float(diff)))
        f.write("\n")

    return differences


if __name__ == "__main__":
    # Run functions in order
    file_tags = ["model"]
    # c = 0
    # clean_dir()
    # for ftag in file_tags:
    #     c_ = pre_organize(ftag)
    #     c += c_
    # run_xcombine_vol_data_vtk(c)
    # post_organize()
    # difference_models(outdir="diffs")
    make_models(tags=["diffs"])


