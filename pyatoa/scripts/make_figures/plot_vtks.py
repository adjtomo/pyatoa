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
from pyatoa import VTKModeler
from pyatoa.utils.images import imgs_to_pdf


def save_pdf(fids, fid_out, clean=True):
    """
    Convenience function to collect pngs into a pdf and then remove pngs

    :type fids: list
    :param fids: list of .png files to convert to .pdf
    :type fid_out: str
    :param fid_out: name of output .pdf file
    :type clean: bool
    :pararm clean: delete the original .png files
    """
    imgs_to_pdf(fids, fid_out)
    if clean:
        for fid in fids:
            os.remove(fid)


if __name__ == "__main__":
    # Set file ids here
    vtk_fid = "vs_chkbd.vtk"
    src_fid = "srcs.vtk"
    rcv_fid = "rcvs.vtk"
    coast_fid = "coast.npy"
    path_out = "./"

    # Get the file tag based on the file id of the VTK file
    tag = os.path.splitext(os.path.basename(vtk_fid))[0]

    # Initiate the class with some preset keyword arguments
    vm = VTKModeler(cmap="RdYlBu", reverse=True, num_clabels=3, num_colors=30,
                    round_to=1, scale_axes=1E-3, zero_origin=True,
                    xlabel="E [km]", ylabel="N [km]", zlabel="Z [km]",
                    figsize=(1000, 1000), offscreen=True,
                    )
    # Load in the VTK files
    vm.load(fid=vtk_fid, src_fid=src_fid, rcv_fid=rcv_fid, coast_fid=coast_fid)

    # Create depth slices and combine into a single PDF
    save_fids = []
    for depth_km in ["surface", 5, 10, 15, 20, 25, 50]:
        save_fid = os.path.join(path_out, f"{tag}_Z_{depth_km}.png")
        save_fids.append(save_fid)

        vm.depth_slice(depth_km=depth_km, show=False, save=save_fid)

    save_pdf(save_fids, os.path.join(path_out, f"{tag}_Z.pdf"))

    # Create cross-sections for X and Y axes, combine into separate PDFs
    for axis in ["X", "Y"]:
        save_fids = []
        for pct in [0.25, 0.5, 0.75]:
            save_fid = os.path.join(path_out, f"{tag}_{axis}_{pct*1E2}pct.png")
            save_fids.append(save_fid)
            vm.cross_section(axis=axis, pct=pct, show=False, save=save_fid)

        save_pdf(save_fids, os.path.join(path_out, f"{tag}_{axis}.pdf"))






