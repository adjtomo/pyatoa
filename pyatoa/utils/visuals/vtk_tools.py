"""
Tools to use the Python package mayavi to interact with .vtk files.
These were created because I was tired of spending all my time fine tuning
plots in Paraview, only to remake them all over again. Hopefull this set of
tools allows for quick, standard looking plots of vtk files for easy
visualization of kernels, models, etc.
"""
import os
import numpy as np
from numpy import array
from mayavi import mlab
from mayavi.modules.surface import Surface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane
from mayavi.modules.slice_unstructured_grid import SliceUnstructuredGrid

from ipdb import set_trace
from IPython import embed


def startup(fid):
    """
    Open a VTK file and return the engine that is visualizing it

    :type fid: str
    :param fid: file id of the .vtk file
    :rtype fig: mlab.figure
    :return fig: figure object from mlab
    :rtype engine: mayavi.core.engine.Engine
    :return engine: engine that is visualizing the vtk file'
    :rtype vtk_file_reader: mayavi.sources.vtk_file_reader.VTKFileReader
    :return vtk_file_reader: the opened vtk file
    """
    # Instantiate mlab
    fig = mlab.figure()
    engine = mlab.get_engine()
    engine.scenes[0].scene.background = (1.0, 1.0, 1.0)  # background to white
    vtk_file_reader = engine.open(fid)

    return fig, engine, vtk_file_reader


def set_colorscale(cmap='RdYlBu', reverse=False, min_max=None,
                   colorbar=True, title=None, nb_labels=None, num_col=20):
    """
    Set the colorscale for the plot and also create a colorbar

    :type cmap: str
    :param cmap: colorscale to use, matches names from ParaView
    :type reverse: bool
    :param reverse: reverse the colorscale
    :type min_max: list
    :param min_max: min and max values for color scale as [min, max]
    :type colorbar: bool
    :param colorbar: create and show a colorbar
    :type title: str
    :param title: colorbar title
    :type nb_labels: int
    :param nb_labels: number of labels to put on the colorbar
    """
    # Set the colorscale
    cbar = mlab.colorbar(title=title, orientation="horizontal",
                         label_fmt="%-#.2f", nb_labels=nb_labels)
    cbar.lut_mode = cmap
    cbar.reverse_lut = reverse
    cbar.number_of_colors = num_col
    cbar.number_of_labels = num_col // 2 + 1

    # Bound the colorscale
    cbar.use_default_range = False
    if not min_max:
        val = np.ceil(max(abs(cbar.data_range)))
        min_max = [-val, val]
    cbar.data_range = array(min_max)

    # Create colorbar
    cbar.show_scalar_bar = colorbar
    cbar.show_legend = colorbar

    if colorbar:
        # These are for the values
        cbar.label_text_property.bold = False
        cbar.label_text_property.italic = False
        cbar.label_text_property.orientation = 0.  # rotate labels
        cbar.label_text_property.color = (0.0, 0.0, 0.0)  # black font
        cbar.label_text_property.font_size = 10

        # These are for the title
        cbar.title_text_property.italic = False
        cbar.title_text_property.bold = False
        cbar.title_text_property.color = (0.0, 0.0, 0.0)  # black font
        cbar.title_text_property.font_size = 5

        # Create some uniform look for colorbar
        cbar.scalar_bar.orientation = 'horizontal'  # horizontal scalebar
        cbar.scalar_bar.text_position = 'precede_scalar_bar'  # title below
        cbar.scalar_bar.title_ratio = 0.36  # smaller title size
        cbar.scalar_bar.bar_ratio = 0.325  # thickness of colorbar

        # Makes the colorbar interactive and changes its size
        cbar.scalar_bar_representation.proportional_resize = True
        cbar.scalar_bar_representation.position = array([0.33, .0005])  # location
        cbar.scalar_bar_representation.position2 = array([0.33, 0.062])  # size
    return cbar


def set_axes(xlabel="E (m)", ylabel="N (m)", zlabel="Z (m)",
             xyz=[True, True, True]):
    """
    Create an axis object to wrap around data
    :return:
    """
    axes = mlab.axes(color=(0., 0., 0.,), line_width=3., xlabel=xlabel,
                     ylabel=ylabel, zlabel=zlabel, x_axis_visibility=xyz[0],
                     y_axis_visibility=xyz[1], z_axis_visibility=xyz[2])

    # Edit full properties
    axes.axes.label_format = '%-#3.2E'
    axes.axes.number_of_labels = 5
    axes.axes.font_factor = .75

    # Edit label attributes
    axes.label_text_property.color = (0.0, 0.0, 0.0)  # black font
    axes.label_text_property.bold = False
    axes.label_text_property.italic = False
    axes.label_text_property.orientation = 0.0
    # axes.label_text_property.font_size = 14  # doesnt work

    # Edit title attributes
    axes.title_text_property.color = (0.0, 0.0, 0.0)  # black font
    axes.title_text_property.italic = False
    axes.title_text_property.bold = False
    axes.title_text_property.orientation = 0.0
    # axes.title_text_property.font_size = 14  # doesn't work

    return axes


def coastline(coast_fid, z_value=1000, color=(0., 0., 0.,), ):
    """
    Plot coastline on top of plot. Coastline should be an npy file that is an
    N x 3 array with the columns representing x, y, z

    :type coast_fid: str
    :param coast_fid: fid for npy file
    :type z_value: float
    :param z_value: height of the coastline in meters, negative down
    :type color: tuple of floats
    :param color: tuple representation of color, defaults to black
    """
    coast = np.load(coast_fid)
    x = coast[:, 0]
    y = coast[:, 1]
    z = np.ones(len(x)) * z_value

    p3d = mlab.points3d(x, y, z, color=color, mode="2dcircle",
                        scale_factor=800.)

    return p3d


def show_surface(engine, vtk_file_reader, view="top-down"):
    """
    Show the surface of the vtk file from the top down

    :param engine:
    :param vtk_file_reader:
    :return:
    """
    surface = Surface()
    engine.add_filter(surface, vtk_file_reader)


def cut_plane(engine, vtk_file_reader, depth_m):
    """
    Slice the data at a certain depth, format is depth is positive

    :type engine: mayavi.core.engine.Engine
    :param engine: engine that is visualizing the vtk file'
    :type vtk_file_reader: mayavi.sources.vtk_file_reader.VTKFileReader
    :param vtk_file_reader: the opened vtk file
    :type depth_m: float
    :param depth_m: depth value in m
    """
    cut = ScalarCutPlane()
    engine.add_filter(cut, vtk_file_reader)

    # Set some attributes of the cut plane
    cut.implicit_plane.origin[-1] = depth_m
    cut.implicit_plane.normal = array([0, 0, 1])
    cut.implicit_plane.widget.enabled = False


    # cut.implicit_plane.widget.normal_to_z_axis = True
    # cut.implicit_plane.widget.origin[-1] = depth_m
    # cut.implicit_plane.plane.origin[-1] = depth_m

    # set_trace()
    # embed(colors="neutral")

    return cut


@mlab.show
def plot_topdown(depth_km="surface"):
    """
    Plot a topdown view of the model, allow specification of

    :return:
    """
    try:
        fid = "diff_vp_nz_tall_north_and_vp_cc.vtk"
        assert(os.path.exists(fid))

        # Initiate the engine and reader
        f, engine, vtkfr = startup(fid)

        # Plot the coastline
        if depth_km == "surface":
            z_value = 1000
        else:
            depth_m = -1 * abs(depth_km) * 1E3
            z_value = depth_m
        p3d = coastline("nz_coast_utm60H_43-173_37-179_xyz.npy", z_value)

        # Plot the VTK file
        if depth_km == "surface":
            show_surface(engine, vtkfr, "top-down")
        else:
            cut = cut_plane(engine, vtkfr, depth_m)

        scene = engine.scenes[0]
        scene.scene.z_plus_view()
        cbar = set_colorscale(title="Difference Vs", nb_labels=3)
        # axes = set_axes(xyz=[True, True, False])

        cut.implicit_plane.origin[-1] = depth_m


    except Exception as e:
        print(e)

    mlab.show()

    # if isinstance(depth_km, float):
    #     sug = SliceUnstructuredGrid()
    #     engine.add_filter(sug, module_manager)
    #     sug.implicit_plane.plane.normal = array([0., 0., -1.])
    #     sug.implicit_plane.plane.origin = array([4.02390000e+05,
    #                                              5.59551500e+06,
    #                                              -1 * depth_km * 1E3])


if __name__ == "__main__":
    plot_topdown(depth_km="surface")

