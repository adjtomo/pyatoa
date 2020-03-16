"""
Tools to use the Python package mayavi to interact with .vtk files.
These were created because I was tired of spending all my time fine tuning
plots in Paraview, only to remake them all over again. Hopefull this set of
tools allows for quick, standard looking plots of vtk files for easy
visualization of kernels, models, etc.
"""
import os
import logging
import numpy as np
from numpy import array
from mayavi import mlab
from mayavi.modules.surface import Surface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane
from mayavi.modules.slice_unstructured_grid import SliceUnstructuredGrid

from pyatoa.utils.tools.calculate import myround

from ipdb import set_trace
from IPython import embed

# Set the logger as a globally accessible variable
logger = logging.getLogger("vtk_plotter")
logger.setLevel("info".upper())  # Default level
logger.propagate = 0  # Prevent propagating to higher loggers
ch = logging.StreamHandler()  # Console log handler
FORMAT = "%(message)s"
formatter = logging.Formatter(FORMAT)  # Set format of logging messages
ch.setFormatter(formatter)
logger.addHandler(ch)


def startup(fid, figsize=None):
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
    logger.debug(f"Reading {fid} and creating figure size {figsize}")
    # Instantiate mlab
    fig = mlab.figure(size=figsize)
    engine = mlab.get_engine()
    engine.scenes[0].scene.background = (1.0, 1.0, 1.0)  # background to white
    vtk_file_reader = engine.open(fid)

    return fig, engine, vtk_file_reader


def set_colorscale(cmap='RdYlBu', reverse=False, default_range=False,
                   min_max=None,  colorbar=True, title=None, nb_labels=None,
                   orientation="horizontal", num_col=20, round_to=0,
                   **kwargs):
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
    logger.debug(f"Creating colorbar with colormap {cmap}")
    cbar.reverse_lut = reverse
    if reverse:
        logger.debug(f"Reversing colorscale")
    cbar.number_of_colors = num_col
    cbar.number_of_labels = num_col // 2 + 1

    # Bound the colorscale
    logger.info(f"Data bounds are set to {cbar.data_range}")
    if default_range:
        # Round the default min and max bounds of the dataset for cleaner values
        if round_to:
            min_max = cbar.data_range
            cbar.data_range = array(
                [myround(min_max[0], base=round_to, choice="down"),
                 myround(min_max[1], base=round_to, choice="up")]
            )
        # Use the straight up min and max bounds of the dataset
        else:
            cbar.use_default_range = default_range
    # If not using the default range
    else:
        if not min_max:
            # If no min_max values, set the absolute max value as the +/- bounds
            val = np.ceil(max(abs(cbar.data_range)))
            min_max = [-val, val]
        cbar.data_range = array(min_max)
    logger.info(f"Data bounds have been set to to {cbar.data_range}")

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

        cbar.scalar_bar.orientation = orientation

        # Horizontal scalebar sits at the bottom
        if orientation == "horizontal":
            # Create some uniform look for colorbar
            cbar.scalar_bar.text_position = 'precede_scalar_bar'  # title below
            cbar.scalar_bar.title_ratio = 0.36  # smaller title size
            cbar.scalar_bar.bar_ratio = 0.325  # thickness of colorbar

            # Makes the colorbar interactive and changes its size
            cbar.scalar_bar_representation.proportional_resize = True
            cbar.scalar_bar_representation.position = array([0.33, .0005])  # location
            cbar.scalar_bar_representation.position2 = array([0.33, 0.062])  # size

        # Vertical scalebar sits to the right of the figure
        elif orientation == "vertical":
            cbar.scalar_bar.title_ratio = 0.36  # smaller title size
            cbar.scalar_bar.bar_ratio = 0.325  # thickness of colorbar
            cbar.scalar_bar_representation.position = array([0.8, 0.1])
            cbar.scalar_bar_representation.position2 = array([0.125, .8])

    return cbar


def set_axes(xlabel="E (m)", ylabel="N (m)", zlabel="Z (m)",
             xyz=[True, True, True], **kwargs):
    """
    Create an axis object to wrap around data
    :return:
    """
    # Font size as a factor
    font_factor = kwargs.get("font_factor", 1.)

    axes = mlab.axes(color=(0., 0., 0.,), line_width=3., xlabel=xlabel,
                     ylabel=ylabel, zlabel=zlabel, x_axis_visibility=xyz[0],
                     y_axis_visibility=xyz[1], z_axis_visibility=xyz[2])

    # Edit full properties
    axes.axes.label_format = '%-#3.2E'
    axes.axes.number_of_labels = 5
    axes.axes.font_factor = font_factor

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


def show_surface(engine, vtk_file_reader):
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
def plot_model_surface(vtkfr, engine, fid, coastline_fid, save=False, show=True,
                       **kwargs):
    """
    Plot a topdown view of the model, with the projection at the surface
    """
    # Put the coastline at some height above the topography
    if coastline_fid:
        p3d = coastline(coastline_fid, 1000)

    # Plot the VTK file
    show_surface(engine, vtkfr)

    # Set the camera with top down view
    scene = engine.scenes[0]
    scene.scene.z_plus_view()

    # Colorbar and axes
    set_colorscale(**kwargs)
    set_axes(xyz=[True, True, False], **kwargs)

    if save:
        mlab.savefig(save.format(tag="surface"))

    if show:
        mlab.show()


# def plot_model_cuts(fid, depth_list=[2, 5, 10, 15, 20, 25, 30, 50])
#     # Initiate the engine and reader
#     f, engine, vtkfr = startup(fid)
#     for depth in depth_list:
#         depth_m = -1 * abs(depth_km) * 1E3
#         z_value = depth_m
#     p3d = coastline("nz_coast_utm60H_43-173_37-179_xyz.npy", z_value)
#
#     cut = cut_plane(engine, vtkfr, depth_m)


if __name__ == "__main__":
    """                         Set your parameters below                    """
    # ID's for files to plot
    fid = "vs_nz_tall_north.vtk"
    coast = "nz_coast_utm60H_43-173_37-179_xyz.npy"

    # Figure parameters
    show = True
    save_fid = "vs_nz_tall_north_{tag}.png"
    figsize = (1000, 1000)

    # Axes parameters
    font_factor = 1.1

    # Colormap parameters
    colormap = "jet"
    reverse = True
    cbar_title = "Vs (m/s)"
    cbar_orientation = "vertical"
    number_of_colors = 35
    default_range = False
    round_to = 50  # only if default_range == True
    min_max = [1200, 3500]  # only if default_range == False, []

    """                         Set your parameters above                    """

    # Initiate the engine and reader
    fig, engine, vtkfr = startup(fid, figsize)

    plot_model_surface(vtkfr, engine, fid=fid, coastline_fid=coast,
                       cmap=colormap, reverse=reverse, title=cbar_title,
                       num_col=number_of_colors, min_max=min_max,
                       default_range=default_range,
                       font_factor=font_factor, orientation=cbar_orientation,
                       round_to=round_to, save=save_fid, show=show)

