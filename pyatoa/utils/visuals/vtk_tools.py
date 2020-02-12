"""
Tools to use the Python package mayavi to interact with .vtk files.
These were created because I was tired of spending all my time fine tuning
plots in Paraview, only to remake them all over again. Hopefull this set of
tools allows for quick, standard looking plots of vtk files for easy
visualization of kernels, models, etc.
"""
import numpy as np
from numpy import array
from mayavi import mlab
from mayavi.modules.surface import Surface
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
    :return engine: engine that is visualizing the vtk file
    """
    # Instantiate mlab
    fig = mlab.figure()
    engine = mlab.get_engine()
    engine.scenes[0].scene.background = (1.0, 1.0, 1.0)  # background to white
    vtk_file_reader = engine.open(fid)

    # Add data and standard objects to engine
    surface = Surface()
    engine.add_filter(surface, vtk_file_reader)

    return fig, engine


def set_colorscale(lut_mode='RdYlBu', reverse=False, min_max=None,
                   colorbar=True, title=None, nb_labels=None):
    """
    Set the colorscale for the plot and also create a colorbar

    :type lut_mode: str
    :param lut_mode: colorscale to use, matches names from ParaView
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
    slm = mlab.colorbar(title=title, orientation="horizontal",
                        label_fmt="%-#3.2e", nb_labels=nb_labels)
    slm.lut_mode = lut_mode
    slm.reverse_lut = reverse

    # Bound the colorscale
    slm.use_default_range = False
    if not min_max:
        val = max(abs(slm.data_range))
        min_max = [-val, val]
    slm.data_range = array(min_max)

    # Create colorbar
    slm.show_scalar_bar = colorbar
    slm.show_legend = colorbar

    if colorbar:
        # These are for the values
        slm.label_text_property.bold = False
        slm.label_text_property.italic = False
        slm.label_text_property.orientation = 25.  # rotate labels slightly
        slm.label_text_property.color = (0.0, 0.0, 0.0)  # black font

        # These are for the title
        slm.title_text_property.italic = False
        slm.title_text_property.bold = False
        slm.title_text_property.color = (0.0, 0.0, 0.0)  # black font

        # Create some uniform look for colorbar
        slm.scalar_bar.orientation = 'horizontal'  # horizontal scalebar
        slm.scalar_bar.text_position = 'precede_scalar_bar'  # title below
        slm.scalar_bar.title_ratio = 0.2  # smaller title size
        slm.scalar_bar.bar_ratio = 0.25  # size of colorbar

    return slm


def set_axes(xlabel="Easting (m)", ylabel="Northing (m)", zlabel="Depth (m)",
             xyz=[True, True, True]):
    """
    Create an axis object to wrap around data
    :return:
    """
    axes = mlab.axes(color=(0., 0., 0.,), line_width=3., xlabel=xlabel,
                     ylabel=ylabel, zlabel=zlabel, x_axis_visibility=xyz[0],
                     y_axis_visibility=xyz[1], z_axis_visibility=xyz[2])

    # Edit label attributes
    axes.axes.label_format = '%-#3.2e'
    axes.axes.number_of_labels = 5
    axes.label_text_property.bold = False
    axes.label_text_property.italic = False
    axes.label_text_property.orientation = 45.0

    # Edit title attributes
    axes.title_text_property.color = (0.0, 0.0, 0.0)  # black font
    axes.title_text_property.italic = False
    axes.title_text_property.bold = False
    axes.title_text_property.orientation = 45.0
    axes.title_text_property.font_size = 45


def coastline(coast_fid):
    """
    Plot coastline on top of plot

    :param coast_fid:
    :return:
    """
    coast = np.load(coast_fid)

@mlab.show
def plot_topdown(depth_km="surface"):
    """
    Just plot the surface of the

    :return:
    """
    f, engine = startup("diff_vs_nz_tall_north_and_vs_checker.vtk")
    set_colorscale(title="Difference Vs")
    set_axes(xyz=[True, True, False])
    mlab.show()

    set_trace()
    embed(colors="neutral")

    # if isinstance(depth_km, float):
    #     sug = SliceUnstructuredGrid()
    #     engine.add_filter(sug, module_manager)
    #     sug.implicit_plane.plane.normal = array([0., 0., -1.])
    #     sug.implicit_plane.plane.origin = array([4.02390000e+05,
    #                                              5.59551500e+06,
    #                                              -1 * depth_km * 1E3])


if __name__ == "__main__":
    plot_topdown()

