"""
Static helper functions to create model visualization objects in Mayavi pipeline
"""
import logging
import numpy as np
from numpy import array
from mayavi import mlab
from pyatoa.utils.calculate import myround


# Set the logger as a globally accessible variable
logger = logging.getLogger("vtk_plotter")
logger.setLevel("debug".upper())  # Default level
logger.propagate = 0  # Prevent propagating to higher loggers
ch = logging.StreamHandler()  # Console log handler
FORMAT = "%(message)s"
formatter = logging.Formatter(FORMAT)  # Set format of logging messages
ch.setFormatter(formatter)
logger.addHandler(ch)


# Mayavi takes colors in RGB color space with values from 0 to 1,
# Lookup table for common colors based on Python naming conventions
colors = {"k": (0., 0., 0.), "w": (1., 1., 1.), "r": (1., 0., 0.),
          "b": (0., 0., 1.), "g": (0., 1., 0.), "y": (1., 1., 0.),
          "o": (1., .5, 0.), "c": (0., 1., 1.),
          "gray": (.5, .5, .5)
          }


def get_coordinates(fid):
    """
    I couldn't, for the life of me, figure out how to get the mesh dimensions
    out of the mayavi objects, so instead just get it from using the VTK lib

    :type fid: str
    :param fid: file id of the .vtk file
    :rtype: numpy array
    :return: Nx4 numpy array with columns corresponding to x, y, z, v
    """
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(fid)
    reader.Update()
    coords = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
    values = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())

    arr_out = np.c_[coords, values]

    return arr_out


def colorscale(orientation, **kwargs):
    """
    Utiltiy function to set the colorscale for the plot and also create a
    colorbar

    :type orientation: str
    :param orientation: 'vertical' or 'horizontal'

    Keyword arguments:
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
    cmap = kwargs.get("cmap", "RdYlBu")
    reverse = kwargs.get("reverse", True)
    default_range = kwargs.get("default_range", False)
    min_max = kwargs.get("min_max", None)
    colorbar = kwargs.get("colorbar", True)
    title = kwargs.get("title", None)
    num_labels = kwargs.get("num_clabels", None)
    num_colors = kwargs.get("num_colors", 20)
    round_to = kwargs.get("round_to", 0)

    # Set the colorscale
    cbar = mlab.colorbar(title=title, orientation=orientation,
                         label_fmt="%-#.2f", nb_labels=num_labels)
    cbar.lut_mode = cmap
    logger.debug(f"Creating colorbar with colormap {cmap}")
    cbar.reverse_lut = reverse
    if reverse:
        logger.debug(f"Reversing colorscale")
    cbar.number_of_colors = num_colors

    # Bound the colorscale
    logger.info(f"Data bounds are set to {cbar.data_range}")
    if default_range:
        # Round the default min and max bounds of the dataset
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
            # If no min_max values, set absolute max as the +/- bounds
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

        # Horizontal scalebar sits at the bottom, for depth slice plots
        if orientation == "horizontal":
            # Create some uniform look for colorbar
            cbar.scalar_bar.text_position = 'precede_scalar_bar'  # below
            cbar.scalar_bar.title_ratio = 0.5  # smaller title size
            cbar.scalar_bar.bar_ratio = 0.4  # thickness of colorbar

            # Makes the colorbar interactive and changes its size
            cbar.scalar_bar_representation.proportional_resize = True
            cbar.scalar_bar_representation.position = array([0.12, .03])
            cbar.scalar_bar_representation.position2 = array([0.75, 0.15])

        # Vertical scalebar sits to the right of the figure, for top down
        elif orientation == "vertical":
            cbar.scalar_bar.title_ratio = 0.36  # smaller title size
            cbar.scalar_bar.bar_ratio = 0.325  # thickness of colorbar
            cbar.scalar_bar_representation.position = array([0.8, 0.1])
            cbar.scalar_bar_representation.position2 = array([0.125, .8])

    return cbar


def set_axes(xlabel="E (m)", ylabel="N (m)", zlabel="Z (m)",
             xyz=[True, True, True], **kwargs):
    """
    Utility function to create an axis object that wraps around data

    :type xlabel: str
    :param xlabel: label for the X-axis
    :type ylabel: str
    :param ylabel: label for the Y-axis
    :type zlabel: str
    :param zlabel: label for the Z-axis
    :type xyz: list of bool
    :param xyz: visibility for the x, y, and z axes in a list
    """
    # Font size as a factor
    font_factor = kwargs.get("font_factor", 1.)

    # Remove labels if axis is turned off
    if not xyz[0]:
        xlabel = ""
    if not xyz[1]:
        ylabel = ""
    if not xyz[2]:
        zlabel = ""

    axes = mlab.axes(color=(0., 0., 0.,), line_width=3., xlabel=xlabel,
                     ylabel=ylabel, zlabel=zlabel, x_axis_visibility=xyz[0],
                     y_axis_visibility=xyz[1], z_axis_visibility=xyz[2],
                     )

    # Edit full properties
    axes.axes.label_format = '%-#.2f'
    axes.axes.number_of_labels = 3
    axes.axes.font_factor = font_factor

    # Edit label attributes
    axes.label_text_property.color = (0.0, 0.0, 0.0)  # black font
    axes.label_text_property.bold = False
    axes.label_text_property.italic = False
    axes.label_text_property.orientation = 0.0

    # Edit title attributes
    axes.title_text_property.color = (0.0, 0.0, 0.0)  # black font
    axes.title_text_property.italic = False
    axes.title_text_property.bold = False
    axes.title_text_property.orientation = 0.0

    return axes


def coastline(coast_fid, z_value=1000, color="k"):
    """
    Plot coastline on top of plot. Coastline should be an npy file that is an
    N x 3 array with the columns representing x, y, z

    :type coast_fid: str
    :param coast_fid: fid for npy file
    :type z_value: float
    :param z_value: height of the coastline in meters, negative down
    :type color: str
    :param color: color of the coastline
    """
    coast = np.load(coast_fid)
    x = coast[:, 0]
    y = coast[:, 1]
    z = np.ones(len(x)) * z_value

    p3d = mlab.points3d(x, y, z, color=colors[color], mode="2dcircle",
                        scale_factor=800.)

    return p3d


def srcrcv(fid, x_value=None, y_value=None, z_value=None, color="w",
           marker="2ddiamond"):
    """
    Take a receivers VTK file, outputted by Pyatoa, and plot it ontop of
    the current projection. Allows for condensing all stations into a single
    plane for e.g. depth slice plots

    Available markers:
     ‘2darrow’ or ‘2dcircle’ or ‘2dcross’ or ‘2ddash’ or ‘2ddiamond’ or
     ‘2dhooked_arrow’ or ‘2dsquare’ or ‘2dthick_arrow’ or ‘2dthick_cross’ or
     ‘2dtriangle’ or ‘2dvertex’ or ‘arrow’ or ‘axes’ or ‘cone’ or ‘cube’ or
     ‘cylinder’ or ‘point’ or ‘sphere’.

    :type fid: str
    :param fid: file id of the source or receiver .vtk files
    :type x_value: float
    :param x_value: if plotting on a plane, collapses axis to a single value
    :type y_value: float
    :param y_value: if plotting on a plane, collapses axis to a single value
    :type z_value: float
    :param z_value: if plotting on a plane, collapses axis to a single value
    :type color: str
    :param color: color of the points
    :type marker: str
    :param marker: marker to use for the points
    """
    coords = np.loadtxt(fid, skiprows=5)
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]

    # If preset values given, overwrite coordinates
    if x_value:
        x = x_value * np.ones(len(x))
    if y_value:
        y = y_value * np.ones(len(y))
    if z_value:
        z = z_value * np.ones(len(z))

    # Plot the stations as points
    p3d = mlab.points3d(x, y, z, color=colors[color], mode=marker,
                        scale_factor=12500., line_width=3.,)

    return p3d


def annotate(s, x=0.225, y=0.825, c="k", width=0.3):
    """
    Annotate text onto the axis. Font size doesnt work, width to control size

    :type s: str
    :param s: string to annotate
    :type x: float
    :param x: part of the X-axis to plot from [0, 1]
    :type y: float
    :param y: part of the Y-axis to plot from [0, 1]
    :type c: str
    :param c: color of the text
    :type width: float
    :param width: size of the font
    """
    return mlab.text(x, y, s, width=width, color=colors[c])

