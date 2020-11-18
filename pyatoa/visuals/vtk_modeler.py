"""
Make figures of models from .vtk files using Mayavi.

A class used to plot model visualizations of .vtk files in an automatable
fashion, using the Mayavi library. Standard figure templates include horizontal
and vertical cross sections which show depth slices of the model.
Repeatedly used auxiliary functionality contained in model_tools
"""
import os
import logging
import numpy as np
from numpy import array
from mayavi import mlab
from mayavi.modules.surface import Surface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane
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


class VTKModeler:
    """
    A Class to automatically visualize slices of a model saved in .vtk
    """
    def __init__(self, figsize=(500, 500), zero_origin=False, scale_axes=1,
                 offscreen=False, logging=True, **kwargs):
        """
        :type figsize: tuple of float
        :param figsize: figure size in pixels
        :type scale_axes: float
        :param scale_axes: constant value to multiply against all ranges,
            to e.g. change units from 'm' to 'km'.
        :type zero_origin: bool
        :param zero_origin: sets X and Y origin to 0, Z stays the same
        :type offscreen: bool
        :param offscreen: render mlab offscreen to avoid creating new windows
        :type logging: bool
        :param logging: turn on logging to see what VTKModeler is doing under
            the hood

        Keyword arguments
        ::
            str cmap:
                colorscale to use, matches names from Matplotlib defaults
                to 'RdYlBu' for Red-Yellow-Blue
            bool reverse:
                reverse the colorscale, defaults to True
            list min_max:
                min and max values for color scale as [min, max], defaults to
                None, use the values set by the model
            bool colorbar:
                create and show a colorbar, defaults to True
            str title:
                colorbar title, defaults to None
            int num_clabels:
                number of labels to put on the colorbar, defaults to Mayavi
                default value
            int num_colors:
                number of colors to use for the colorbar, defaults to 20
            int round_to:
                Round the range values of the colorbar to a given base value
                `round_to` to get rid of decimals or weird values.
            str xlabel:
                label for the X-axis, default "E"
            str ylabel:
                label for the Y-axis, default "N"
            str zlabel:
                label for the Z-axis, default "Z"
            list of bool xyz:
                visibility for the x, y, and z axes in a list [x, y, z]
                e.g. [True, True, False] turns off Z axis.
            float font_factor:
                multiply font size by this value for larger or smaller axis
                font sie. Default 1.
            str src_marker:
                marker to be used to show sources. defaults to '2dcircle' for
                depth slice and 'sphere' for cross-section
            str rcv_marker:
                marker to be used to show receivers. defaults to '2ddiamond'
            str rcv_color:
                color to plot receiver markers, defaults 'w'
            str src_color:
                color to plot source markers, defaults 'g'
        """
        # Figure setup
        self.figsize = figsize
        self.scale_axes = scale_axes
        self.zero_origin = zero_origin

        # Empty initiation objects
        self.fig = None
        self.engine = None
        self.vtkfr = None
        self.fid = None
        self.srcs = None
        self.rcvs = None
        self.coast = None
        self.coords = None
        self.ranges = None
        self.axes_ranges = None
        self.kwargs = kwargs

        # Rendering offscreen to avoid popping up windows
        if offscreen:
            mlab.options.offscreen = True

        # Toggle the logger
        if not logging:
            logger.setLevel("INFO")

    def __str__(self):
        """
        Simple string representation
        """
        str_out = ""
        for key, val in self.__dict__.items():
            if not isinstance(val, str) and hasattr(val, "len") and \
                    len(val) > 3:
                val = len(val)
            str_out += f"{key:>12}: {val}\n"
        return str_out

    def load(self, fid=None, src_fid=None, rcv_fid=None, coast_fid=None):
        """
        Load a .VTK file and any auxiliary files like source and receiver
        locations. Uses argument names to determine how to read each file.
        Assumes standard .VTK file structure

        :type fid: str
        :param fid: the file id of the .vtk file to be plotted
        :type src_fid: str
        :param src_fid: (optional) .vtk file id for sources to be plotted
        :type rcv_fid: str
        :param rcv_fid: (optiona) .vtk file id for receivers to be plotted
        :type coast_fid: str
        :param coast_fid: (optional) .npy file id for a coastline to be plotted
        """
        # Load up the main vtk file and read in information from the file
        if fid:
            assert (os.path.exists(fid)), f"File {fid} does not exist"
            self.fid = fid
            self.coords = get_coordinates(self.fid)
            self.ranges = get_ranges(self.coords)
            self.axes_ranges = get_ranges(coords=self.coords,
                                          scale_axes=self.scale_axes,
                                          zero_origin=self.zero_origin)

        # Load in auxiliary data like source and receiver locations
        if src_fid:
            self.srcs = np.loadtxt(src_fid, skiprows=5)
        if rcv_fid:
            self.rcvs = np.loadtxt(rcv_fid, skiprows=5)
        if coast_fid:
            self.coast = np.load(coast_fid)

    def _startup(self):
        """
        Open a VTK file and return the engine that is visualizing it.
        Whenever the figure is shown, startup needs to be called again,
        similar to how showing a matplotlib figure will destroy the instance.
        """
        assert self.fid, "File ID must be specified before generating figure"

        # Instantiate mlab
        self.fig = mlab.figure(size=self.figsize)
        self.engine = mlab.get_engine()
        self.engine.scenes[0].scene.background = colors["w"]
        self.vtkfr = self.engine.open(self.fid)

    def _cut_plane(self, axis, slice_at=None, ratio=None):
        """
        Slice the data at a certain depth, or for a given depth cross section,
        dependent on choice. Format of depth is positive

        :type axis: str
        :param axis: choice of axis to slice along, 'X', 'Y' or 'Z'
        :type slice_at: float
        :param slice_at: depth or distance value in meters to slice at
        :type ratio: float
        :param ratio: ratio of the axis bounds to slice at, from 0 to 1
        """
        assert(slice_at is not None or ratio is not None)
        cut = ScalarCutPlane()
        self.engine.add_filter(cut, self.vtkfr)

        # Set some attributes of the cut plane
        origin = cut.implicit_plane.origin
        if axis == "X":
            cut.implicit_plane.origin = array([slice_at, origin[1], origin[2]])
            cut.implicit_plane.normal = array([1, 0, 0])
        elif axis == "Y":
            cut.implicit_plane.origin = array([origin[0], slice_at, origin[2]])
            cut.implicit_plane.normal = array([0, 1, 0])
        elif axis == "Z":
            cut.implicit_plane.origin = array([origin[0], origin[1], slice_at])
            cut.implicit_plane.normal = array([0, 0, 1])

        cut.implicit_plane.widget.enabled = False

        return cut

    def depth_slice(self, depth_km="surface", save=None, show=True,
                    startup=True, anno_text=None):
        """
        Plot a topdown view of the model, either with a surface projection
        (map view) or at some depth slice determined by `depth_km`

        :type depth_km: float or None
        :param depth_km: depth to show the model at in units of km, by default
            plots a map view of the 'surface'
        :type save: str
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        :type startup: bool
        :param startup: run the _startup() function which closes all instances
            of mlab and creates a new mlab figure and engine. Normally the
            necessary thing to do so defaults to True
        :type anno_text: str
        :param anno_text: text to annotate into the corner of the figure, 
            defaults to annotating the depth value
        """
        src_color = self.kwargs.get("src_color", "g")
        src_marker = self.kwargs.get("src_marker", "2dcircle")
        rcv_color = self.kwargs.get("rcv_color", "w")
        rcv_marker = self.kwargs.get("rcv_marker", "2ddiamond")

        if startup:
            self._startup()

        logger.info(f"Depth Slice (Z) of '{self.fid}' at {depth_km} [km]")
        coastline_z = self.ranges[-1]

        if depth_km != "surface":
            depth_m = -1 * abs(depth_km) * 1E3
            coastline_z = depth_m

        # Put the coastline at some height above the topography, or at depth
        if self.coast is not None:
            coastline(self.coast, coastline_z)

        # Plot the stations and receivers
        if self.rcvs is not None:
            srcrcv(self.rcvs, color=rcv_color, z_value=coastline_z,
                   marker=rcv_marker)
        if self.srcs is not None:
            srcrcv(self.srcs, color=src_color, z_value=coastline_z,
                   marker=src_marker)

        # Plot the model at the given height
        # Axes behave weirdly when cutting planes so turn off Y-axis, not sure
        # why this works... sorry
        if depth_km == "surface":
            # Show a surface projection
            self.engine.add_filter(Surface(), self.vtkfr)
            tag = depth_km
            set_axes(xyz=[True, True, False], ranges=self.axes_ranges,
                     **self.kwargs)
        else:
            # Show a slice (plane) at a given depth
            self._cut_plane(axis="Z", slice_at=depth_m)
            tag = f"depth_{int(abs(depth_km))}km"
            set_axes(xyz=[True, False, True], ranges=self.axes_ranges,
                     **self.kwargs)

        # Set the camera with top down view
        scene = self.engine.scenes[0]
        scene.scene.z_plus_view()

        # Plot extras
        colorscale(orientation="vertical", **self.kwargs)
        if anno_text is None:
            anno_text = tag.replace("_", " ")
        annotate(s=anno_text, c="k", width=0.175)  # was 0.2

        # Finalize
        save_tag = None
        if save:
            mlab.savefig(save.format(tag=tag))
            save_tag = save.format(tag=tag)
        if show:
            mlab.show()
        else:
            mlab.close(self.fig)

        return save_tag

    def cross_section(self, axis, pct=0.5, show_axis=True,
                      show=True, save=None, startup=True):
        """
        Plot a cross-section of the model, with the Y-axis plotting depth

        If axis == X, the slice is normal to the Y-axis
        if axis == Y, the slice is normal to the X-axis

        :type axis: str
        :param axis: axis of axis to slice along, 'X', 'Y'
        :type pct: float
        :param pct: percentage to slice the model at based on the model range.
            defaults to 50%, or the middle of the model
        :type show_axis: bool
        :param show_axis: show the axis of the plot
        :type save: str
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        :type startup: bool
        :param startup: run the _startup() function which closes all instances
            of mlab and creates a new mlab figure and engine. Normally the
            necessary thing to do so defaults to True
        """
        src_color = self.kwargs.get("src_color", "g")
        src_marker = self.kwargs.get("src_marker", "sphere")
        rcv_color = self.kwargs.get("rcv_color", "w")
        rcv_marker = self.kwargs.get("rcv_marker", "2ddiamond")

        if startup:
            self._startup()

        axis = axis.upper()

        def slice_range(ranges_):
            """
            Convenience function to calculate where to slice by checking the
            min and max of the range and multiplying by a percentage

            :type ranges_: list
            :param ranges_: [xmin, xmax, ymin, ymax, zmin, zmax]
            """
            if axis == "X":
                return pct * (ranges_[1] - ranges_[0]) + ranges_[0]
            elif axis == "Y":
                return pct * (ranges_[3] - ranges_[2]) + ranges_[2]

        # Determine what part of the axis we are slicing
        slice_val = slice_range(self.ranges)
        tag = slice_range(self.axes_ranges)

        logger.info(
                f"Cross section [{axis}] of '{self.fid}'  at {slice_val*1E-3}")
        if show_axis:
            set_axes(xyz=[True, True, False], ranges=self.axes_ranges,
                     **self.kwargs)

        self._cut_plane(axis=axis, slice_at=slice_val)

        # Differentiate which way were slicing to put the sources and receivers
        slice_x, slice_y = None, None
        if axis == "X":
            slice_x = slice_val
        elif axis == "Y":
            slice_y = slice_val

        # Plot the sources and receivers along a given slice
        if self.rcvs is not None:
            srcrcv(self.rcvs, color=rcv_color, x_value=slice_x, y_value=slice_y,
                   marker=rcv_marker)
        if self.srcs is not None:
            # Sources default to spheres because I didn't want to figure out how
            # to rotate 2d glyphs
            srcrcv(self.srcs, color=src_color, x_value=slice_x, y_value=slice_y,
                   marker=src_marker)

        # Set the camera with a side on view. No preset so set to the preset and
        # rotate to get to the proper frame of reference
        scene = self.engine.scenes[0]
        if axis == "X":
            scene.scene.x_plus_view()
        elif axis == "Y":
            scene.scene.y_plus_view()
            scene.scene.camera.roll(90)

        # Plot extras
        colorscale(orientation="horizontal", **self.kwargs)

        # Annotation location changes based on which axis you slice
        anno_tag = f"{axis.lower()}={pct*1E2}%={tag:.2f}km"
        if axis == "X":
            annotate(x=.2, y=.3, s=anno_tag, c="k", width=0.15)
        elif axis == "Y":
            annotate(x=.27, y=.3, s=anno_tag, c="k", width=0.15)

        # Finalize
        save_tag = None
        if save:
            save_tag = save.format(tag=f"{tag:.0f}km")
            mlab.savefig(save_tag)
        if show:
            mlab.show()
        else:
            mlab.close(self.fig)

        return save_tag


def get_coordinates(fid):
    """
    I couldn't, for the life of me, figure out how to get the mesh dimensions
    out of the mayavi objects, so instead just get it from using the VTK lib

    :type fid: str
    :param fid: file id of the .vtk file
    :rtype: numpy array
    :return: Nx3 numpy array with columns corresponding to x, y, z
    """
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(fid)
    reader.Update()
    coords = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())

    return coords


def get_ranges(coords, scale_axes=False, zero_origin=False):
    """
    Get a list of ranges to use for setting Axes elements by generating an empty
    axis object. A bit hacky but it works within the context of Mayavi.

    Option to scale the axes to the units of range by some constant value.
    Option to set the origin to zero.

    :type coords: np.array
    :param coords: Nx4 numpy array with columns representing x, y, z, v, the
        output of get_coordinates()
    :type scale_axes: float
    :param scale_axes: constant value to multiply against all ranges, to e.g.
        change units from 'm' to 'km'
    :type zero_origin: bool
    :param zero_origin: sets X and Y origin to 0, Z stays the same
    :rtype: list of floats
    :return: [xmin, xmax, ymin, ymax, zmin, zmax]
    """
    xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
    ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
    zmin, zmax = coords[:, 2].min(), coords[:, 2].max()

    # Set the ranges
    if zero_origin:
        ranges = [0, xmax-xmin, 0, ymax-ymin, zmin, zmax]
    else:
        ranges = [xmin, xmax, ymin, ymax, zmin, zmax]

    # Convert the range
    if scale_axes:
        ranges = [_ * scale_axes for _ in ranges]

    return ranges


def colorscale(orientation, **kwargs):
    """
    Utility function to set the colorscale for the plot and also create a
    colorbar

    :type orientation: str
    :param orientation: 'vertical' or 'horizontal'

    Keyword arguments
    ::
        str cmap: 
            colorscale to use, matches names from Matplotlib
            defaults to 'RdYlBu' for Red-Yellow-Blue
        bool reverse: 
            reverse the colorscale, defaults to True
        list min_max: 
            min and max values for color scale as [min, max],
            defaults to None, use the values set by the model
        bool colorbar: 
            create and show a colorbar, defaults to True
        str cbar_title: 
            colorbar title, defaults to None
        int num_clabels: 
            number of labels to put on the colorbar, defaults to Mayavi default
        int round_to: 
            Round the range values of the colorbar to a given
            base value `round_to` to get rid of decimals or weird values.
    """
    cmap = kwargs.get("cmap", "RdYlBu")
    reverse = kwargs.get("reverse", True)
    default_range = kwargs.get("default_range", False)
    min_max = kwargs.get("min_max", None)
    colorbar = kwargs.get("colorbar", True)
    cbar_title = kwargs.get("cbar_title", None)
    num_labels = kwargs.get("num_clabels", None)
    num_colors = kwargs.get("num_colors", 20)
    round_to = kwargs.get("round_to", 0)

    # Set the colorscale
    cbar = mlab.colorbar(title=cbar_title, orientation=orientation,
                         label_fmt="%-#.2f", nb_labels=num_labels)
    cbar.lut_mode = cmap
    logger.debug(f"Creating colorbar with colormap: '{cmap}'")
    if reverse:
        cbar.reverse_lut = reverse
        logger.debug(f"Reversing colorscale")
    cbar.number_of_colors = num_colors

    # Out of bounds colors !!! This doesnt work, porque no?
    # cbar.lut.use_below_range_color = 1
    # cbar.lut.above_range_color = array([-1., -1., -1., -1.,])
    # cbar.lut.use_above_range_color = 1
    # cbar.lut.below_range_color = array([1., 1., 1., 1.,])

    # Bound the colorscale
    logger.debug(f"Colorbar bounds currently set to: {cbar.data_range}")
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
            maxval = max(abs(cbar.data_range))
            minval = min(cbar.data_range)
            if minval > 0:
                min_max = [0, maxval]
            elif maxval < 0:
                min_max = [minval, 0]
            else:
                val = max([abs(_) for _ in [minval, maxval]])
                min_max = [-val, val]
        cbar.data_range = array(min_max)

    # Use 'E' notation for small values like gradients
    if max(cbar.data_range) < .01:
        cbar.scalar_bar.label_format = "%-#.2E"

    logger.debug(f"Colorbar bounds have been set to: {cbar.data_range}")

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
            cbar.scalar_bar_representation.position = array([0.2, .04])
            cbar.scalar_bar_representation.position2 = array([0.6, 0.12])

        # Vertical scalebar sits to the right of the figure, for top down
        elif orientation == "vertical":
            cbar.scalar_bar.title_ratio = 0.36  # smaller title size
            cbar.scalar_bar.bar_ratio = 0.36  # thickness of colorbar
            cbar.scalar_bar_representation.position = array([0.8, 0.1])
            cbar.scalar_bar_representation.position2 = array([0.125, .8])

    return cbar


def set_axes(xlabel="E", ylabel="N", zlabel="Z", ranges=None,
             xyz=[True, True, True], **kwargs):
    """
    Utility function to create an axis object that wraps around data.
    Changing the ranges on the axes only changes labels, doesn't affect e.g.
    plotting other points onto the axis

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
                     ranges=ranges
                     )

    # Edit full properties
    axes.axes.label_format = '%-#.2f'
    axes.axes.number_of_labels = 6
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


def coastline(coords, z_value=1000, color="k"):
    """
    Plot coastline on top of plot. Coastline should be an npy file that is an
    N x 3 array with the columns representing x, y, z

    :type coords: np.array
    :param coords: Nx3 array with columns relating to x, y z
    :type z_value: float
    :param z_value: height of the coastline in meters, negative down
    :type color: str
    :param color: color of the coastline
    """
    x = coords[:, 0]
    y = coords[:, 1]
    z = np.ones(len(x)) * z_value

    p3d = mlab.points3d(x, y, z, color=colors[color], mode="2dcircle",
                        scale_factor=800.)

    return p3d


def srcrcv(coords, x_value=None, y_value=None, z_value=None, color="w",
           marker="2ddiamond"):
    """
    Take a receivers VTK file, outputted by Pyatoa, and plot it ontop of
    the current projection. Allows for condensing all stations into a single
    plane for e.g. depth slice plots

    .. note::
        Available markers: 2darrow, 2dcircle, 2dcross, 2ddash, 2ddiamond,
        2dhooked_arrow, 2dsquarem, 2dthick_arrow, 2dthick_cross, 2dtriangle,
        2dvertex, arrow, axes, cone, cube, cylinder, point, sphere

    :type coords: np.array
    :param coords: Nx3 array with columns relating to x, y z
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
