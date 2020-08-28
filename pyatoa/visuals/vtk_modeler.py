"""
Make figures of models from .vtk files using Mayavi.

A class used to plot model visualizations of .vtk files in an automatable
fashion, using the Mayavi library. Standard figure templates include horizontal
and vertical cross sections which show depth slices of the model.
Repeatedly used auxiliary functionality contained in model_tools
"""
import os
import time
import numpy as np
from numpy import array
from mayavi import mlab
from mayavi.modules.surface import Surface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane
from pyatoa.visuals.model_tools import (logger, get_coordinates, set_axes,
                                        annotate, coastline, colorscale,
                                        srcrcv, colors, get_ranges)


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

        COLORSCALE keyword arguments:
            :type cmap: str
            :param cmap: colorscale to use, matches names from Matplotlib
                defaults to 'RdYlBu' for Red-Yellow-Blue
            :type reverse: bool
            :param reverse: reverse the colorscale, defaults to True
            :type min_max: list
            :param min_max: min and max values for color scale as [min, max],
                defaults to None, use the values set by the model
            :type colorbar: bool
            :param colorbar: create and show a colorbar, defaults to True
            :type title: str
            :param title: colorbar title, defaults to None
            :type num_clabels: int
            :param num_clabels: number of labels to put on the colorbar,
                defaults to Mayavi default
            :type num_colors: int
            :param num_colors: number of colors to use for the colorbar,
                defaults to 20
            :type round_to: int
            :param round_to: Round the range values of the colorbar to a given
                base value `round_to` to get rid of decimals or weird values.

        AXIS keyword arguments:
            :type xlabel: str
            :param xlabel: label for the X-axis, default "E"
            :type ylabel: str
            :param ylabel: label for the Y-axis, default "N"
            :type zlabel: str
            :param zlabel: label for the Z-axis, default "Z"
            :type xyz: list of bool
            :param xyz: visibility for the x, y, and z axes in a list [x, y, z]
                e.g. [True, True, False] turns off Z axis.
            :type font_factor: float
            :param font_factor: multiply font size by this value for larger or
                smaller axis font sie. Default 1.

        PLOTTING keyword arguments:
            :type src_marker: str
            :param src_marker: marker to be used to show sources. defaults to
                '2dcircle' for depth slice and 'sphere' for cross-section
            :type rcv_marker: str
            :param rcv_marker: marker to be used to show receivers.
                defaults to '2ddiamond'
            :type rcv_color: str
            :param rcv_color: color to plot receiver markers, defaults 'w'
            :type src_color: str
            :param src_color: color to plot source markers, defaults 'g'
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

        if self.fig is not None:
            time.sleep(1)
            mlab.close(self.fig)

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
        assert(slice_at or ratio)
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

    def depth_slice(self, depth_km="surface", save=False, show=True,
                    startup=True):
        """
        Plot a topdown view of the model, either with a surface projection
        (map view) or at some depth slice determined by `depth_km`

        Y|
         |
         |_______
                X

        :type depth_km: float or None
        :param depth_km: depth to show the model at in units of km, by default
            plots a map view of the 'surface'
        :type save: bool
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        :type startup: bool
        :param startup: run the _startup() function which closes all instances
            of mlab and creates a new mlab figure and engine. Normally the
            necessary thing to do so defaults to True
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
        annotate(s=f"{tag.replace('_', ' ')}", c="k", width=0.2)

        # Finalize
        save_tag = None
        if save:
            mlab.savefig(save.format(tag=tag))
            save_tag = save.format(tag=tag)
        if show:
            mlab.show()

        return save_tag

    def cross_section(self, axis, pct=0.5, show_axis=True,
                      show=True, save=False, startup=True):
        """
        Plot a cross-section of the model, with the Y-axis plotting depth

        Z|
         |
         |________
                 X or Y

        If axis == X, the slice is normal to the Y-axis
        if axis == Y, the slice is normal to the X-axis

        :type axis: str
        :param axis: axis of axis to slice along, 'X', 'Y'
        :type pct: float
        :param pct: percentage to slice the model at based on the model range.
            defaults to 50%, or the middle of the model
        :type show_axis: bool
        :param show_axis: show the axis of the plot
        :type save: bool
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

        return save_tag

