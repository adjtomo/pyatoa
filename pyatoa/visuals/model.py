"""
A class used to plot model visualizations of .vtk files in an automatable
fashion, using the Mayavi library. Standard figure templates include horizontal
and vertical cross sections which show depth slices of the model.
Repeatedly used auxiliary functionality contained in model_tools
"""
import os
import numpy as np
from numpy import array
from mayavi import mlab
from mayavi.modules.surface import Surface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane
from pyatoa.visuals.model_tools import (logger, get_coordinates, set_axes,
                                        annotate, coastline, colorscale,
                                        srcrcv, colors, get_ranges)


class Model:
    """
    A Class to automatically visualize slices of a model saved in .vtk
    """
    def __init__(self, fid, srcs=None, rcvs=None, coast=None,
                 figsize=(1000, 1000), zero_origin=False, convert=1,
                 offscreen=False, logging=True):
        """
        :type fid: str
        :param fid: the file id of the .vtk file to be plotted
        :type srcs: str
        :param srcs: (optional) .vtk file id for sources to be plotted
        :type rcvs: str
        :param rcvs: (optiona) .vtk file id for receivers to be plotted
        :type coast: str
        :param coast: (optional) .npy file id for a coastline to be plotted
        :type figsize: tuple of float
        :param figsize: figure size
        """
        self.coords = get_coordinates(fid)
        self.ranges = get_ranges(self.coords)
        self.convert = convert
        self.axes_ranges = get_ranges(self.coords, convert=self.convert,
                                      zero_origin=zero_origin)

        # Distribute pathnames
        self.fid = fid
        assert(os.path.exists(self.fid)), f"File {fid} does not exist"
        self.srcs, self.rcvs, self.coast = None, None, None
        if srcs:
            self.srcs = np.loadtxt(srcs, skiprows=5)
        if rcvs:
            self.rcvs = np.loadtxt(rcvs, skiprows=5)
        if coast:
            self.coast = np.load(coast)

        # Initiate the plotting engine
        self.figsize = figsize
        self.fig, self.engine, self.vtkfr = None, None, None

        # Rendering offscreen to avoid popping up windows
        if offscreen:
            mlab.options.offscreen = True

        # Toggle the logger
        if not logging:
            logger.setLevel("INFO")

    def startup(self):
        """
        Open a VTK file and return the engine that is visualizing it.
        Whenever the figure is shown, startup needs to be called again,
        similar to how showing a matplotlib figure will destroy the instance.
        """
        mlab.close(all=True)
        # Instantiate mlab
        self.fig = mlab.figure(size=self.figsize)
        self.engine = mlab.get_engine()
        self.engine.scenes[0].scene.background = colors["w"]
        self.vtkfr = self.engine.open(self.fid)

    def cut_plane(self, choice, slice_at=None, ratio=None):
        """
        Slice the data at a certain depth, or for a given depth cross section,
        dependent on choice. Format of depth is positive

        :type choice: str
        :param choice: choice of axis to slice along, 'X', 'Y' or 'Z'
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
        if choice == "X":
            cut.implicit_plane.origin = array([slice_at, origin[1], origin[2]])
            cut.implicit_plane.normal = array([1, 0, 0])
        elif choice == "Y":
            cut.implicit_plane.origin = array([origin[0], slice_at, origin[2]])
            cut.implicit_plane.normal = array([0, 1, 0])
        elif choice == "Z":
            cut.implicit_plane.origin = array([origin[0], origin[1], slice_at])
            cut.implicit_plane.normal = array([0, 0, 1])

        cut.implicit_plane.widget.enabled = False

        return cut

    def plot_model_topdown(self, depth_km=None, save=False, show=True,
                           **kwargs):
        """
        Plot a topdown view of the model, with the projection at the surface

        kwargs are passed to colorscale()

        :type depth_km: float or None
        :param depth_km: depth to show the model at in units of km, if None is
            given, will plot a surface projection
        :type save: bool
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        """
        logger.info(f"slicing '{self.fid}' across Z at {depth_km}km")
        coastline_z = self.ranges[-1]
        if depth_km != "surface":
            depth_m = -1 * abs(depth_km) * 1E3
            coastline_z = depth_m

        # Put the coastline at some height above the topography, or at depth
        if self.coast is not None:
            coastline(self.coast, coastline_z)

        # Plot the stations and receivers
        if self.rcvs is not None:
            srcrcv(self.rcvs, color="w", z_value=coastline_z,
                   marker="2ddiamond")
        if self.srcs is not None:
            srcrcv(self.srcs, color="g", z_value=coastline_z,
                   marker="2dcircle")

        # Plot the model at the given height
        # Axes behave weirdly when cutting planes so turn off Y-axis, not sure
        # why this works... sorry
        if depth_km == "surface":
            # Show a surface projection
            self.engine.add_filter(Surface(), self.vtkfr)
            tag = depth_km
            set_axes(xyz=[True, True, False], ranges=self.axes_ranges, **kwargs)
        else:
            # Show a slice (plane) at a given depth
            self.cut_plane(choice="Z", slice_at=depth_m)
            tag = f"depth_{int(abs(depth_km))}km"
            set_axes(xyz=[True, False, True], ranges=self.axes_ranges, **kwargs)

        # Set the camera with top down view
        scene = self.engine.scenes[0]
        scene.scene.z_plus_view()

        # Plot extras
        colorscale(orientation="vertical", **kwargs)
        annotate(s=f"{tag.replace('_', ' ')}", c="k", width=0.2)

        # Finalize
        if save:
            mlab.savefig(save.format(tag=tag))
            return save.format(tag=tag)
        if show:
            mlab.show()
        return None

    def plot_depth_cross_section(self, choice, slice_at, show_axis=True,
                                 show=True, save=False, **kwargs):
        """
        Plot a depth cross-section of the model, with the Y-axis plotting depth

        If choice is X, the slice is normal the Y-axis
        if choice is Y, the slice is normal the X-axis

        kwargs passed to colorscale() and set_axes()

        :type choice: str
        :param choice: choice of axis to slice along, 'X', 'Y'
        :type slice_at: float
        :param slice_at: depth or distance value in meters to slice at
        :type show_axis: bool
        :param show_axis: show the axis of the plot
        :type save: bool
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        """
        choice = choice.upper()

        def slice_range(ranges_, slice_at_, choice_):
            """
            Convenience function to calculate where to slice by checking the
            min and max of the range and multiplying by a percentage

            :type ranges_: list
            :param ranges_: [xmin, xmax, ymin, ymax, zmin, zmax]
            :type slice_at_: float
            :param slice_at_: percentage of axis from [0,1]
            :type choice_: str
            :param choice_: choice of X or Y axis
            """
            if choice_ == "X":
                return slice_at_ * (ranges_[1] - ranges_[0]) + ranges_[0]
            elif choice_ == "Y":
                return slice_at_ * (ranges_[3] - ranges_[2]) + ranges_[2]

        # Determine what part of the axis we are slicing
        slice_val = slice_range(self.ranges, slice_at, choice)
        tag = slice_range(self.axes_ranges, slice_at, choice)
        logger.info(
                f"slicing '{self.fid}' across {choice} at {slice_val*1E-3}")
        if show_axis:
            set_axes(xyz=[True, True, False], ranges=self.axes_ranges,
                     **kwargs)
        self.cut_plane(choice=choice, slice_at=slice_val)

        # Differentiate which way were slicing to put the sources and receivers
        slice_x, slice_y = None, None
        if choice == "X":
            slice_x = slice_val
        elif choice == "Y":
            slice_y = slice_val

        # Plot the sources and receivers along a given slice
        if self.rcvs is not None:
            srcrcv(self.rcvs, color="k", x_value=slice_x, y_value=slice_y,
                   marker="2ddiamond")
        if self.srcs is not None:
            # Sources need to be spheres because I don't want to figure out how
            # to rotate 2d glyphs lol
            srcrcv(self.srcs, color="g", x_value=slice_x, y_value=slice_y,
                   marker="sphere")

        # Set the camera with a side on view. No preset so set to the preset and
        # rotate to get to the proper frame of reference
        scene = self.engine.scenes[0]
        if choice == "X":
            scene.scene.x_plus_view()
        elif choice == "Y":
            scene.scene.y_plus_view()
            scene.scene.camera.roll(90)

        # Plot extras
        colorscale(orientation="horizontal", **kwargs)
        # Annotation location changes based on which axis you slice
        anno_tag = f"{choice.lower()}={slice_at*1E2}%={tag:.2f}km"
        if choice == "X":
            annotate(x=.2, y=.3, s=anno_tag, c="k", width=0.15)
        elif choice == "Y":
            annotate(x=.27, y=.3, s=anno_tag, c="k", width=0.15)

        # Finalize
        if save:
            save_tag = save.format(tag=f"{tag:.0f}km")
            mlab.savefig(save_tag)
            return save_tag
        if show:
            mlab.show()
        return None

