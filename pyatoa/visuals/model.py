"""
A class used to plot model visualizations of .vtk files in an automatable
fashion, using the Mayavi library. Standard figure templates include horizontal
and vertical cross sections which show depth slices of the model.
Repeatedly used auxiliary functionality contained in model_tools
"""
import os
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
                 figsize=(1000, 1000), zero_origin=False, convert=1):
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
        self.srcs = srcs
        self.rcvs = rcvs
        self.coast = coast

        # Initiate the plotting engine
        self.figsize = figsize
        self.fig, self.engine, self.vtkfr = None, None, None

    def startup(self):
        """
        Open a VTK file and return the engine that is visualizing it.
        Whenever the figure is shown, startup needs to be called again,
        similar to how showing a matplotlib figure will destroy the instance.
        """
        logger.debug(f"Reading {self.fid} and creating figure size "
                     f"{self.figsize}")
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
            cut.implicit_plane.origin = array([origin[0], origin[1], origin[2]])
            cut.implicit_plane.normal = array([1, 0, 0])
        elif choice == "Y":
            cut.implicit_plane.origin = array([origin[0], origin[1], origin[2]])
            cut.implicit_plane.normal = array([0, 1, 0])
        elif choice == "Z":
            cut.implicit_plane.origin = array([origin[0], origin[1], slice_at])
            cut.implicit_plane.normal = array([0, 0, 1])

        cut.implicit_plane.widget.enabled = False

        return cut

    @mlab.show
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
        coastline_z = self.ranges[-1]
        if depth_km:
            depth_m = -1 * abs(depth_km) * 1E3
            coastline_z = depth_m

        # Put the coastline at some height above the topography, or at depth
        if self.coast:
            coastline(self.coast, coastline_z)

        # Plot the stations and receivers
        if self.rcvs:
            srcrcv(self.rcvs, color="w", z_value=coastline_z,
                   marker="2ddiamond")
        if self.srcs:
            srcrcv(self.srcs, color="g", z_value=coastline_z,
                   marker="2dcircle")

        # Plot the model at the given height
        # Axes behave weirdly when cutting planes so turn off Y-axis, not sure
        # why this works... sorry
        if depth_km:
            # Show a slice (plane) at a given depth
            self.cut_plane(choice="Z", slice_at=depth_m)
            tag = f"depth_{int(abs(depth_km))}km"
            set_axes(xyz=[True, False, True], ranges=self.axes_ranges, **kwargs)
        else:
            # Show a surface projection
            self.engine.add_filter(Surface(), self.vtkfr)
            tag = "surface"
            set_axes(xyz=[True, True, False], ranges=self.axes_ranges, **kwargs)

        # Set the camera with top down view
        scene = self.engine.scenes[0]
        scene.scene.z_plus_view()

        # Plot extras
        colorscale(orientation="vertical", **kwargs)
        annotate(s=f"{tag.replace('_', ' ')}", c="k", width=0.2)

        # Finalize
        if save:
            mlab.savefig(save.format(tag=tag))
        if show:
            mlab.show()

    @mlab.show
    def plot_depth_cross_section(self, choice, slice_at, show_axis=True,
                                 show=True, save=False, **kwargs):
        """
        Plot a depth cross-section of the model, with the Y-axis plotting depth

        If choice is X, the slice is parallel to Y-axis for a constant value of
        X same for choice of Y

        kwargs passed to colorscale()

        :type choice: str
        :param choice: choice of axis to slice along, 'X', 'Y' or 'Z'
        :type slice_at: float
        :param slice_at: depth or distance value in meters to slice at
        :type show_axis: bool
        :param show_axis: show the axis of the plot
        :type save: bool
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        """
        # Generate axes for the data, Z values are now on the Y-axis so
        # visibility toggling needs to reflect that
        if show_axis:
            if choice == "X":
                set_axes(xyz=[True, True, False], ranges=self.axes_ranges,
                         **kwargs)
            elif choice == "Y":
                set_axes(xyz=[True, True, True], ranges=self.axes_ranges,
                         **kwargs)

        self.cut_plane(choice=choice, slice_at=slice_at)

        # Determine which direction we are slicing
        if self.rcvs or self.srcs:
            slice_x, slice_y = None, None
            if choice == "X":
                # Determine where to slice based on the range and the percentage
                # given in the variable slice_at
                slice_x = (slice_at * (self.ranges[1] - self.ranges[0]) +
                           self.ranges[0])
                logger.info(f"Slicing X at {slice_x}")
            elif choice == "Y":
                slice_y = (slice_at * (self.ranges[3] - self.ranges[2]) +
                           self.ranges[2])
                logger.info(f"Slicing Y at {slice_y}")
        # Plot the sources and receivers along a given slice
        if self.rcvs:
            srcrcv(self.rcvs, color="k", x_value=slice_x, y_value=slice_y,
                   marker="2ddiamond")
        if self.srcs:
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
        if choice == "X":
            tag = f"{choice.lower()}={slice_x*self.convert:.2f}km"
            annotate(x=.2, y=.3, s=tag, c="k", width=0.15)
        elif choice == "Y":
            tag = f"{choice.lower()}={slice_y*self.convert:.2f}km"
            annotate(x=.27, y=.3, s=tag, c="k", width=0.15)

        # Finalize
        if save:
            mlab.savefig(save.format(tag=tag))
        if show:
            mlab.show()

    def top_down(self, depths=[None, 5, 10, 15, 25, 50], show=True, save=False):
        """
        Main function to visualize model from the surface down to a given depth

        :type depths: list of float
        :param depths: iterable list of depths to show the model at
        :type save: bool
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        """
        # Figure parameters
        save_fid = None
        if save:
            save_fid = self.fid + "_Z_{tag}.png"

        # Axes parameters
        font_factor = 1.1
        convert = 1E-3
        zero_origin = True

        # Colormap parameters
        colormap = "RdYlBu"
        reverse = False
        cbar_title = "log(Vs)"
        min_max = [-.25, .25]
        number_of_colors = 51
        number_of_labels = 5
        default_range = False
        round_to = 50  # only if default_range == True

        # Plot for a given set of depths
        for depth in depths:
            self.startup()
            self.plot_model_topdown(depth_km=depth, cmap=colormap,
                                    reverse=reverse, title=cbar_title,
                                    num_colors=number_of_colors,
                                    min_max=min_max, default_range=default_range,
                                    convert=convert, zero_origin=zero_origin,
                                    font_factor=font_factor,
                                    num_clabels=number_of_labels,
                                    round_to=round_to, save=save_fid, show=show)

    def depth_slice(self, x_values=[0.5], y_values=[0.5], show=True,
                    save=False):
        """
        Plot depth cross sections of the model parallel to X and Y axes

        :type x_values: list of float
        :param x_values: iterable list of X values to slice the model at
        :type y_values: list of int or None
        :param y_values: iterable list of Y values to slice the model at
        :type save: bool
        :param save: save the figure with a unique generic identifier
        :type show: bool
        :param show: show the figure after making it
        """
        # ID's for files to plot
        save_fid_x, save_fid_y = None, None
        if save:
            save_fid_x = self.fid + "_X_{tag}.png"
            save_fid_y = self.fid + "_Y_{tag}.png"

        # Axes parameters
        font_factor = 1.1
        convert = 1E-3
        zero_origin = True

        # Colormap parameters
        colormap = "RdYlBu"
        reverse = False
        cbar_title = "Vs (m/s)"
        number_of_colors = 51
        number_of_labels = 3
        default_range = False
        round_to = 50  # only if default_range == True
        min_max = [-.25, .25]

        # Initiate the engine and reader
        for xval in y_values:
            self.startup()
            self.plot_depth_cross_section(choice="X", slice_at=xval,
                                          font_factor=font_factor,
                                          convert=convert,
                                          zero_origin=zero_origin,
                                          cmap=colormap, reverse=reverse,
                                          default_range=default_range,
                                          min_max=min_max, title=cbar_title,
                                          num_colors=number_of_colors,
                                          num_clabels=number_of_labels,
                                          round_to=round_to, save=save_fid_x,
                                          show=show
                                          )

        for yval in x_values:
            self.startup()
            self.plot_depth_cross_section(choice="Y", slice_at=yval,
                                          font_factor=font_factor,
                                          convert=convert,
                                          zero_origin=zero_origin,
                                          cmap=colormap, reverse=reverse,
                                          default_range=default_range,
                                          min_max=min_max, title=cbar_title,
                                          num_colors=number_of_colors,
                                          num_clabels=number_of_labels,
                                          round_to=round_to, save=save_fid_y,
                                          show=show)



