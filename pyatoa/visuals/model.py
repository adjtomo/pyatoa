"""
Tools to use the Python package mayavi to interact with .vtk files.
Allows for creation of standardized looking plots of .vtk files for viewing
models and model slices etc.
"""
import os
from numpy import array
from mayavi import mlab
from mayavi.modules.surface import Surface
from mayavi.modules.scalar_cut_plane import ScalarCutPlane

from pyatoa.visuals.model_tools import (logger, get_coordinates, set_axes,
                                        annotate, coastline, colorscale,
                                        srcrcv, colors)


class Model:
    """
    A class to plot models from .vtk files in an automated fashion
    """
    def __init__(self, fid, srcs=None, rcvs=None, coast=None,
                 figsize=(1000, 1000),):
        """
        :param fid:
        :param srcs:
        :param rcvs:
        :param coast:
        """
        self.coords = get_coordinates(fid)

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
        Open a VTK file and return the engine that is visualizing it
        """
        if self.fig and self.engine and self.vtkfr:
            return
        logger.debug(f"Reading {self.fid} and creating figure size "
                     f"{self.figsize}")
        # Instantiate mlab
        self.fig = mlab.figure(size=self.figsize)
        self.engine = mlab.get_engine()
        self.engine.scenes[0].scene.background = colors["w"]
        self.vtkfr = self.engine.open(self.fid)

    def show_surface(self):
        """
        Show the surface of the vtk file from the top down
        """
        surface = Surface()
        self.engine.add_filter(surface, self.vtkfr)

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
        """
        coastline_z = self.coords[:, 2].max()
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

        # Plot the model at the surface or at some depth
        if depth_km:
            self.cut_plane(choice="Z", slice_at=depth_m)
            tag = f"depth_{int(abs(depth_km))}km"
        else:
            self.show_surface()
            tag = "surface"

        # Set the camera with top down view
        scene = self.engine.scenes[0]
        scene.scene.z_plus_view()

        # Plot extras
        colorscale(orientation="vertical", **kwargs)
        # set_axes(xyz=[True, True, False], **kwargs)
        annotate(s=f"{tag.replace('_', ' ')}", c="k", width=0.2)

        # Finalize
        if save:
            mlab.savefig(save.format(tag=tag))
        if show:
            mlab.show()

    @mlab.show
    def plot_depth_cross_section(self, choice, slice_at, axes=True, save=False,
                                 show=True, **kwargs):
        """
        Plot a depth cross-section of the model, with the Y-axis plotting depth

        If choice is X, the slice is parallel to Y-axis for a constant value of X
        same for choice of Y

        :type choice:
        :param choice:
        :type ratio:
        :param ratio:
        :type save:
        :param save:
        :type show:
        :param show:
        """
        # Generate axes for the data, Z values are now on the Y-axis so
        # visibility toggling needs to reflect that
        if axes:
            if choice == "X":
                set_axes(xyz=[True, True, False], **kwargs)
            elif choice == "Y":
                set_axes(xyz=[True, True, True], **kwargs)

        self.cut_plane(choice=choice, slice_at=slice_at)

        # Determine which direction we are slicing
        if self.rcvs or self.srcs:
            if choice == "X":
                slice_x, slice_y = slice_at, None
            elif choice == "Y":
                slice_x, slice_y = None, slice_at
        # Plot the sources and receivers along a given slice
        if self.rcvs:
            srcrcv(self.rcvs, color="w", x_value=slice_x, y_value=slice_y,
                   marker="2ddiamond")
        if self.srcs:
            srcrcv(self.srcs, color="g", x_value=slice_x, y_value=slice_y,
                   marker="2dcircle")

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
        tag = f"{choice.lower()}={slice_at*1E2:.2f}%"
        annotate(x=.2, y=.675, s=tag, c="k", width=0.1)

        # Finalize
        if save:
            mlab.savefig(save.format(tag=tag))
        if show:
            mlab.show()

    def top_down(self, depths=[None, 5, 10, 15, 25, 50], show=True, save=False):
        """
        Main function to visualize model from the surface down to a given depth
        """
        # Figure parameters
        save_fid = None
        if save:
            save_fid = self.fid + "_Z_{tag}.png"

        # Axes parameters
        font_factor = 1.1

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
                                    font_factor=font_factor,
                                    num_clabels=number_of_labels,
                                    round_to=round_to, save=save_fid, show=show)

    def depth_slice(self, xvalues=[0.5], yvalues=[0.5], show=True,
                      save=False):
        """
        Plot depth cross sections of the model parallel to X and Y axes
        """
        # ID's for files to plot
        save_fid_x, save_fid_y = None, None
        if save:
            save_fid_x = self.fid + "_X_{tag}.png"
            save_fid_y = self.fid + "_Y_{tag}.png"

        # Axes parameters
        font_factor = 1.1

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
        for xval in xvalues:
            self.startup()
            self.plot_depth_cross_section(choice="X", slice_at=xval,
                                          font_factor=font_factor,
                                          cmap=colormap, reverse=reverse,
                                          default_range=default_range,
                                          min_max=min_max, title=cbar_title,
                                          num_colors=number_of_colors,
                                          num_clabels=number_of_labels,
                                          round_to=round_to, save=save_fid_x,
                                          show=show
                                          )

        for yval in yvalues:
            self.startup()
            self.plot_depth_cross_section(choice="Y", slice_at=yval,
                                          font_factor=font_factor,
                                          cmap=colormap, reverse=reverse,
                                          default_range=default_range,
                                          min_max=min_max, title=cbar_title,
                                          num_colors=number_of_colors,
                                          num_clabels=number_of_labels,
                                          round_to=round_to, save=save_fid_y,
                                          show=show)




