"""
Functionality to plot regular xyz grid files
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt


class XYZViewer:
    """
    A class to read, manipulate and plot structured grid files (xyz) that are
    in the format of the Specfem external tomography file.
    """
    def __init__(self, variables=None):
        """
        Initiate the class with some empty variables
        """
        self.variables = variables or ["x", "y", "z", "vp", "vs", 
                                       "rho", "qp", "qs"]
        self.zero_origin = False
        self.xgrid = None
        self.ygrid = None
        self._coast = None

    def read(self, fid, fmt="specfem"):
        """
        A wrapper for np.loadtxt to read in xyz files.

        ..note::
            For now this function assumes that it is reading an external
            tomography file that is formatted how SPECFEM expects it.

        :type fid: str
        :param fid: file id to read from
        """
        if fmt == "specfem":
            values = np.loadtxt(fid, dtype=float, skiprows=4).T
        elif fmt == "semslicer":
            values = np.loadtxt(fid, dtype=float, delimiter=",").T
        else:
            raise ValueError("fmt must be 'specfem' or 'semslicer'")

        # It is possible that attenuation is not included, if so change the
        # default variables
        assert(len(values) == len(self.variables)), \
                (f"Number of columns {len(values)} does not match number of "
                 f"excepted variables {len(self.variables)}")

        # Set each of the columns as an internal variable
        for i, variable in enumerate(self.variables):
            setattr(self, variable, values[i])

        self.check()

    def read_semslicer(self, fid):
        """
        A wrapper for np.loadtxt to read in xyz files that are outputted by the
        semslicer fortran script

        :type fid: str
        :param fid: file id to read from
        """
        values = np.loadtxt(fid, dtype=float, skiprows=4).T

        # It is possible that attenuation is not included
        if len(values) != len(self.variables):
            variables = self.variables[:6]
            assert(len(values) == len(variables)), \
                f"Number of columns ({len(values)}) does not match known format"

        # Set each of the columns as an internal variable
        for i, variable in enumerate(self.variables):
            setattr(self, variable, values[i])

        self.check()

    def compare(self, fid, choice=""):
        """
        Read in a corresponding xyz file and difference values
        :param fid:
        :return:
        """
        pass

    def show(self, **kwargs):
        """
        A wrapper for pyplot.show() to avoid importing pyplot just for show()
        """
        plt.show(**kwargs)

    def savefig(self, fid):
        """
        A wrapper for pyplot.savefig()
        """
        if not os.path.exists(os.path.dirname(fid)):
            os.mkdir(os.path.dirname(fid))

        plt.savefig(fid)

    def close(self, *args, **kwargs):
        """
        A warpper for pyplot.close('all')
        """
        plt.close(*args, **kwargs)

    def grid(self):
        """
        Define a Numpy mesh grid that can be used for contour plotting, using
        internal variables. Will be used to reshape data arrays
        """
        if self.zero_origin:
            self.xgrid, self.ygrid = np.meshgrid(
                np.arange(0, self.x_max - self.x_min + self.dx, self.dx),
                np.arange(0, self.y_max - self.y_min + self.dy, self.dy)
            )
        else:
            self.xgrid, self.ygrid = np.meshgrid(
                np.arange(self.x_min, self.x_max + self.dx, self.dx),
                np.arange(self.y_min, self.y_max + self.dy, self.dy)
            )
        assert (np.shape(self.xgrid) == (self.ny, self.nx)), "Error in gridding"

    def coast(self, fid=None):
        """
        Plot the coastline based on a three colunn
        :return:
        """
        if self._coast is None:
            assert(fid is not None), f"2-column ascii required for coast"
            self._coast = np.loadtxt(fid)

        if self.zero_origin and self._coast[0][0] != 0:
            for i in [0, 1]:
                self._coast[:, i] = self._coast[:, i] - self._coast[:, i].min()

        plt.scatter(self._coast[:, 0], self._coast[:, 1], c="k", s=1, zorder=5)

    def check(self):
        """
        Check the min/max values, grid spacing, npts etc. assign internally
        """
        for variable in self.variables:
            values = getattr(self, variable)
            setattr(self, f"{variable}_min", values.min())
            setattr(self, f"{variable}_max", values.max())
            if variable in ["x", "y", "z"]:
                unique_values = np.unique(values)
                spacing = abs(unique_values[1] - unique_values[0])
                setattr(self, f"unique_{variable}", unique_values)
                setattr(self, f"n{variable}", len(unique_values))
                setattr(self, f"d{variable}", spacing)

        setattr(self, "npts", len(values))
        self.grid()

    def decimate(self, factor):
        """
        Crudely decimate the values, if, e.g. only a low-resolution plot needs
        to be made.

        :type factor: int
        :param factor: factor to decimate npts by
        """
        # Two methods for decimating, one if factor is a factor of npts
        if not self.npts % factor:
            for variable in self.variables:
                values = getattr(self, variable)
                values = values.reshape(-1, factor).max(1)
                setattr(self, variable, values)
        else:
            for variable in self.variables:
                values = getattr(self, variable)
                values = np.maximum.reduceat(values,
                                             np.arange(0, self.npts, factor)
                                             )
                setattr(self, variable, values)

        self.check()

    def depth_slice(self, depth_km, variable, **kwargs):
        """
        Plot a slice across the XY plane at a given depth value Z
        Kwargs are passed to pyplot.countourf

        :param depth_km:
        :return:
        """
        assert(depth_km in self.unique_z), f"{depth_km} is not a valid Z value"
        assert(variable in self.variables), \
            f"{variable} not in {self.variables}"

        idx = np.where(self.z == depth_km)[0]
        value = getattr(self, variable)[idx]

        value = value.reshape(np.shape(self.xgrid))
        assert(idx.any()), "No values found for Z={depth_km}km"

        f, ax = plt.subplots()
        plt.contourf(self.xgrid, self.ygrid, value, **kwargs)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(f"{variable.title()}", rotation=270, labelpad=15)

        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
        ax.set_aspect("equal")

        plt.title(f"{variable.title()} at Z={depth_km} km")

        return f, ax

    def cross_section(self, choice, value):
        """
        Plot a cross section across the XZ or YZ plane given an X or Y value

        :param choice:
        :param value:
        :return:
        """
        assert(choice.lower() in ["x", "y"]), "Choice must be 'x' or 'y'"
        values = getattr(self, "{choice}_unique")
        assert(value in values), f"Value {value} not in {values}"

        idx = np.where(values == value)[0]

        f, ax = plt.subplots()
        plt.contourf()

    def volume(self, variable):
        """
        Plot the 3D volume
        :return:
        """
        pass

        assert(variable in self.variables), f"{variable} not found"
        f = plt.figure()
        ax = plt.axes(projection="3d")
        ax.plot_surface(self.x, self.y, self.z, color=getattr(self, variable))


if __name__ == "__main__":
    from glob import glob

    xyz = XYZViewer()
    xyz.zero_origin = True
    coast = ("/Users/Chow/Documents/academic/vuw/data/carto/coastline/"
             "coast_nznorth_utm60.txt")
    # x, y, z = interface = np.loadtxt(
    #     "williams_hikurangi_interface_utm60_nonan.xyz").T
    for fid in sys.argv[1:]:
        xyz.read(fid, fmt="specfem")
        # xyz.decimate(2)
        for depth in xyz.unique_z:
            xyz.depth_slice(depth, "vs", cmap="RdYlBu", levels=21)
            xyz.coast(coast)
            savefid = f"{os.path.splitext(fid)[0]}_{depth}m.png"
            xyz.savefig(os.path.join("figures", savefid))
            xyz.close()


