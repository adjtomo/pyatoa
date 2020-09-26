"""
Functionality to plot regular xyz grid files
"""
import numpy as np
import matplotlib.pyplot as plt


class GridMaker:
    """
    A class to read, manipulate and plot structured grid files (xyz) that are
    in the format of the Specfem external tomography file.
    """
    def __init__(self):
        """
        Initiate the class with some empty variables
        """
        self.variables = ["x", "y", "z", "vp", "vs", "rho", "qp", "qs"]
        self.zero_origin = False
        self.xgrid = None
        self.ygrid = None
        self._coast = None

    def read(self, fid):
        """
        A wrapper for np.loadtxt to read in xyz files.

        ..note::
            For now this function assumes that it is reading an external
            tomography file that is formatted how SPECFEM expects it.

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
        plt.savefig(fid)

    def close(self):
        """
        A warpper for pyplot.close('all')
        """
        plt.close("all")

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

        plt.title(f"{variable.title()} at Z={depth_km} km")

        return f, ax

    # def volume(self, variable):
    #     """
    #     Plot the 3D volume
    #     :return:
    #     """
    #     assert(variable in self.variables), f"{variable} not found"
    #     f = plt.figure()
    #     ax = plt.axes(projection="3d")
    #     ax.plot_surface(self.x, self.y, self.z, color=getattr(self, variable))



