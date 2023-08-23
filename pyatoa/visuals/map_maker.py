#!/usr/bin/env python3
"""
Map making functionality for the Pyatoa package. Creates Cartopy basemaps
featuring events and stations w/ optional moment tensors
"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from obspy.imaging.beachball import beach
from obspy.geodetics.flinnengdahl import FlinnEngdahl
from obspy.core.event.catalog import Catalog

from pyatoa import logger
from pyatoa.utils.form import format_event_name
from pyatoa.utils.srcrcv import gcd_and_baz
from pyatoa.utils.calculate import enforce_angle_pi

DEGREE_CHAR = u"\N{DEGREE SIGN}"


class MapMaker:
    """
    A class to call on the Basemap package to generate a map with
    source-receiver information
    """
    def __init__(self, cat, inv, dpi=100, figsize=None, figure=None,
                 gridspec=None, corners=None, corner_buffer_deg=2., **kwargs):
        """
        Initiate recurring parameters and parse out a few parameters for
        easier access
        """
        if isinstance(cat, Catalog):
            self.event = cat[0]
        else:
            self.event = cat

        self.inv = inv
        self.kwargs = kwargs

        # To be filled by initiate()
        self.fig = None
        self.ax = None

        # Set up a few useful parameters that will be called repeatedly
        self.ev_lat = self.event.preferred_origin().latitude
        self.ev_lon = self.event.preferred_origin().longitude
        self.sta_lat = inv[0][0].latitude
        self.sta_lon = inv[0][0].longitude

        # Used for coordinate transforms between lat/lon and projection
        self.ref_proj = ccrs.PlateCarree()

        # Extents is a tuple [lon_min, lon_max, lat_min, lat_max]
        self.extents = self.define_bounding_box(corners, corner_buffer_deg)

        if figsize is None:
            figsize = (600 / dpi, 600 / dpi)

        self.fig, self.ax, self.projection = self.initiate_figure(
            figsize, dpi, figure, gridspec
        )

    def define_bounding_box(self, corners=None, corner_buffer_deg=2):
        """
        Distribute the corners provided by the user, or determine corners
        using the event and station locations with a reasonable buffer.

        :type corners: dict
        :param corners: dict containing corner points, if None, lat lon values
            to be determiend by station and receiver locations
        :type buffer: float
        :param buffer: if no corners are given, put a buffer of length 'buffer'
            in units of degrees, around the min and max lat and lon values,
            to ensure that atleast some extra extent of map is covered.
            Defaults to 1 deg or roughly 111.11 km. But, if the distance covered
            between source and receiver is greater than 'buffer', than a quarter
            that distance will be used as the buffer. Confusing?
        :type corner_buffer_deg: float
        :param corner_buffer_deg: size of the bounding box to be generated
            around the source and receiver, units of degrees
        :rtype: tuple of float
        :preturn: [lon_min, lon_max, lat_min, lat_max]
        """
        # If no corners are given, provide a reasonable buffer around the
        # source and receiver locations.
        if corners is None:
            # Latitude is easy, just grab values directly
            lat_min = min(self.ev_lat, self.sta_lat)
            lat_max = max(self.ev_lat, self.sta_lat)

            d_lat = max(corner_buffer_deg, 0.25 * (lat_max - lat_min))
            # !!! BUG? Using __iadd__ (+=) kept returning NoneType not float
            # !!! Related to how ObsPy defining types (core/util/obspy_types)
            lat_min = lat_min - d_lat
            lat_max = lat_max + d_lat
            # lat_min -= d_lat
            # lat_max += d_lat  # <<< only lat_max was being affected

            # Crude scaling for longitude based on max latitude value to get a
            # roughly square buffer domain around the source and receiver
            if abs(lat_max) < 23:
                lon_scale = 1
            elif 23 <= abs(lat_max) < 45:
                lon_scale = 1.4
            else:
                lon_scale = 2.5

            lon_min = min(self.ev_lon, self.sta_lon)
            lon_max = max(self.ev_lon, self.sta_lon)
            d_lon = lon_scale * max(corner_buffer_deg,
                                    0.25 * (lon_max - lon_min)
                                    )
            lon_min = lon_min - d_lon
            lon_max = lon_max + d_lon
        else:
            # Parse the corners into usable values
            assert(isinstance(corners, dict))

            # Lower the corner keys to allow for upper and lower key entries
            corners = {key.lower(): val for key, val in corners.items()}

            lat_min = corners["lat_min"]
            lat_max = corners["lat_max"]
            lon_min = corners["lon_min"]
            lon_max = corners["lon_max"]

        lon_max = enforce_angle_pi(lon_max)
        lon_min = enforce_angle_pi(lon_min)

        return lon_min, lon_max, lat_min, lat_max

    def initiate_figure(self, figsize=None, dpi=100, figure=None, gridspec=None,
                        **kwargs):
        """
        Create a very barebones minimalist (black and white) map to plot data on

        .. note::
            kwargs passed to projection
            https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html
        """
        axis_linewidth = self.kwargs.get("axis_linewidth", 1.5)
        proj_str = self.kwargs.get("projection", "Stereographic")

        # Attempt to grab a plate projection from Cartopy
        try:
            projection = getattr(ccrs, proj_str)
        except AttributeError:
            logger.warning(f"{proj_str} is not a valid Cartopy CRS projection, "
                           f"setting default value: 'Stereographic'")
            projection = getattr(ccrs, "Stereographic")

        # Most major projections use the central_longitude function while only
        # some require central latitude, so we set that outside the init
        projection = projection(central_longitude=self.ev_lon, **kwargs)
        projection.central_latitude = self.ev_lat

        # Start fig either standalone or tandem with a waveform plot in GridSpec
        if figure is None:
            fig = plt.figure(figsize=figsize, dpi=dpi)
        else:
            fig = figure

        if gridspec is None:
            ax = plt.axes(projection=projection)
        else:
            ax = fig.add_subplot(gridspec[1], projection=projection)

        # Since the extents are in Lat/Lon, we set the extent in a lat/lon proj
        ax.set_extent(self.extents, crs=self.ref_proj)
        ax.coastlines(lw=axis_linewidth)
    
        gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                          y_inline=False, linewidth=.25, alpha=0.25, color="k")
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"rotation": 0}

        # Axis linewidth is set differently than in Matplotlib, see:
        ax.spines["geo"].set_linewidth(axis_linewidth)

        scale_bar(ax, 100, ref_proj=self.ref_proj)

        return fig, ax, projection

    def source(self, fm_type="focal_mechanism"):
        """
        Plot the source, either as a focal mechanism, moment tensor, or as a
        simple point, based on the input.

        .. note::
            scale_source kwarg was guessed with trial and error and is based on
            the guessed returned length from the scale_bar() function defined
            at the bottom, which tries to guess a reasonable length of the 
            scale bar based on the dimensions of the map

        :type fm_type: str
        :param fm_type: choice to plot
            focal_mechanism: 6 component focal mechanism
            strike_dip_rake: classic double couple look
        """
        marker = self.kwargs.get("source_marker", "o")
        color = self.kwargs.get("source_color", "indianred")
        size = self.kwargs.get("source_size", 75)

        lw = self.kwargs.get("source_lw", 1.75)
        scale_source = self.kwargs.get("scale_source", .3)

        # No focal mechanism? Just plot a marker
        self.ax.scatter(self.ev_lon, self.ev_lat, transform=self.ref_proj,
                        marker=marker, color=color, edgecolor="k", linewidth=lw,
                        s=size, zorder=10)

        if hasattr(self.event, "focal_mechanisms") and \
                self.event.focal_mechanisms:
            if fm_type == "focal_mechanism":
                fm = self.event.focal_mechanisms[0].moment_tensor.tensor or \
                     self.event.preferred_focal_mechanism().moment_tensor.tensor
                beach_input = [fm["m_rr"], fm["m_tt"], fm["m_pp"],
                               fm["m_rt"], fm["m_rp"], fm["m_tp"]
                               ]
            elif fm_type == "strike_dip_rake":
                nod_plane = self.event.focal_mechanisms[0].nodal_planes or \
                            self.event.preferred_focal_mechanism().nodal_planes
                # try determine the preferred nodal plane, default to 1
                try:
                    sdr = nod_plane[f"nodal_plane_{nod_plane.preferred_plane}"]
                except AttributeError:
                    sdr = nod_plane.nodal_plane_1
                beach_input = [sdr.strike, sdr.dip, sdr.rake]
            else:
                raise ValueError("fm_type must be 'focal_mechanism' or "
                                 "'strike_dip_rake")

            # Transform lon/lat to given projection to give to beachball
            x, y = self.projection.transform_point(
                x=self.ev_lon, y=self.ev_lat, src_crs=self.ref_proj
            )
            # Guess a width of the focal mechanism that fits nicely on the plot
            width = scale_source * scale_bar(ax=self.ax, length=None,
                                             return_length=True,
                                             ref_proj=self.ref_proj)

            b = beach(beach_input, xy=(x, y), width=width,
                      linewidth=lw, facecolor=color)
            b.set_zorder(11)
            self.ax.add_collection(b)

    def receiver(self):
        """
        Plot the receiver with a standard look
        """
        marker = self.kwargs.get("station_marker", "v")
        color = self.kwargs.get("station_color", "forestgreen")
        size = self.kwargs.get("station_size", 75)
        lw = self.kwargs.get("station_lw", 1.5)

        self.ax.scatter(self.sta_lon, self.sta_lat,
                        transform=self.ref_proj, marker=marker,
                        color=color, linewidth=lw, s=size,
                        edgecolor="k", zorder=10)

    def connect(self):
        """
        Plot a connecting line between source and receiver
        """
        ls = self.kwargs.get("srcrcv_linestyle", "--")
        lw = self.kwargs.get("srcrcv_linewidth", 1.5)
        lc = self.kwargs.get("srcrcv_color", "k")

        self.ax.plot([self.ev_lon, self.sta_lon], [self.ev_lat, self.sta_lat],
                     transform=self.ref_proj, linestyle=ls, linewidth=lw,
                     c=lc, zorder=8)

    def annotate(self, location="lower-right", anno_latlon=False):
        """
        Annotate event receiver information into bottom right corner of the map

        TODO figure out where to put definition of 'location'

        :type location: str
        :param location: location of the annotation block, available:
            'upper-right', 'lower-right', 'upper-left', 'lower-left', 'center'
        :type anno_latlon: bool
        :param anno_latlon: annotate the latitude and longitude values of the
            source and receiver next to their markers. Not always very clean
            so defaults to off.
        """
        fontsize = self.kwargs.get("anno_fontsize", 8)
        location = self.kwargs.get("anno_location", "lower-right")

        ok_locs = ["lower-right", "upper-right", "lower-left",
                   "upper-left", "center"]
        assert location in ok_locs, f"location must be in {ok_locs}"

        # Determine the location of the annotation
        if location == "center":
            x = y = 0.5
            ha = va = ma = "center"
        if "lower" in location:
            y = 0.01
            va = "bottom"
        elif "upper" in location:
            y = 0.97
            va = "top"
        if "left" in location:
            x = 0.05
            ha = ma = "left"
        elif "right" in location:
            x = 0.95
            ha = ma = "right"

        # Collect some useful information
        event_id = format_event_name(self.event)
        sta_id = f"{self.inv[0].code}.{self.inv[0][0].code}"
        gc_dist, baz = gcd_and_baz(self.event, self.inv[0][0])
        origin_time = self.event.origins[0].time
        depth = self.event.preferred_origin().depth * 1E-3
        if self.event.preferred_magnitude() is not None:
            magnitude = self.event.preferred_magnitude().mag
            mag_type = self.event.preferred_magnitude().magnitude_type
            mag_str = f"{mag_type} {magnitude:.2f}\n"
        else:
            mag_str = "Magnitude info unavailable\n"

        region = FlinnEngdahl().get_region(self.ev_lon, self.ev_lat)

        self.ax.text(s=(f"{region.title()}\n"
                        f"{'-'*len(region)}\n"
                        f"{event_id} / {sta_id}\n"
                        f"{origin_time.format_iris_web_service()}\n"
                        f"{mag_str}"
                        f"Depth: {depth:.2f} km\n"
                        f"Dist: {gc_dist:.2f} km\n"
                        f"BAz: {baz:.2f}{DEGREE_CHAR}\n"
                        ),
                     x=x, y=y, ha=ha, va=va, ma=ma, zorder=5,
                     fontsize=fontsize, transform=self.ax.transAxes
                     )

        if anno_latlon:
            # Annotate the lat lon values next to source and receiver
            self.ax.text(s=f"\t({self.ev_lat:.2f}, {self.ev_lon:.2f})",
                         x=self.ev_lon, y=self.ev_lat, fontsize=fontsize)
            self.ax.text(s=f"\t({self.sta_lat:.2f}, {self.sta_lon:.2f})",
                         x=self.sta_lon, y=self.sta_lat, fontsize=fontsize)

    def plot(self, show=True, save=None, **kwargs):
        """
        Main function to generate the basemap, plot all the components, and
        show or save the figure

        :type show: bool
        :param show: show the figure in the gui
        :type save: str
        :param save: if not None, save to the given path stored in this var.
        :type corners: dict
        :param corners: dict containing corner points, if None, lat lon values
            to be determiend by station and receiver locations
        """
        self.source()
        self.receiver()
        self.connect()
        self.annotate()
        
        if save:
            plt.savefig(save)
        if show:
            plt.show()


def scale_bar(ax, length=None, location=(0.85, 0.95), linewidth=3, 
              return_length=False, ref_proj=ccrs.PlateCarree()):
    """
    Create a scale bar on a Cartopy plot.
    Modifiedd from: https://stackoverflow.com/questions/32333870/
                          how-can-i-show-a-km-ruler-on-a-cartopy-matplotlib-plot

    :type ax: matplotlib.pyplot.axes
    :param ax: axes to draw the scalebar on.
    :type length: float
    :param length: length of the scalebar in km.
    :type location: tuple of floats
    :param location: center of the scalebar in axis coordinates.
        (ie. 0.5 is the middle of the plot)
    :type linewidth: float
    :param linewidth: the thickness of the scalebar.
    :type return_length: bool
    :param return_length: Simply returns the scaled length of the bar, added
        to use for scaling of the moment tensor
    """
    # Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ref_proj)

    # Make tmc horizontally centred on the middle of the map,
    # vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly, approx=True)

    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)

    # Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    # Calculate a scale bar length if none has been given
    # Theres probably a more pythonic way of rounding the number but this works
    if not length:
        length = (x1 - x0) / 5000  # in km
        ndim = int(np.floor(np.log10(length)))  # number of digits in number
        length = round(length, -ndim)  # round to 1sf

        def scale_number(x):
            """Returns numbers starting with the list"""
            if str(x)[0] in ["1", "2", "5"]:
                return int(x)
            else:
                return scale_number(x - 10 ** ndim)
        length = scale_number(length)

    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    
    if return_length:
        return bar_xs[1] - bar_xs[0]
    else:
        # Plot the scalebar
        ax.plot(bar_xs, [sby, sby], transform=tmc, color="k", 
                linewidth=linewidth)
        # Plot the scalebar label
        ax.text(sbx, sby, f"{length} km", transform=tmc,
                horizontalalignment="center", verticalalignment="bottom")
