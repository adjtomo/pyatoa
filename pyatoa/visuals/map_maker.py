#!/usr/bin/env python3
"""
Create Basemaps featuring stations and moment tensors.

Note:
    The Basemap package is deprecated so these functionalities may have to be
    replaced, e.g. with Cartopy, at some point.

Functions used to produce maps using Basemap that have a standard look across
the Pyatoa workflow.
"""
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from obspy.geodetics.flinnengdahl import FlinnEngdahl
from obspy.core.event.catalog import Catalog

from pyatoa.utils.srcrcv import gcd_and_baz


class MapMaker:
    """
    A class to call on the Basemap package to generate a map with
    source-receiver information
    """
    def __init__(self, cat, inv, **kwargs):
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

        # To be filled by plot()
        self.fig = None
        self.ax = None
        self.m = None

        # To be filled in by initiate()
        self.ev_x = None
        self.ev_y = None
        self.sta_x = None
        self.sta_y = None

        # Set up a few useful parameters that will be called repeatedly
        self.ev_lat = self.event.preferred_origin().latitude
        self.ev_lon = self.event.preferred_origin().longitude
        self.sta_lat = inv[0][0][0].latitude
        self.sta_lon = inv[0][0][0].longitude

        # To be filled in by plot()
        self.lat_min = None
        self.lat_max = None
        self.lon_min = None
        self.lon_max = None

    def check_corners(self, corners=None):
        """
        Distribute the corners provided by the user, or determine corners
        using the event and station locations

        :type corners: dict
        :param corners: dict containing corner points, if None, lat lon values
            to be determiend by station and receiver locations
        """
        if corners is None:
            # If no corners are given, provide a reasonable buffer around the
            # source and receiver locations
            lon_min = min(self.ev_lon, self.sta_lon)
            lon_max = max(self.ev_lon, self.sta_lon)
            d_lon = 0.5 * (lon_max - lon_min)
            self.lon_min = lon_min - d_lon
            self.lon_max = lon_max + d_lon

            lat_min = min(self.ev_lat, self.sta_lat)
            lat_max = max(self.ev_lat, self.sta_lat)
            d_lat = 0.5 * (lat_max - lat_min)
            self.lat_min = lat_min - d_lat
            self.lat_max = lat_max + d_lat
        else:
            # Parse the corners into usable values
            assert(isinstance(corners, dict))
            self.lat_min = corners["lat_min"]
            self.lat_max = corners["lat_max"]
            self.lon_min = corners["lon_min"]
            self.lon_max = corners["lon_max"]

    def initiate(self, dpi, figsize):
        """
        Set up the basemap object with a certain defined look
        """
        # Optional mpl kwargs to allow placing the map inside another figure
        figure = self.kwargs.get("figure", None)
        ax = self.kwargs.get("ax", None)

        # Basemap kwargs
        area_thresh = self.kwargs.get("area_thresh", None)
        continent_color = self.kwargs.get("contininent_color", "w")
        lake_color = self.kwargs.get("lake_color", "w")
        coastline_zorder = self.kwargs.get("coastline_zorder", 5)
        coastline_linewidth = self.kwargs.get("coastline_linewidth", 2.0)
        axis_linewidth = self.kwargs.get("axis_linewidth", 2.0)
        axis_fontsize = self.kwargs.get("axis_fontsize", 8)
        fill_color = self.kwargs.get("fill_color", "w")
        plw = self.kwargs.get("parallel_linewidth", 0.)
        mlw = self.kwargs.get("meridian_linewidth", 0.)
        projection = self.kwargs.get("projection", "stere")
        resolution = self.kwargs.get("resolution", "l")

        # Initiate matplotlib instances
        if figure is None:
            self.fig = plt.figure(figsize=figsize, dpi=dpi)
        else:
            self.fig = figure

        # Initiate map and draw in style
        self.m = Basemap(projection=projection, resolution=resolution,
                         rsphere=6371200,
                         lat_0=(self.lat_min + self.lat_max)/2,
                         lon_0=(self.lon_min + self.lon_max)/2,
                         llcrnrlat=self.lat_min, urcrnrlat=self.lat_max,
                         llcrnrlon=self.lon_min, urcrnrlon=self.lon_max,
                         area_thresh=area_thresh, ax=ax,
                         )

        # By default, no meridan or parallel lines
        self.m.drawparallels(np.arange(int(self.lat_min), int(self.lat_max), 1),
                             labels=[1, 0, 0, 0], linewidth=plw,
                             fontsize=axis_fontsize
                             )
        self.m.drawmeridians(
            np.arange(int(self.lon_min), int(self.lon_max) + 1, 1),
            labels=[0, 0, 0, 1], linewidth=mlw, fontsize=axis_fontsize
        )

        # Create auxiliary parts of the map for clarity
        self.m.drawcoastlines(linewidth=coastline_linewidth,
                              zorder=coastline_zorder)
        self.m.fillcontinents(color=continent_color, lake_color=lake_color)
        self.m.drawmapboundary(fill_color=fill_color)
        self.scalebar()
        for axis in ["top", "bottom", "left", "right"]:
            plt.gca().spines[axis].set_linewidth(axis_linewidth)

        # Calculate the source-receiver locations based on map coordinates
        self.ev_x, self.ev_y = self.m(self.ev_lon, self.ev_lat)
        self.sta_x, self.sta_y = self.m(self.sta_lon, self.sta_lat)

    def scalebar(self):
        """
        Put the scale bar in a corner at a reasonable distance from each edge
        """
        loc = self.kwargs.get("scalebar_location", "upper-right")
        fontsize = self.kwargs.get("scalebar_fontsize", 8)
        lw = self.kwargs.get("scalebar_linewidth", 1)

        if loc == "upper-right":
            lat_pct = 0.94
            lon_pct = 0.875
        elif loc == "lower-right":
            lat_pct = 0.04
            lon_pct = 0.9

        # Place the scalebar based on the the corners given in lat/lon
        lat = self.lat_min + (self.lat_max - self.lat_min) * lat_pct
        lon = self.lon_min + (self.lon_max - self.lon_min) * lon_pct

        self.m.drawmapscale(lon, lat, lon, lat, 100,
                            yoffset=0.01 * (self.m.ymax - self.m.ymin),
                            zorder=100, linewidth=lw, fontsize=fontsize
                            )

    def source(self, fm_type="focal_mechanism"):
        """
        Plot the source, either as a focal mechanism, moment tensor, or as a
        simple point, based on the input.

        :type fm_type: str
        :param fm_type: choice to plot
            focal_mechanism: 6 component focal mechanism
            strike_dip_rake: classic double couple look
        """
        marker = self.kwargs.get("source_marker", "o")
        color = self.kwargs.get("source_color", "indianred")
        lw = self.kwargs.get("source_lw", 1.75)
        width = self.kwargs.get("source_width", 35)

        # No focal mechanism? Just plot a marker
        self.m.scatter(self.ev_x, self.ev_y, marker=marker,
                       color=color, edgecolor="k", linewidth=lw)

        if hasattr(self.event, "focal_mechanisms"):
            if fm_type == "focal_mechanism":
                fm = self.event.focal_mechanisms[0].moment_tensor.tensor or \
                     self.event.preferred_focal_mechanism().moment_tensor.tensor
                beach_input = [fm['m_rr'], fm['m_tt'], fm['m_pp'],
                               fm['m_rt'], fm['m_rp'], fm['m_tp']
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

            b = beach(beach_input, xy=(self.ev_x, self.ev_y), width=width,
                      linewidth=lw, facecolor=color, axes=plt.gca())
            b.set_zorder(10)
            plt.gca().add_collection(b)

    def receiver(self):
        """
        Plot the receiver with a standard look
        """
        marker = self.kwargs.get("station_marker", "v")
        color = self.kwargs.get("station_color", "forestgreen")
        size = self.kwargs.get("station_size", 90)
        lw = self.kwargs.get("station_lw", 1.5)

        self.m.scatter(self.sta_x, self.sta_y, marker=marker, color=color,
                       linewidth=lw, s=size, edgecolor="k", zorder=10)

    def connect(self):
        """
        Plot a connecting line between source and receiver
        """
        ls = self.kwargs.get("srcrcv_linestyle", "--")
        lw = self.kwargs.get("srcrcv_linewidth", 1.5)
        lc = self.kwargs.get("srcrcv_color", "k")

        self.m.plot([self.ev_x, self.sta_x], [self.ev_y, self.sta_y],
                    linestyle=ls, linewidth=lw, c=lc, zorder=8)

    def annotate(self, location="lower-right", anno_latlon=False):
        """
        Annotate event receiver information into bottom right corner of the map

        :type location: str
        :param location: location of the annotation block, available:
            'upper-right', 'lower-right', 'upper-left', 'lower-left', 'center'
        :type anno_latlon: bool
        :param anno_latlon: annotate the latitude and longitude values of the
            source and receiver next to their markers. Not always very clean
            so defaults to off.
        """
        fontsize = self.kwargs.get("anno_fontsize", 8)
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
        event_id = self.event.resource_id.id.split('/')[1]
        sta_id = f"{self.inv[0].code}.{self.inv[0][0].code}"
        gc_dist, baz = gcd_and_baz(self.event, self.inv[0][0])
        origin_time = self.event.origins[0].time
        depth = self.event.preferred_origin().depth * 1E-3
        magnitude = self.event.preferred_magnitude().mag
        mag_type = self.event.preferred_magnitude().magnitude_type
        region = FlinnEngdahl().get_region(self.ev_lon, self.ev_lat)

        # Need to use plot because basemap object has no annotate method
        plt.gca().text(s=(f"{region.title()}\n"
                          f"{'-'*len(region)}\n"
                          f"{event_id} / {sta_id}\n"
                          f"{origin_time.format_iris_web_service()}\n"
                          f"{mag_type} {magnitude:.2f}\n"
                          f"Depth: {depth:.2f} km\n"
                          f"Dist: {gc_dist:.2f} km\n"
                          f"BAz: {baz:.2f} deg\n"
                          ),
                       x=x, y=y, ha=ha, va=va, ma=ma,
                       transform=plt.gca().transAxes, zorder=5,
                       fontsize=fontsize,
                       )

        if anno_latlon:
            # Annotate the lat lon values next to source and receiver
            plt.gca().text(s=f"\t({self.ev_lat:.2f}, {self.ev_lon:.2f})",
                           x=self.ev_x, y=self.ev_y, fontsize=fontsize)
            plt.gca().text(s=f"\t({self.sta_lat:.2f}, {self.sta_lon:.2f})",
                           x=self.sta_x, y=self.sta_y, fontsize=fontsize)

    def plot(self, show=True, save=None, corners=None, **kwargs):
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
        # Allow kwarg updating in the plot call
        self.kwargs.update(kwargs)

        # Matplotlib kwargs
        dpi = self.kwargs.get("dpi", 100)
        figsize = self.kwargs.get("figsize", (600 / dpi, 600 / dpi))
        location = self.kwargs.get("anno_location", "lower-right")

        self.check_corners(corners)
        self.initiate(dpi, figsize)
        self.source()
        self.receiver()
        self.connect()
        self.annotate(location)

        if show:
            plt.show()
        if save:
            plt.savefig(save, dpi=dpi, figsize=figsize)


