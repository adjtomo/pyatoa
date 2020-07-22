#!/usr/bin/env python3
"""
Create a simple source-receiver plot, providing geographical information
relating the source and receiver. Meant to replace the map making functionality
originally written into Pyatoa to avoid the deprecated Basemap package.
"""
import matplotlib.pyplot as plt

from pyatoa.utils.srcrcv import gcd_and_baz
from pyatoa.visuals.manager_plotter import pretty_grids

from obspy.core.event.catalog import Catalog
from obspy.imaging.beachball import beach
from obspy.geodetics.flinnengdahl import FlinnEngdahl


class SourceReceiver:
    """
    A class controlling source-receiver plotting
    """
    def __init__(self, cat, inv, **kwargs):
        """
        Source receiver plot only requires an event an station object

        :type event: obspy.core.event.Event or obspy.core.event.Catalog
        :param event: a catalog or event object which may or may not contain
            focal mechanism or strike dip rake information.
        :type station: obspy.core.inventory.Inventory
        :param station: an inventory containing relevant network and station
            information. Should only contain a single station as the class
            will only look at the zeroth component
        :type kwargs: dict
        :param kwargs: list of optional key word arguments that can be
            passed in to change the specific look of the figure
        """
        if isinstance(cat, Catalog):
            self.event = cat[0]
        else:
            self.event = cat

        self.inv = inv
        self.kwargs = kwargs
        self.fig = None
        self.ax = None

        # Parse out the event and station information
        self.ev_lat = self.event.preferred_origin().latitude
        self.ev_lon = self.event.preferred_origin().longitude

        self.sta_lat = inv[0][0][0].latitude
        self.sta_lon = inv[0][0][0].longitude

    def setup(self):
        """
        Create the plot and assign the figure and axes to internal variables
        :return:
        """
        figsize = self.kwargs.get("figsize", (5, 4))
        dpi = self.kwargs.get("dpi", 100)
        fontsize = self.kwargs.get("fontsize", 8)
        axes_linewidth = self.kwargs.get("axes_linewidth", 1)
        buffer = self.kwargs.get("buffer", 0.1)

        self.fig, self.ax = plt.subplots(figsize=figsize, dpi=dpi)
        pretty_grids(self.ax, fontsize=fontsize, linewidth=axes_linewidth,
                     sci_format=False)

        # Push the bounds slightly outside the source and receiver locations
        x_min = min(self.ev_lon, self.sta_lon)
        x_max = max(self.ev_lon, self.sta_lon)
        dx = buffer * (x_max - x_min)
        self.ax.set_xlim([x_min - dx, x_max + dx])

        y_min = min(self.ev_lat, self.sta_lat)
        y_max = max(self.ev_lat, self.sta_lat)
        dy = buffer * (y_max - y_min)
        self.ax.set_ylim([y_min - dy, y_max + dy])

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
        color = self.kwargs.get("source_color", "r")
        lw = self.kwargs.get("source_lw", 1.75)
        width = self.kwargs.get("source_width", 60)

        # No focal mechanism? Just plot a marker
        self.ax.scatter(self.ev_lon, self.ev_lat, marker=marker,
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

            b = beach(beach_input, xy=(self.ev_lon, self.ev_lat), width=width,
                      linewidth=lw, facecolor=color, axes=self.ax)
            b.set_zorder(10)
            self.ax.add_collection(b)

    def receiver(self):
        """
        Plot the receiver with a standard look
        """
        marker = self.kwargs.get("station_marker", "v")
        color = self.kwargs.get("station_color", "w")
        size = self.kwargs.get("station_size", 80)
        lw = self.kwargs.get("station_lw", 1.25)

        self.ax.scatter(self.sta_lon, self.sta_lat, marker=marker, color=color,
                        linewidth=lw, s=size, edgecolor="k", zorder=10)

    def connecting_line(self):
        """
        Plot a connecting line between source and receiver
        """
        ls = self.kwargs.get("srcrcv_linestyle", "--")
        lw = self.kwargs.get("srcrcv_linewidth", 1.5)
        lc = self.kwargs.get("srcrcv_color", "k")

        self.ax.plot([self.ev_lon, self.sta_lon], [self.ev_lat, self.sta_lat],
                     linestyle=ls, linewidth=lw, c=lc, zorder=8)

    def annotations(self):
        """
        Provide annotations
        :return:
        """
        # Collect some useful information
        event_id = self.event.resource_id.id.split('/')[1]
        sta_id = f"{self.inv[0].code}.{self.inv[0][0].code}"
        gc_dist, baz = gcd_and_baz(self.event, self.inv[0][0])
        origin_time = self.event.origins[0].time
        depth = self.event.preferred_origin().depth * 1E-3
        magnitude = self.event.preferred_magnitude().mag
        mag_type = self.event.preferred_magnitude().magnitude_type
        region = FlinnEngdahl().get_region(self.ev_lon, self.ev_lat)

        # Find the mid-point of the line for annotation
        # x_mid = (self.ev_lon + self.sta_lon) / 2
        # y_mid = (self.ev_lat + self.sta_lat) / 2

        self.ax.text(s=(f"{region.title()}\n"
                        f"Dist: {gc_dist:.2f} km / BAz: {baz:.2f}\n"
                        f"{'-'*10}\n"
                        f"{event_id} ({self.ev_lat:.2f}, {self.ev_lon:.2f})\n"
                        f"   Start: {origin_time.format_iris_web_service()}\n"
                        f"   Mag: {mag_type} {magnitude:.2f}\n"
                        f"   Depth: {depth:.2f} km\n"
                        f"{'-'*10}\n"
                        f"{sta_id} ({self.sta_lat:.2f}, {self.sta_lon:.2f})\n"
                        ),
                     x=0.105, y=0.95, horizontalalignment="left",
                     verticalalignment="top", transform=self.ax.transAxes,
                     zorder=5
                     )

    def auxiliaries(self):
        """
        Plot any additional auxiliary information
        """
        pass

    def plot(self):
        """
        Convenience function to plot all the requiste parts of the plot
        """
        self.setup()
        self.source()
        self.receiver()
        self.connecting_line()
        self.annotations()
        self.auxiliaries()

        plt.show()
