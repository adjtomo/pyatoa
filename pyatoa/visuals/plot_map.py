"""misfit visualization tool to be called through adjointBuilder
produces a basemap with beachball and all available stations as well as the
relevant station highlighted. important information annotated (e.g.
misift information, distance, BAz etc.)
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach

from pyatoa.utils.operations.source_receiver import gcd_and_baz
from pyatoa.utils.operations.calculations import myround


mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def build_colormap(array):
    """build a custom range colormap
    """
    vmax = myround(np.nanmax(array), base=1, choice='up')
    vmin = myround(np.nanmin(array), base=1, choice='down')
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.plasma
    colormap = cm.ScalarMappable(norm=norm, cmap=cmap)

    return colormap


def plot_hikurangi_trench(m, path):
    """trace the hikurangi trench on a basemap object 'm'
    """
    trenchcoords = np.load(path)
    lats = trenchcoords['LAT']
    lons = trenchcoords['LON']
    x, y = m(lons, lats)

    # interpolate points to make a smoother curve
    xprime = np.flip(x, axis=0)
    yprime = np.flip(y, axis=0)
    xprimenew = np.linspace(x.min(), x.max(), 100)
    yprimenew = np.interp(xprimenew, xprime, yprime)

    m.plot(xprimenew, yprimenew, '--', linewidth=1.25, color='k', zorder=2)


def plot_active_faults(m, path):
    """plot onshore and offshore fault coordinate files
    """
    active_faults = np.load(path)
    lats = active_faults['LAT']
    lons = active_faults['LON']
    faults = active_faults['FAULT']

    for i in range(faults.min(), faults.max()+1, 1):
        indices = np.where(faults==i)
        x, y = m(lons[indices], lats[indices])
        m.plot(x, y, '--', linewidth=0.5, color='k', zorder=2, alpha=0.25)


def event_beachball(m, moment_tensor):
    """plot event beachball on basemap 'm' object for a given geonet event_id
    """
    eventx, eventy = m(moment_tensor['Longitude'], moment_tensor['Latitude'])
    focal_mechanism = [moment_tensor['strike2'], moment_tensor['dip2'],
                       moment_tensor['rake2']]
    b = beach(focal_mechanism, xy=(eventx, eventy), width=2.5E4, linewidth=1,
              facecolor='r')
    b.set_zorder(1000)
    ax = plt.gca()
    ax.add_collection(b)


def initiate_basemap(map_corners):
    """
    set up local map
    """
    continent_color = 'w'
    lake_color = 'w'

    # initiate map
    m = Basemap(projection='stere', resolution='h', rsphere=6371200,
                lat_0=np.mean(map_corners[:2]), lon_0=np.mean(map_corners[2:]),
                llcrnrlat=map_corners[0], llcrnrlon=map_corners[2],
                urcrnrlat=map_corners[1], urcrnrlon=map_corners[3],
                )
    m.drawcoastlines(linewidth=1.5)
    m.fillcontinents(color=continent_color, lake_color=lake_color)
    m.drawparallels(np.arange(int(map_corners[0]),int(map_corners[1]), 1),
                    labels=[1, 0, 0, 0], linewidth=0.0)
    m.drawmeridians(np.arange(int(map_corners[2]),int(map_corners[3])+1, 1),
                    labels=[0, 0, 0, 1], linewidth=0.0)

    return m


def populate_basemap(m, lats, lons, names=None, nets=None):
    """fill map with latitude/longitude pairs, i.e. stations, events
    """
    x, y = m(lons, lats)

    # manual set colors based on network
    colordict = {"NZ": 'k', "XX": 'k', "X?": 'gray'}
    edgedict = {"NZ": '-', "XX": '--', "X?": '--'}
    edgestyle = '-'
    if nets:
        colors, edgestyle = [], []
        for net in nets:
            colors.append(colordict[net])
            edgestyle.append(edgedict[net])
    m.scatter(x,y, marker='v', color='w', edgecolor='k', linestyle=edgestyle,
              s=60, zorder=5)
    if names:
        for n_, x_, y_ in zip(names, x, y):
            plt.annotate(n_, xy=(x_, y_), xytext=(x_, y_), zorder=6, fontsize=7,
                         bbox=dict(facecolor='w', edgecolor='k',
                                   boxstyle='round')
                         )


def annotate_information(m, event_id, event, inv=None):
    """annotate event information into hard coded map area
    """
    event.origins[0].time.precision = 0
    for magni in event.magnitudes:
        if magni.magnitude_type == "M":
            magnitude = magni.mag
            magnitude_type = magni.magnitude_type
    if inv:
        annotemplate = (
            "{id} / {sta}\n{date}\n{type}={mag:.2f}"
            "\nDepth(km)={depth:.2f}\nDist(km)={dist:.2f}\nBAz(deg)={baz:.2f}")
        gcdist, baz = gcd_and_baz(event, inv)
    else:
        annotemplate = ("{id} / {sta}\n{date}\n"
                        "{type}={mag:.2f}\nDepth(km)={depth:.2f}")

    plt.annotate(
        s=annotemplate.format(id=event_id, sta=inv[0][0].code,
                              date=event.origins[0].time,
                              stalat=inv[0][0].latitude,
                              stalon=inv[0][0].latitude,
                              evlat=event.origins[0].latitude,
                              evlon=event.origins[0].longitude,
                              depth=event.origins[0].depth*1E-3,
                              type=magnitude_type, mag=magnitude,
                              dist=gcdist, baz=baz
                              ),
        xy=((m.xmax-m.xmin)*0.725, (m.ymax-m.ymin)*0.0375),
        multialignment='right', fontsize=10)


def generate_map(config, event, inv,
                 map_corners=[-42.5007,-36.9488,172.9998,179.5077],
                 show_faults=False):
    """
    TODO: change map corners to reflect the new mesh created in August

    initiate and populate a basemap object for New Zealands north island.
    Functionality to manually ignore stations based on user quality control
    Takes station coordinates and coloring from npz files
    Choice to annotate two stations which correspond to waveforms
    Calls beachball and trench tracer to populate basemap
    :type corners: list of floats
    :param corners: values for map corners to set bounds
     e.g. [lat_bot,lat_top,lon_left,lon_right]
    """
    m = initiate_basemap(map_corners=map_corners)
    annotate_information(m, config.event_id, event, inv)
    if show_faults:
        plot_hikurangi_trench(m,
            os.path.join(config.paths['faults'], "hikurangiTrenchCoords.npz")
                              )
        plot_active_faults(m,
            os.path.join(config.paths['faults'], "onshoreFaultCoords.npz")
                           )
        plot_active_faults(m,
            os.path.join(config.paths['faults'], "offshoreFaultCoords.npz")
                           )

    stationfile = (pathnames()['stationxml'] + 
                                    'COORDINATE/masterlist/masterlist.npz')    
    stationlist = np.load(stationfile)

    populate_basemap(m,lats=stationlist['LAT'],
                       lons=stationlist['LON'],
                       nets=stationlist['NET'])    
    scalelon,scalelat = 178.75,-37.2
    m.drawmapscale(scalelon,scalelat,scalelon,scalelat,100,
                                                yoffset=0.01*(m.ymax-m.ymin))

    return m

    return m


if __name__ == "__main__":

