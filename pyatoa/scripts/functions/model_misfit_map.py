import os
import glob
import json
import pyasdf
import numpy as np
import matplotlib.pyplot as plt

from pyatoa.utils.visuals.map_tools import initiate_basemap, \
    interpolate_and_contour, event_beachball
from pyatoa.utils.asdf.extractions import sum_misfits


def read_datasets(path_to_ds="./", model="m00", save="./misfitinfo.json"):
    """
    Read the datasets that were created by the Pyatoa workflow and skim the
    necessary information from them

    :type path_to_ds: str
    :param path_to_ds: path to the Pyasdf datsets
    :type model: str
    :param model: model number to read in, e.g. "m00"
    :type save: str
    :param save: file id to save the .npz file for easier retrieval
    :return:
    """
    misdict = {}
    if not os.path.exists(save):
        for dataset in glob.glob(os.path.join(path_to_ds, '*.h5')):
            with pyasdf.ASDFDataSet(dataset) as ds:
                event_id = ds.events[0].resource_id.id.split('/')[-1]
                ev_lat = ds.events[0].preferred_origin().latitude
                ev_lon = ds.events[0].preferred_origin().longitude

                misdict[event_id] = {"lat": ev_lat, "lon": ev_lon}
                if hasattr(ds, 'waveforms'):
                    for sta in ds.waveforms.list():
                        if 'AdjointSources' in ds.auxiliary_data:
                            misfit = sum_misfits(ds, model, station=sta)
                            if not misfit:
                                continue
                            lat = ds.waveforms[sta].StationXML[0][0].latitude
                            lon = ds.waveforms[sta].StationXML[0][0].longitude
                            misdict[event_id][sta] = {"misfit": misfit,
                                                      "lat": lat, "lon": lon}

        if save:
            with open(save, 'w') as f:
                json.dumps(misdict, f, indent=4)

    else:
        with open(save, 'r') as f:
            misdict = json.load(f)

    return misdict


if __name__ == "__main__":
    model = "m00"

    markersize = 75
    linewidth = 1.75
    alpha = 1.0
    map_corners = {'lat_min': -42.5007, 'lat_max': -36.9488,
                   'lon_min': 172.9998, 'lon_max': 178.6500}

    misdict = read_datasets(path_to_ds=os.getcwd(), model=model)

    # Initiate matplotlib instances
    f = plt.figure(figsize=(8, 10), dpi=100)
    m = initiate_basemap(map_corners=map_corners, scalebar=True)

    # Plot all the event locations
    sta_misfit = {}
    for event_id in misdict.keys():
        ev_x, ev_y = m(misdict[event_id]['lon'], misdict[event_id]['lat'])
        m.scatter(ev_x, ev_y, marker='o', color='r', alpha=alpha,
                  edgecolor='k', linestyle='-', s=markersize,
                  linewidth=linewidth, zorder=100
                  )

        for sta in misdict[event_id].keys():
            if sta in ['lat', 'lon']:
                continue
            # put misfits in dict object
            if sta in sta_misfit:
                sta_misfit[sta].append(misdict[event_id][sta]['misfit'])
            else:
                sta_misfit[sta] = [misdict[event_id][sta]['misfit']]

    x_plot, y_plot, misfit_plot = [], [], []
    for sta in sta_misfit.keys():
        sta_x, sta_y = m(misdict[event_id][sta]['lon'],
                         misdict[event_id][sta]['lat'])
        x_plot.append(sta_x)
        y_plot.append(sta_y)

        # get some variables from the misfit
        mean_misfit = sta_misfit[sta].mean()
        total_misfit = sta_misfit[sta].sum()
        max_misfit = sta_misfit[sta].max()

        misfit_plot.append(total_misfit)

    # plot the misfit map
    interpolate_and_contour(m=m, x=x_plot, y=y_plot, z=misfit_plot,
                            len_xi=100, len_yi=150, colormap='viridis',
                            interpolation_method='cubic',
                            fill=True, cbar_label=''
                            )

    plt.title()

    plt.show()
