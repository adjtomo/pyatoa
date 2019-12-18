"""
Cartopy based mapping tools for regional map projections
"""
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


def init_map(map_corners):
    """
    Initiate the map in orthographic projection based on map corners.
    Central latitude and longitude values will be inferred from map corners.

    :type map_corners: list of floats
    :param map_corners: [lon min, lon max, lat min, lat max]
    :return:
    """
    central_lon = (map_corners[1] - map_corners[0]) / 2
    central_lat = (map_corners[3] - map_corners[2]) / 2

    ax = plt.axes(projection=ccrs.Orthographic(central_lon, central_lat))
    ax.set_extent(map_corners)
    ax.gridlines()
    ax.coastlines("50m")

    return ax
