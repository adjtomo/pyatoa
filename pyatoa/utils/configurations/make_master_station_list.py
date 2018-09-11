"""
Pyatoa uses a master station list to determine a few parameters such as station plotting, data availablility etc.
This script will help generate a station list in the correct format, such that pyatoa can easily retrieve information
"""
from obspy import read_inventory


def make_master_station_list():
    """
    master function to create station list
    :return:
    """
    template = ("{NETWORK},{STATIONCODE},{LOCATIONCODE},{CHANNELCODE},"
            "{LATITUDE},{LONGITUDE},{STARTTIME},{ENDTIME}"
                )
    adfadf


