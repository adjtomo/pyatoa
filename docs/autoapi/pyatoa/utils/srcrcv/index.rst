:py:mod:`pyatoa.utils.srcrcv`
=============================

.. py:module:: pyatoa.utils.srcrcv

.. autoapi-nested-parse::

   Functions used in determining information related to sources and receivers,
   or their corresponding representations in ObsPy



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.srcrcv.lonlat_utm
   pyatoa.utils.srcrcv.utm_zone_from_lat_lon
   pyatoa.utils.srcrcv.gcd_and_baz
   pyatoa.utils.srcrcv.merge_inventories
   pyatoa.utils.srcrcv.seismogram_length
   pyatoa.utils.srcrcv.sort_by_backazimuth
   pyatoa.utils.srcrcv.event_by_distance



.. py:function:: lonlat_utm(lon_or_x, lat_or_y, utm_zone=None, inverse=False)

   Convert latitude and longitude coordinates to UTM projection using PyProj

   :type lon_or_x: float or int
   :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
   :type lat_or_y: float or int
   :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
   :type utm_zone: int
   :param utm_zone: UTM zone for conversion from WGS84
   :type inverse: bool
   :param inverse: if inverse == False, latlon => UTM, vice versa.
   :rtype: tuple (float, float)
   :return: (x in UTM or longitude in WGS84, y in UTM or latitude in WGS84)


.. py:function:: utm_zone_from_lat_lon(lat, lon)

   Calculate the UTM zone longitude value using quick maffs.
   Get the sign of the UTM zone based on the latitude value.

   :type lat: float
   :param lat: latitude coordinate in degrees
   :type lon: float
   :param lon: longitude coordinate in degrees
   :rtype: int
   :return: UTM zone number


.. py:function:: gcd_and_baz(event, sta)

   Calculate great circle distance and backazimuth values for a given
   station and event configuration

   :type event: obspy.core.event.Event
   :param event: event object
   :type sta: obspy.core.inventory.station.Station
   :param sta: station object
   :rtype: tuple (float, float)
   :return: (great circle distance in km, backazimuth in degrees)


.. py:function:: merge_inventories(inv_a, inv_b)

   Adding inventories together duplicates network and station codes, which is
   kind of annoying for looping. This function will add two inventories
   together while minimizing the amount of redundant networks, stations,
   channels inside the merged inventory.

   :type inv_a: obspy.core.inventory.Inventory
   :param inv_a: inventory to merge into, will be returned
   :type inv_b: obspy.core.inventory.Inventory
   :param inv_b: inventory to merge into inv_a
   :rtype: obspy.core.inventory.Inventory
   :return: merged inventories


.. py:function:: seismogram_length(distance_km, slow_wavespeed_km_s=2, binsize=50, minimum_length=100)

   Dynamically determine the length of the seismogram based on source-receiver
   distance. Bin into lengths to keep some uniformity in window lengths

   :type distance_km: float
   :param distance_km: source-receiver distance in km
   :type slow_wavespeed_km_s: int
   :param slow_wavespeed_km_s: slowest wavespeed in model, in km/s
   :type binsize: int
   :param binsize: bin size for rounding the length to the nearest value
   :type minimum_length: int
   :param minimum_length: the shortest a seismogram can be
   :rtype: int
   :return: expected seismogram length


.. py:function:: sort_by_backazimuth(ds, clockwise=True)

   Its illustrative to show misfit information for an event, sorted by
   backazimuth. Stations with misfit information are generally sorted
   alphabetically, so this function just calcualtes backazimuth and returns a
   sorted list of station names. Can go clockwise or counter clockwise,
   starting from 0 degrees.

   :type ds: pyasdf.ASDFDataSet()
   :param ds: dataset containing event and station information
   :type clockwise: bool
   :param clockwise: False = counter clockwise
   :rtype: list
   :return: list of stations in order from 0deg to 360deg in direction


.. py:function:: event_by_distance(cat, filter_type=False, filter_bounds=None, random=False)

   Sort through an obspy catalog by interevent distance. If we have a lot of
   events in a catalog, it's best to spatially vary them such that we don't
   redundantly oversample a spatial region.
   Returns an index list for events that are most distant from one another,
   without repeating any used events.

   Catalog filter parameters can be found here:
   https://docs.obspy.org/packages/autogen/obspy.core.event.Catalog.filter.html

   :type cat: obspy.event.Catalog
   :param cat: catalog to sort through
   :type filter_type: str
   :param filter_type: filter to be passed to the Catalog filter
   :type filter_bounds: list of floats
   :param filter_bounds: (min filter bound, max filter bound)
   :type random: bool
   :param random: randomly determined starting point
   :rtype: obspy.event.Catalog
   :return: filtered catalog object


