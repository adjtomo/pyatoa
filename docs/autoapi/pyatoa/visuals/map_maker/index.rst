:py:mod:`pyatoa.visuals.map_maker`
==================================

.. py:module:: pyatoa.visuals.map_maker

.. autoapi-nested-parse::

   Map making functionality for the Pyatoa package. Creates Cartopy basemaps
   featuring events and stations w/ optional moment tensors



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.visuals.map_maker.MapMaker



Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.map_maker.scale_bar



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.visuals.map_maker.DEGREE_CHAR


.. py:data:: DEGREE_CHAR
   :annotation: = °

   

.. py:class:: MapMaker(cat, inv, dpi=100, figsize=None, figure=None, gridspec=None, corners=None, corner_buffer_deg=2.0, **kwargs)

   A class to call on the Basemap package to generate a map with
   source-receiver information

   .. py:method:: define_bounding_box(corners=None, corner_buffer_deg=2)

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


   .. py:method:: initiate_figure(figsize=None, dpi=100, figure=None, gridspec=None, **kwargs)

      Create a very barebones minimalist (black and white) map to plot data on

      .. note::
          kwargs passed to projection
          https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html


   .. py:method:: source(fm_type='focal_mechanism')

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


   .. py:method:: receiver()

      Plot the receiver with a standard look


   .. py:method:: connect()

      Plot a connecting line between source and receiver


   .. py:method:: annotate(location='lower-right', anno_latlon=False)

      Annotate event receiver information into bottom right corner of the map

      TODO figure out where to put definition of 'location'

      :type location: str
      :param location: location of the annotation block, available:
          'upper-right', 'lower-right', 'upper-left', 'lower-left', 'center'
      :type anno_latlon: bool
      :param anno_latlon: annotate the latitude and longitude values of the
          source and receiver next to their markers. Not always very clean
          so defaults to off.


   .. py:method:: plot(show=True, save=None, **kwargs)

      Main function to generate the basemap, plot all the components, and
      show or save the figure

      :type show: bool
      :param show: show the figure in the gui
      :type save: str
      :param save: if not None, save to the given path stored in this var.
      :type corners: dict
      :param corners: dict containing corner points, if None, lat lon values
          to be determiend by station and receiver locations



.. py:function:: scale_bar(ax, length=None, location=(0.85, 0.95), linewidth=3, return_length=False, ref_proj=ccrs.PlateCarree())

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


