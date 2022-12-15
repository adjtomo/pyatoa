SeisFlows Inversion Assessment
==============================

The ``Inspector`` class uses the [Pandas library](https://pandas.pydata.org/>)
to aggregate, manipulate and visualize misfit information collected in
ASDFDataSets during a seismic inversion run with
`SeisFlows <https://github.com/adjtomo/seisflows>`__.

Each entry corresponds to a single misfit window, such that statistical
analysis can be run on all misfit windows in aggregate. Users can compare
misfit statistics between iterations, events, stations etc. to determine
how misfit evolves throughout the course of an inversion.

An example of the Inspector is shown below, where the ``windows`` atributes
stores metadata information with respect to a given misfit window.

.. code:: python

    insp.windows


.. parsed-literal::

             event iteration step network station channel component    misfit  length_s      dlnA  window_weight  max_cc_value  relative_endtime  relative_starttime  cc_shift_in_seconds           absolute_starttime             absolute_endtime
    0  2018p130600       i01  s00      NZ     BFZ     HHE         E  0.365397     62.76 -0.709653       5.469764      0.871537             77.07               14.31                 1.08  2018-02-18T07:43:42.437644Z  2018-02-18T07:44:45.197644Z
    1  2018p130600       i01  s00      NZ     BFZ     HHN         N  1.620000     39.15 -0.828518       3.882748      0.991762             77.07               37.92                 1.89  2018-02-18T07:44:06.047644Z  2018-02-18T07:44:45.197644Z
    2  2018p130600       i01  s00      NZ     BFZ     HHZ         Z  0.004050     21.21 -0.903363       2.101535      0.990823             41.46               20.25                 0.00  2018-02-18T07:43:48.377644Z  2018-02-18T07:44:09.587644Z


.. code:: python

    insp.windows.iloc[0]


.. parsed-literal::

    event                                  2018p130600
    iteration                                      i01
    step                                           s00
    network                                         NZ
    station                                        BFZ
    channel                                        HHE
    component                                        E
    misfit                                    0.365397
    length_s                                     62.76
    dlnA                                     -0.709653
    window_weight                             5.469764
    max_cc_value                              0.871537
    relative_endtime                             77.07
    relative_starttime                           14.31
    cc_shift_in_seconds                           1.08
    absolute_starttime     2018-02-18T07:43:42.437644Z
    absolute_endtime       2018-02-18T07:44:45.197644Z
    Name: 0, dtype: object

Dataset Discovery
~~~~~~~~~~~~~~~~~

During a SeisFlows Inversion workflow run with preprocessing module
``Pyaflowa``, ASDFDataSets containing waveform data, misfit windows and
adjoint sources are generated.

The Inspector requires a path to this director to gather misfit information.
Use the `discover` function to aggregate data.

.. code:: python

    from pyatoa import Inspector

    insp = Inspector(tag="example")
    insp.discover(path="{path_to_asdf_datasets}")

.. note::

    The `discover` function will wildcard search the given `path` for any files
    ending in `.h5`, the default file extension for ASDFDataSets

Example Data
------------

For those that would like example data to test the Inspector without running
a SeisFlows inversion, you can generate an example dataset with the following
command:

.. code:: python

    from pyatoa import Manager, Inspector
    from pyasdf import ASDFDataSet

    ds = ASDFDataSet("example")

    mgmt = Manager(ds=ds)
    mgmt.load()  # loads example data
    mgmt.flow()
    mgmt.write()  # writes data into dataset

    insp = Inspector("example")
    insp.discover(path="./")

A real-inversion Inspector can also be generated from test data stored within
the Pyatoa repository

.. code:: python

    from pyatoa import Inspector

    insp.read(path="{path_to_pyatoa}/pyatoa/tests/test_data")

Data Attribute Access
~~~~~~~~~~~~~~~~~~~~~

Pandas Dataframes are like spreadsheets, storing data in row-column format.
During data discovery, the Inspector retrieves source and receiver
metadata, misfit windows information (e.g., starttimes, time shifts, etc.), and
adjoint source information (e.g., total misfit).

Source and receiver metadata
----------------------------

A list of event ids and station names can be accessed through the
``events`` and ``stations`` attributes.

.. code:: python

    insp.events  # returns list of event ids
    insp.stations  # returns list of station ids


Source and receiver metadata like hypocentral location are accesible through
the ``sources`` and ``receivers`` attributes.

.. code:: python

    insp.sources

.. parsed-literal::

                                        time  magnitude   depth_km   latitude   longitude
    event_id
    2018p130600  2018-02-18T07:43:48.127644Z   5.156706  20.594599 -39.948975  176.299515


The ``srcrcv`` attribute provides relative information for each source-receiver
pair, including epicentral distance and backazimuth

.. code:: python

    insp.srcrcv


Misfit Windows
~~~~~~~~~~~~~~

Misfit window information is stored in the ``windows`` attribute. Each row in
the window dataframe attribute corresponds to a single misfit window and
contains metadata for the source and receiver used to generate it.

.. code:: python

    insp.windows

Users can query a single column of each dataframe to gather information in
array format. For example, to get the max cross correlation value of all
windows in your inversion:


.. code:: python

    insp.windows["max_cc_value"]

Data Access Functions
~~~~~~~~~~~~~~~~~~~~~

Using Pandas syntax, the User should be able to get at any permutation
of data that they want to analyze, however the Inspector has some built-in
data access functions.

Misfit Information
-------------------

The ``misfit`` function calculates total misfit for various levels (e.g., per
iteration, per station, per event)

.. code:: python

    insp.misfit(level="station")

Available misfit `levels` are: 'step', 'event', and 'station'


Window Number
-------------
The ``nwin`` attribute provides information about the number of misfit windows,
and overall window length (in units of time) gathered during a single iteration.
This is useful for understanding how much of your waveforms have been windowed
during an inversion.

.. code:: python

    insp.nwin(level="step")

Available `levels` are: 'step', 'event', and 'station'

Window Statistics
-----------------

The ``stats`` function aggregates all the columns into a per-evaluation,
per-event calculation. That is, for every event in every iteration, column
values like overall misfit, or total window number, will be averaged.

The default ‘stat’ function takes the mean, but other NumPy statistical
functions like mean, std (standard deviation) or var (variance) can also
be applied.

.. code:: python

    insp.stats(choice="max", key="length_s")

The above code snippet will return the maximum window length for each event in
your inversion and for each iteration.

Minmax
~~~~~~

This simple argument simple prints out minimum and maximum values for each
column entry for the entire inversion.

.. code:: python

    insp.minmax()

Compare Iterations
~~~~~~~~~~~~~~~~~~~

The ``compare`` function allows the User to compare different iterations in
their inversion. This is useful when comparing misfit of your initial and final
models.

.. note::

    By default, ``compare`` will compare the first and last iteration in an
    inversion (assumed to be initial and final models)

.. code:: python

    insp.compare(iteration_a="i01", step_count_a="s00",
                 iteration_b="i01", step_count_b="s01"
                 )


Comparing Windows
~~~~~~~~~~~~~~~~~

The compare windows function finds differences between matching misfit windows
for two iterations in your inversion.

.. note::
    Using this function requires that the two compared iterations have the same
    window choices, i.e., the windows from evaluation A must have been used
    during evaluation B.


.. code:: python

    insp.compare_windows(iteration_a="i01", step_count_a="s00",
                         iteration_b="i01", step_count_b="s01"
                         )


Manipulating Inspector Objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following section will show you how to manipulate the Inspector object
itself, e.g., read/write to disk, add data from new datasets, merge two
inspectors

Read/Write From Disk
--------------------

The Inspector can be written to disk as a collection of comma separated value
(.csv) files, which can be opened with spreadsheet software (e.g., Excel).

.. code:: python

    insp.write(path="./", fmt="csv", tag="example")

Inspectors can also re-read these output files

.. code:: python

    insp.read(path="./", fmt="csv", tag="example")


Add New Data to Inspector
--------------------------

New datasets can be added to an existing Inspector class with the ``append``
function.

.. code:: python

    insp.append(dsfid="{path_to_asdfdataset}")

Merging Two Inspectors
----------------------

During very large inversions, it may be useful to split the inversion into
various stages or legs, resulting in multiple sets of related ASDFDataSets.

The ``extend`` function allows you to aggregate measurements from all
inversion stages.

.. code:: python

    from pyatoa import Inspector

    insp_stage1 = Inspector(tag="stage1")
    insp_stage1.discover("{path_to_stage1_datasets}")

    insp_stage2 = Inspector(tag="stage2")
    insp_stage2.discover("{path_to_stage2_datasets}")

    insp_stage1.extend(insp_stage2.windows)


Plotting Routines
~~~~~~~~~~~~~~~~~

The Inspector comes with a suite of plotting routines to visualize the dataset.
Check out the `gallery <gallery.html>`__ for examples and code snippets to
generate them.

Source-receiver Metadata
-------------------------

A very simple source-receiver scatter plot can be created with the ``map``
function

.. code:: python

    insp.map()

The ``event_depths`` functions plots a 2D cross section of all events at depth

.. code:: python

    insp.event_depths(xaxis="longitude")

The ``raypaths`` function shows connecting lines for any source-receiver pair
that has atleast one measurement

.. code:: python

    insp.raypaths(iteration="i01", step_count="s00")

The ``raypath_density`` function provides a more detailed raypath plot, which is
colored by the density of overlapping raypaths

.. code:: python

    insp.raypath_density(iteration="i01", step_count="s00")

The ``event_hist`` function creates a simple event histogram based on event
information such as magnitude.

.. code:: python

    insp.sources.keys()  #  <- Use to check available choices
    insp.event_hist(choice="magnitude")



Misfit Window Timing
---------------------

The following plotting functions are concerned with visualizing the time
dependent part of the measurements

The ``travel_times`` function plots a proxy for phase arrivals, similar to
a seismic record section.

.. code:: python

    insp.travel_times(markersize=2, t_offset=-20, constants=[2, 4, 6, 8, 10])

The ``plot_windows`` function plots time windows (as bars) against source
receiver distance, illustrating seismic phases included in the inversion.


.. code:: python

    insp.plot_windows(iteration="i01", step_count="s00")


Inversion Statistics
--------------------

The following plotting functions help the user understand how an inversion is
progressing by comparing iterations against one another. These are common
inversion statistics plots shown in many tomography publications.


The ``convergence`` function plots total misfit per iteration over the course
of an inversion. An additional Y axis is used to plot the number of windows for
each iteration (or the overall length of the time windows)

.. code:: python

    insp.convergence(windows="nwin")


The ``hist`` function generates histograms for a given measurement column,
such as overall cross correlation or amplitude anomaly.

.. code:: python

    insp.windows.keys()  # <- To see available choices
    insp.hist(choice="cc_shift_in_seconds")

The ``hist`` function can also be used to generate two sets of histograms that
compare one iteration to another:

.. code:: python

    insp.hist(iteration="i01", step_count="s00", iteration_comp="i01",
              step_count_comp="s01", choice="dlnA")


Measurement Statistics
-----------------------

These plotting functions allow the user to plot measurements for a given
evaluation in order to better understand the statistical distribution of
measurements, or comparisons against one another.

The ``scatter`` function compares any two attributes in the `windows` dataframe

.. code:: python

    insp.scatter(x="relative_starttime", y="max_cc_value")

The ``measurement_hist`` function generates histograms of source or receiver
metadata. Useful for identifying events or stations which may be outliers in
terms of overall measurements.

.. code:: python

    insp.measurement_hist(iteration="i01", step_count="s00", choice="station")

The ``station_event_misfit_map`` creates a map for a single station. All other
points correspond to events which the station has recorded. Colors of these
markers correspond to given measurement criteria.

.. code:: python

    insp.station_event_misfit_map(station="BFZ", iteration="i01",
                                  step_count="s00", choice="misfit")

The ``station_event_misfit_map`` creates a map for a single event. All other
points correspond to stations which have recorded the event. Colors of these
markers correspond to given measurement criteria.

.. code:: python

    insp.event_station_misfit_map(event="2014p952799", iteration="i01",
                                  step_count="s00", choice="nwin", cmap="jet_r")


The ``event_misfit_map`` plots all events on a map and their corresponding
scaled misfit value for a given evaluation (defaults to last evaluation in the
Inspector).

.. code:: python

    insp.event_misfit_map(choice="misfit")
