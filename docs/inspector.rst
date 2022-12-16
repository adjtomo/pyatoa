Bulk Inversion Assessment
==============================

The ``Inspector`` class uses the `Pandas library <https://pandas.pydata.org/>`__
to aggregate, manipulate and visualize misfit information collected in
ASDFDataSets during a seismic inversion run with
`SeisFlows <https://github.com/adjtomo/seisflows>`__.

.. note::

    The inversion does not explicitely have to be run with SeisFlows, but the
    Inspector requires datasets formatted in a specific fashion, which is
    automatically done in SeisFlows

Each entry corresponds to a single misfit window, such that statistical
analysis can be run on all misfit windows in aggregate. Users can compare
misfit statistics between iterations, events, stations etc. to determine
how misfit evolves throughout the course of an inversion.

An example of the Inspector is shown below, where the ``windows`` atributes
stores metadata information with respect to a given misfit window.

.. code:: python

    >>> insp.windows
             event iteration step network station channel component    misfit  length_s      dlnA  window_weight  max_cc_value  relative_endtime  relative_starttime  cc_shift_in_seconds           absolute_starttime             absolute_endtime
    0  2018p130600       i01  s00      NZ     BFZ     HHE         E  0.365397     62.76 -0.709653       5.469764      0.871537             77.07               14.31                 1.08  2018-02-18T07:43:42.437644Z  2018-02-18T07:44:45.197644Z
    1  2018p130600       i01  s00      NZ     BFZ     HHN         N  1.620000     39.15 -0.828518       3.882748      0.991762             77.07               37.92                 1.89  2018-02-18T07:44:06.047644Z  2018-02-18T07:44:45.197644Z
    2  2018p130600       i01  s00      NZ     BFZ     HHZ         Z  0.004050     21.21 -0.903363       2.101535      0.990823             41.46               20.25                 0.00  2018-02-18T07:43:48.377644Z  2018-02-18T07:44:09.587644Z


.. code:: python

    >>> insp.windows.iloc[0]
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

    from pyatoa.scripts.load_example_data import load_example_inspector

    insp = load_example_inspector()

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

    >>> insp.sources
                                           time  magnitude  depth_km  latitude  longitude
    event_id
    2014p952799  2014-12-19T12:51:22.480000Z       4.76   33.0664  -38.2520   177.9926
    2013p617227  2013-08-17T08:58:40.320000Z       5.21   15.0781  -41.7326   174.0643

The ``srcrcv`` attribute provides relative information for each source-receiver
pair, including epicentral distance and backazimuth

.. code:: python

    >>> insp.srcrcv
              event network station  distance_km  backazimuth
    0   2014p952799      NZ     BFZ   308.576683    29.701984
    1   2014p952799      NZ     BKZ   165.256199    52.610982
    2   2014p952799      NZ    ETVZ   221.435082    64.421412
    3   2014p952799      NZ    FWVZ   239.506726    63.067781
    4   2014p952799      NZ     HAZ    58.051017   161.539674
    ..          ...     ...     ...          ...          ...
    63  2013p617227      NZ     WAZ   233.019584   199.206495
    64  2013p617227      NZ     WEL    77.038723   229.477199
    65  2013p617227      NZ    WHVZ   301.173355   204.909266
    66  2013p617227      NZ     WIZ   538.670140   208.880622
    67  2013p617227      NZ    WSRZ   538.802628   208.756817

    [68 rows x 5 columns]


Misfit Windows
--------------

Misfit window information is stored in the ``windows`` attribute. Each row in
the window dataframe attribute corresponds to a single misfit window and
contains metadata for the source and receiver used to generate it.

.. code:: python

    >>> insp.windows
               event iteration  ...           absolute_starttime             absolute_endtime
    0    2014p952799       i01  ...  2014-12-19T12:51:49.315000Z  2014-12-19T12:53:15.880000Z
    1    2014p952799       i01  ...  2014-12-19T12:51:29.812500Z  2014-12-19T12:53:12.110000Z
    2    2014p952799       i01  ...  2014-12-19T12:51:34.380000Z  2014-12-19T12:53:12.110000Z
    3    2014p952799       i01  ...  2014-12-19T12:51:44.675000Z  2014-12-19T12:52:47.532500Z
    4    2014p952799       i01  ...  2014-12-19T12:51:28.362500Z  2014-12-19T12:52:47.532500Z
    ..           ...       ...  ...                          ...                          ...
    709  2013p617227       i01  ...  2013-08-17T08:59:06.937500Z  2013-08-17T09:00:24.295000Z
    710  2013p617227       i01  ...  2013-08-17T08:59:02.660000Z  2013-08-17T09:00:19.727500Z
    711  2013p617227       i01  ...  2013-08-17T08:58:42.577500Z  2013-08-17T09:00:08.997500Z
    712  2013p617227       i01  ...  2013-08-17T08:58:56.135000Z  2013-08-17T08:59:51.452500Z
    713  2013p617227       i01  ...  2013-08-17T08:58:52.872500Z  2013-08-17T08:59:54.062500Z

    [714 rows x 17 columns]


Users can query a singvale column of each dataframe to gather information in
array format. For example, to get the max cross correlation value of all
windows in your inversion:


.. code:: python

    >>> insp.windows["max_cc_value"]
    0      0.918584
    1      0.904824
    2      0.919713
    3      0.967400
    4      0.969720
             ...
    709    0.871737
    710    0.737115
    711    0.923461
    712    0.928504
    713    0.925466
    Name: max_cc_value, Length: 714, dtype: float64

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

   >>> insp.misfit(level="station")
                                        unscaled_misfit  nwin     misfit
    iteration step event       station
    i01       s00  2013p617227 BFZ            51.934378     3  17.311459
                               BKZ            82.107746     3  27.369249
                               FWVZ           87.990781     4  21.997695
                               HAZ            21.898865     3   7.299622
                               HIZ            71.910561     3  23.970187
    ...                                             ...   ...        ...
              s03  2014p952799 WAZ            42.896701     3  14.298900
                               WEL            49.171292     2  24.585646
                               WHVZ           38.019649     3  12.673216
                               WIZ             8.919856     3   2.973285
                               WSRZ            8.919856     3   2.973285

    [221 rows x 3 columns]

.. note::

    Available misfit `levels` are: 'step', 'event', and 'station'


Window Number
-------------
The ``nwin`` attribute provides information about the number of misfit windows,
and overall window length (in units of time) gathered during a single iteration.
This is useful for understanding how much of your waveforms have been windowed
during an inversion.

.. code:: python

    insp.nwin(level="step")
                    nwin    length_s
    iteration step
    i01       s00    206  14760.4200
              s01    204  14682.1200
              s02    207  14612.4475
              s03     97   7265.7325

.. note::

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

    >>> insp.stats(choice="max", key="length_s")
    iteration  step  event
    i01        s00   2013p617227    121.5100
                     2014p952799    125.8600
               s01   2013p617227    121.6550
                     2014p952799    125.9325
               s02   2013p617227    121.8000
                     2014p952799    127.6725
               s03   2013p617227     86.4200
                     2014p952799    125.9325
    Name: length_s, dtype: float64

The above code snippet will return the maximum window length for each event in
your inversion and for each iteration.

Minmax
------

This simple argument simple prints out minimum and maximum values for each
column entry for the entire inversion.

.. code:: python

    >>> insp.minmax(pprint=True)
    nwin:                      97.0000
    len:                       7265.7325
    misfit_min:                0.0096
    misfit_max:                35.9764
    misfit_mean:               12.2199
    misfit_median:             12.5125
    misfit_std:                8.4696
    length_s_min:              18.4150
    length_s_max:              125.9325
    length_s_mean:             74.9045
    length_s_median:           81.3450
    length_s_std:              23.1415
    dlnA_min:                  -0.1859
    dlnA_max:                  1.0395
    dlnA_mean:                 0.3850
    dlnA_median:               0.3664
    dlnA_std:                  0.2340
    max_cc_value_min:          0.7075
    max_cc_value_max:          0.9968
    max_cc_value_mean:         0.9239
    max_cc_value_median:       0.9286
    max_cc_value_std:          0.0560
    cc_shift_in_seconds_min:   -0.1450
    cc_shift_in_seconds_max:   9.4975
    cc_shift_in_seconds_mean:  4.7783
    cc_shift_in_seconds_median: 4.9300
    cc_shift_in_seconds_std:   1.9101

Compare Iterations
------------------

The ``compare`` function allows the User to compare different iterations in
their inversion. This is useful when comparing misfit of your initial and final
models.

.. note::

    By default, ``compare`` will compare the first and last iteration in an
    inversion (assumed to be initial and final models)

.. code:: python

    >>> insp.compare(iteration_a="i01", step_count_a="s00",
    >>>              iteration_b="i01", step_count_b="s01")
                 nwin_i01s00  misfit_i01s00  nwin_i01s01  misfit_i01s01  diff_misfit  diff_nwin
    event
    2014p952799           93       6.114807           92       5.697945    -0.416862         -1
    2013p617227          113       7.134975          112       7.173721     0.038745         -1


Comparing Windows
-----------------

The compare windows function finds differences between matching misfit windows
for two iterations in your inversion.

.. note::
    Using this function requires that the two compared iterations have the same
    window choices, i.e., the windows from evaluation A must have been used
    during evaluation B.


.. code:: python

    >>> insp.compare_windows(iteration_a="i01", step_count_a="s00",
    >>>                      iteration_b="i01", step_count_b="s00")
               event network  ... diff_max_cc_value diff_cc_shift_in_seconds
    0    2014p952799      NZ  ...               0.0                      0.0
    1    2014p952799      NZ  ...               0.0                      0.0
    2    2014p952799      NZ  ...               0.0                      0.0
    3    2014p952799      NZ  ...               0.0                      0.0
    4    2014p952799      NZ  ...               0.0                      0.0
    ..           ...     ...  ...               ...                      ...
    201  2013p617227      NZ  ...               0.0                      0.0
    202  2013p617227      NZ  ...               0.0                      0.0
    203  2013p617227      NZ  ...               0.0                      0.0
    204  2013p617227      NZ  ...               0.0                      0.0
    205  2013p617227      NZ  ...               0.0                      0.0

    [206 rows x 17 columns]

.. note::

    The above code block compares an evaluation with itself because the test
    data does not have the same window selection for compared evaluations.


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
Check out the `Gallery <gallery.html>`__ for figure examples and example code
snippets used to generate them.
