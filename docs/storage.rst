Saving Data with ASDF
=====================

Pyatoa stores data and processing results to `PyASDF
ASDFDataSets <https://seismicdata.github.io/pyasdf/asdf_data_set.html>`__,
which are seismological data structures built upon the HDF5 file format.

Data are stored as ObsPy or NumPy objects within the ASDFDataSet and can be
easily retrieved for later processing.

Collections of ASDFDataSets can be read by the Inspector class to perform
`bulk misfit assessment <inspector.html>`__.

To load an example dataset to play around with:

.. code:: python

    from pyatoa.scripts.load_example_data import load_example_asdfdataset

    ds = load_example_asdfdataset()

Writing Data to a Dataset
-------------------------

Pyatoa stores data, metadata and processing results in ASDFDataSets. This can
either be done manually, or automatically during a processing workflow.

Writing Data Manually
~~~~~~~~~~~~~~~~~~~~~

The :meth:`write <pyatoa.core.manager.Manager.write>` function of the Manager
writes data to ASDFDataSets, including: observed waveforms, synthetic waveforms,
station metadata, event metadata, misfit windows and adjoint sources.


.. code:: python

    from pyasdf import ASDFDataSet
    from pyatoa import Manager

    ds = ASDFDataSet("example.h5")
    mgmt = Manager(ds=ds)
    # ... some processing steps

    mgmt.write(ds=ds)

The :class:`Config <pyatoa.core.config.Config>` class can also write itself to
a dataset, this must be done separate from the
:class:`Manager <pyatoa.core.manager.Manager>` write.

.. code:: python

    from pyatoa import Config

    cfg = Config()
    cfg.write(write_to=ds)

Writing Data Automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~

Prior to data gathering and processing with the
:class:`Manager <pyatoa.core.manager.Manager>`, you can set the
:class:`Config <pyatoa.core.config.Config>` parameter ``save_to_ds`` to tell a
Manager to automatically save any data and processing results to the dataset
(this is set `True` by default).

The Manager must be supplied a valid ADSFDataSet. See the `data discovery
<discovery.html>`__ page for automated data discovery routines.

.. code::

    from pyatoa import Config, Manager
    from pyasdf import ASDFDataSet

    ds = ASDFDataSet("example.h5")
    cfg = Config(save_to_ds=True, paths={...})  # paths set to local data
    mgmt = Manager(ds=ds, config=cfg)
    mgmt.gather()  # <- automatically stores gathered data to dataset
    mgmt.standardize().preprocess()
    mgmt.window()  # <- automatically stores misfit windows to dataset
    mgmt.measure()  # <- automatically stores adjoint sources to dataset

Accessing Data from Datasets
----------------------------

This section details how to access waveforms, misfit results and metadata stored
inside an ASDFDataSet.

See the `PyASDF documentation
<https://seismicdata.github.io/pyasdf/tutorial.html#reading-an-existing-asdf-data-set>`__
for more information.

Event and Station Metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~

To access ``event`` metadata, stored as an ObsPy Event object

.. note::

    By design, Pyatoa only stores one event per ASDFDataSet, to avoid file
    sizes getting too large;

.. code:: python

    >>> ds.events[0]
    Event:	2018-02-18T07:43:48.130000Z | -39.949, +176.299 | 4.86 mw

                          resource_id: ResourceIdentifier(id="smi:local/cmtsolution/2018p130600/event")
                           event_type: 'earthquake'
                  preferred_origin_id: ResourceIdentifier(id="smi:local/cmtsolution/2018p130600/origin#cmt")
               preferred_magnitude_id: ResourceIdentifier(id="smi:local/cmtsolution/2018p130600/magnitude#mw")
         preferred_focal_mechanism_id: ResourceIdentifier(id="smi:local/cmtsolution/2018p130600/focal_mechanism")
                                 ---------
                   event_descriptions: 1 Elements
                             comments: 1 Elements
                     focal_mechanisms: 1 Elements
                              origins: 2 Elements
                           magnitudes: 3 Elements

To access the ``station`` list, which stores data and metadata for all stations
in the dataset:

.. code:: python

    >>> ds.waveforms.list()
    ['NZ.BFZ']


Waveforms are stored alongside metadata coded by the the network and station
code of each receiver.

.. code:: python

    >>> ds.waveforms.NZ_BFZ
    Contents of the data set for station NZ.BFZ:
        - Has a StationXML file
        - 2 Waveform Tag(s):
            observed
            synthetic_i01s00

To access station metadata, stored as an ObsPy Inventory object

.. code:: python

    >>> ds.waveforms.NZ_BFZ.StationXML
    Inventory created at 2020-02-02T22:21:59.000000Z
        Created by: Delta
                None
        Sending institution: GeoNet (WEL(GNS_Test))
        Contains:
            Networks (1):
                NZ
            Stations (1):
                NZ.BFZ (Birch Farm)
            Channels (3):
                NZ.BFZ.10.HHZ, NZ.BFZ.10.HHN, NZ.BFZ.10.HHE

Observed and Synthetic Waveforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Observed waveforms are tagged by Pyatoa with the ``Config.observed_tag``
attribute, which is 'observed' by default. Waveforms are stored as Stream
objects.

.. code:: python

    >>> ds.waveforms.NZ_BFZ.observed
    3 Trace(s) in Stream:
    NZ.BFZ..BXE | 2018-02-18T07:43:28.130000Z - 2018-02-18T07:49:30.557500Z | 13.8 Hz, 5000 samples
    NZ.BFZ..BXN | 2018-02-18T07:43:28.130000Z - 2018-02-18T07:49:30.557500Z | 13.8 Hz, 5000 samples
    NZ.BFZ..BXZ | 2018-02-18T07:43:28.130000Z - 2018-02-18T07:49:30.557500Z | 13.8 Hz, 5000 samples

Synthetic waveforms are tagged by Pyatoa with the ``Config.synthetic_tag``
attribute.

.. code:: python

    ds.waveforms.NZ_BFZ.synthetic

During a SeisFlows inversion, the ``synthetic_tag`` may reflect the iteration
and step count assigned by SeisFlows.

.. note::

    See the `naming standards page <standards.html>`__ for further explanation
    on tagging for inversions.

For iteration 1, step count 0, synthetics will be saved as:

.. code:: python

    >>> ds.waveforms.NZ_BFZ.synthetics_i01s00
    3 Trace(s) in Stream:
    NZ.BFZ..BXE | 2018-02-18T07:43:28.130000Z - 2018-02-18T07:49:30.557500Z | 13.8 Hz, 5000 samples
    NZ.BFZ..BXN | 2018-02-18T07:43:28.130000Z - 2018-02-18T07:49:30.557500Z | 13.8 Hz, 5000 samples
    NZ.BFZ..BXZ | 2018-02-18T07:43:28.130000Z - 2018-02-18T07:49:30.557500Z | 13.8 Hz, 5000 samples

This tagging system allows Pyatoa to save multiple sets of synthetic waveforms
to a single ASDFDataSet.

Misfit Windows
~~~~~~~~~~~~~~

Misfit windows, adjoint sources and configuration parameters are stored in the
``auxiliary_data`` attribute of the ASDFDataSet.

.. code:: python

    >>> ds.auxiliary_data
    Data set contains the following auxiliary data types:
        AdjointSources (1 item(s))
        Configs (2 item(s))
        MisfitWindows (1 item(s))

The ``MisfitWindows`` attribute stores information about misfit windows

.. code:: python

    ds.auxiliary_data.MisfitWindows

During an inversion, misfit windows are tagged by the ``iter_tag`` and
``step_tag`` attributes of :class:`Config <pyatoa.core.config.Config>`

.. code:: python

    >>> ds.auxiliary_data.MisfitWindows
    1 auxiliary data sub group(s) of type 'MisfitWindows' available:
        i01
    >>> ds.auxiliary_data.MisfitWindows.i01
    1 auxiliary data sub group(s) of type 'MisfitWindows/i01' available:
        s00
    >>> ds.auxiliary_data.MisfitWindows.i01.s00
    3 auxiliary data item(s) of type 'MisfitWindows/i01/s00' available:
        NZ_BFZ_E_0
        NZ_BFZ_N_0
        NZ_BFZ_Z_0

Accessing each misfit window provides a dictionary of window parameters, same
as the information that is outputted by Pyflex.

.. code:: python

    >>> ds.auxiliary_data.MisfitWindows.i01.s00.NZ_BFZ_E_0
    Auxiliary Data of Type 'MisfitWindows'
        Path: 'i01/s00/NZ_BFZ_E_0'
        Data shape: '(2,)', dtype: 'int64'
        Parameters:
            absolute_endtime: 2018-02-18T07:44:59.915000Z
            absolute_starttime: 2018-02-18T07:43:57.130000Z
            cc_shift_in_samples: 97
            cc_shift_in_seconds: 7.0325
            center_index: 833
            channel_id: NZ.BFZ..BXE
            dlnA: 0.8178943677509113
            dt: 0.0725
            left_index: 400
            max_cc_value: 0.9260584412126905
            min_period: 8.0
            phase_arrival_P: 15.262235117775926
            phase_arrival_Pn: 15.131536549180034
            phase_arrival_S: 25.700988089152666
            phase_arrival_Sn: 25.674453184025445
            phase_arrival_p: 14.045597727214647
            phase_arrival_s: 23.62091920350575
            phase_arrival_sP: 18.77953271333086
            relative_endtime: 91.785
            relative_starttime: 28.999999999999996
            right_index: 1266
            time_of_first_sample: 2018-02-18T07:43:28.130000Z
            window_weight: 7.267822403942347


Adjoint Sources
~~~~~~~~~~~~~~~~~~~~~~~~~

Adjoint sources can be accessed in the same manner as misfit windows, through
the ``AdjointSources`` attribute of auxiliary data.

.. code:: python

    ds.auxiliary_data.AdjointSources

During an inversion, adjoint sources are tagged by the ``iter_tag`` and
``step_tag`` attributes of :class:`Config <pyatoa.core.config.Config>`

.. code:: python

    >>> ds.auxiliary_data.AdjointSources.i01.s00
    3 auxiliary data item(s) of type 'AdjointSources/i01/s00' available:
        NZ_BFZ_BXE
        NZ_BFZ_BXN
        NZ_BFZ_BXZ

Adjoint sources are stored as dictionaries with relevant creation information:

.. code:: python

    >>> ds.auxiliary_data.AdjointSources.default.NZ_BFZ_BXE
    Auxiliary Data of Type 'AdjointSources'
        Path: 'i01/s00/NZ_BFZ_BXE'
        Data shape: '(5000, 2)', dtype: 'float64'
        Parameters:
            adj_src_type: cc_traveltime_misfit
            component: BXE
            dt: 0.0725
            location:
            max_period: 20.0
            min_period: 8.0
            misfit: 24.220799999999993
            network: NZ
            starttime: 2018-02-18T07:43:28.130000Z
            station: BFZ


The actual data array of the adjoint source is also stored here in two column
format (time, amplitude):

.. code:: python

    >>> ds.auxiliary_data.AdjointSources.i01.s00.NZ_BFZ_BXE.data[:]
    array([[-20.    ,   0.    ],
           [-19.9275,   0.    ],
           [-19.855 ,   0.    ],
           ...,
           [342.2825,   0.    ],
           [342.355 ,   0.    ],
           [342.4275,   0.    ]])

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

Users can access :class:`Config <pyatoa.core.config.Config>` parameters from
the auxiliary data attribute. This is useful for understanding how windows
and adjoint sources were generated.

.. code:: python

    >>> ds.auxiliary_data.Configs.i01.s00
    Auxiliary Data of Type 'Configs'
        Path: 'i01/s00'
        Data shape: '(1,)', dtype: 'bool'
        Parameters:
            _synthetic_tag: None
            adj_src_type: cc_traveltime_misfit
            client: None
            component_list: ['Z' 'N' 'E']
            end_pad: 350
            event_id: 2018p130600
            filter_corners: 4
            iteration: 1
            max_period: 20.0
            min_period: 8.0
            observed_tag: observed
            ...


Loading Data From a Dataset
----------------------------

Data previously saved to an ``ASDFDataSet`` can be loaded back into a
:class:`Manager <pyatoa.core.manager.Manager>` class using the the
:meth:`load <pyatoa.core.manager.Manager.load>` function. This is useful for
repeating measurements, re-using misfit windows on new data, or running
seismic inversions.

Config Parameters
~~~~~~~~~~~~~~~~~

To load the :class:`Config <pyatoa.core.config.Config>` class from an
ASDFDataSet, you need to specify a ``path`` which was generated from the
``iter_tag`` and ``step_tag`` attributes of the saved
:class:`Config <pyatoa.core.config.Config>`.

.. code:: python

    cfg = Config()
    cfg.read(read_from=ds, path="i01/s00", fmt="asdf")

Data and Metadata
~~~~~~~~~~~~~~~~~~

The Managers :meth:`load <pyatoa.core.manager.Manager.load>` function searches
for metadata, waveforms and configuration parameters, based on the ``code``
and ``path`` arguments.

The ``path`` attribute is specified by the ``iter_tag`` and
``step_tag`` attributes of the saved
:class:`Config <pyatoa.core.config.Config>`.


.. note::

    Waveforms stored in the ASDFDataSet are **unprocessed**. Users will have
    to re-run the :meth:`standardize <pyatoa.core.manager.Manager.standardize>`
    and :meth:`preprocess <pyatoa.core.manager.Manager.preprocess>` functions
    to retrieve the waveforms used to generate saved windows/adjoint sources.

.. code:: python

    mgmt = Manager(ds=ds)
    mgmt.load(code="NZ.BFZ", path="i01/s00")


Windows and Adjoint Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    Misfit windows and adjoint sources are **not** explicitely re-loaded when
    calling the load function.

To re-load windows, you can call the
:meth:`window <pyatoa.core.manager.Manager.window>` function, setting the
``fix_windows`` argument to True and specifying the ``iteration`` and
``step_count`` to retrieve windows from:

.. code:: python

    mgmt.window(fix_windows=True, iteration="i01", step_count="s00")

The Manager does not currently have the capability to re-load adjoint sources,
but given a loaded Config and set of windows, you can re-calculate adjoint
sources with the :meth:`measure <pyatoa.core.manager.Manager.measure>` function:

.. code:: python

    mgmt.measure()


