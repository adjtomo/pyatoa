Misfit Quantification
=====================

Pyatoa is comprised of a two core classes which take care of the main package
functionality. These are the :class:`Config <pyatoa.core.config.Config>`,
:class:`Manager <pyatoa.core.manager.Manager>` classes. The following
page explains each class and their functionalities.

Logging
~~~~~~~

Pyatoa comes with a detailed logger; most processes will output log messages.
Log levels in decreasing verbosity are: DEBUG, INFO, WARNING, CRITICAL

.. code:: python

    from pyatoa import logger

    logger.setLevel("DEBUG")

An example of logging to stdout during misfit quantification:

.. code-block:: text

    [2022-03-03 11:01:32] - pyatoa - INFO: preprocessing observation data
    [2022-03-03 11:01:32] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-03 11:01:32] - pyatoa - DEBUG: removing response, units to DISP
    [2022-03-03 11:01:32] - pyatoa - DEBUG: rotating from generic coordinate system to ZNE
    [2022-03-03 11:01:32] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-03 11:01:32] - pyatoa - INFO: preprocessing synthetic data
    [2022-03-03 11:01:32] - pyatoa - INFO: adjusting taper to cover time offset -20.0
    [2022-03-03 11:01:32] - pyatoa - DEBUG: no response removal, synthetic data or requested not to
    [2022-03-03 11:01:32] - pyatoa - DEBUG: bandpass filter: 10.0 - 30.0s w/ 2.0 corners
    [2022-03-03 11:01:32] - pyatoa - DEBUG: convolving data w/ Gaussian (t/2=0.70s)
    [2022-03-03 11:01:32] - pyatoa - INFO: running Pyflex w/ map: default

Config
~~~~~~

The :class:`Config <pyatoa.core.config.Config>` class is how Users set
parameters in Pyatoa. Config parameters determine how waveforms are gathered,
processed, windowed and measured.

Internal bookkeeping parameters ensure that all data is maintained to the same
standard during an inversion.

Printing the Config instance will provide a useful display of Config
parameters.

.. code:: python

    >>> from pyatoa import Config
    >>> cfg = Config(min_period=10., max_period=30.)
    >>> print(cfg)
    CONFIG
        iteration:               None
        step_count:              None
        event_id:                None
    GATHER
        start_pad:               20
        end_pad:                 500
        save_to_ds:              True
    PROCESS
        min_period:              10.0
        max_period:              30.0
        filter_corners:          2
        unit_output:             DISP
        rotate_to_rtz:           False
        win_amp_ratio:           0.0
        synthetics_only:         False
    LABELS
        component_list:          ['E', 'N', 'Z']
        observed_tag:            observed
        synthetic_tag:           synthetic
        paths:                   {'waveforms': [], 'synthetics': [], 'responses': [], 'events': []}
    EXTERNAL
        pyflex_preset:           default
        adj_src_type:            cc_traveltime_misfit
        pyflex_config:           <pyflex.config.Config object at 0x142c4fbd0>
        pyadjoint_config:        <pyadjoint.config.Config object at 0x135e16b10>


Reading/Writing Config
``````````````````````

A set of :class:`Config <pyatoa.core.config.Config>` parameters can be read to
and written from YAML files and ASDFDataSets. To write the Config to a YAML
file:

.. code:: python

    cfg.write(write_to="config.yaml", fmt="yaml")

See the `Saving data with ASDF <storage.html>`__ doc page to see how
the Config object is written to ASDFDataSets.


File Naming Convention
``````````````````````
The :class:`Config <pyatoa.core.config.Config>` object includes parameters that
are used to keep track of files during an inversion.

Also see the `Standards <standards.html>`__ page for more details on file
naming conventions.

Iteration and Step Count
++++++++++++++++++++++++

The ``iteration`` and line search ``step_count`` parameters are used to tag
synthetic waveform data and output figures.

Users can access the string representations used to tag files through the
``iter_tag`` and ``step_tag`` attributes.

.. code:: python

    >>> cfg = Config(iteration=1, step_count=0)
    >>> print(cfg.iter_tag)
    i01
    >>> print(cfg.step_tag)
    s00

Waveform Tags
+++++++++++++

The ``observed_tag`` and ``synthetic_tag`` parameters are used to save waveforms
in ASDFDataSets. See the `Saving data with ASDF <storage.html>`__ doc page
to see how to access waveforms within an ASDFDatASet.

The `synthetic_tag` distinguishes which model they were
created with and is derived directly from the iteration and step count tags.

.. code:: python

    >>> print(cfg.observed_tag)
    observed
    >>> print(cfg.synthetic_tag)
    synthetic_i01s00


See the `standards <standards.html>`__ docs page for more information on
the standards that Pyatoa uses for internal and external file naming.

Windowing and Measurement Parameters
````````````````````````````````````

Under the hood, Config controls the
`Pyflex Config <http://adjtomo.github.io/pyflex/#config-object>`__ and
`Pyadjoint Config
<https://github.com/krischer/pyadjoint/blob/master/src/pyadjoint/config.py>`__
objects. Valid parameters of those :class:`Config <pyatoa.core.config.Config>` objects
can be passed directly to Config.

The ``pyflex_preset`` and ``adj_src_type`` parameter lets the User define the
misfit function.

- `Click here to see available Pyflex presets <https://github.com/adjtomo/pyatoa/blob/master/pyatoa/plugins/pyflex_presets.py>`__
- `Click here to see available adjoint source types <http://adjtomo.github.io/pyadjoint/adjoint_sources/index.html>`__

.. code:: python

    >>> from pyatoa import Config
    >>> cfg = Config(pyflex_preset="default",
    >>>              adj_src_type="cc_traveltime_misfit",
    >>>              tshift_acceptance_level=8.0,  # Pyflex parameter,
    >>>              min_cycle_in_window=1.0       # Pyadjoint parameter
    >>>              )
    >>> print(cfg.pyflex_config.tshift_acceptance_level)
    8.0
    >>> print(cfg.pyadjoint_config.min_cycle_in_window)
    1.0


Manager
~~~~~~~

The :class:`Manager <pyatoa.core.manager.Manager>` is the main workhorse of
Pyatoa. Its job is to group waveforms and metadata, process misfit, and output
misfit windows and adjoint sources.

The Manager takes the :class:`Config <pyatoa.core.config.Config>` object as
input, which allows the User to control internal processing. Printing the
Manager shows available data and processing status.

.. note::

    If no Config object is provided, the Manager will instantiate its own with
    default parameters.

.. code:: python

    >>> from pyatoa import Config, Manager
    >>> cfg = Config()
    >>> mgmt = Manager(config=cfg)
    Manager Data
        dataset   [ds]:        None
        quakeml   [event]:     None
        station   [inv]:       None
        observed  [st_obs]:    None
        synthetic [st_syn]:    None
    Stats & Status
        half_dur:              None
        time_offset_sec:       None
        standardized:          False
        obs_processed:         False
        syn_processed:         False
        nwin   [windows]:      None
        misfit [adjsrcs]:      None

Loading Example Data
````````````````````

To load some example data and play around with Manager, you can use the
:meth:`load <pyatoa.core.manager.Manager.load>`
function which will grab data from the local test directory.

Test data includes an event, station response, and observed and synthetic
waveforms. Printing the Manager shows the loaded data available.

.. code:: python

    >>> mgmt.load()
    >>> print(mgmt)
    Manager Data
        dataset   [ds]:        None
        quakeml   [event]:     smi:nz.org.geonet/2018p130600
        station   [inv]:       NZ.BFZ
        observed  [st_obs]:    3
        synthetic [st_syn]:    3
    Stats & Status
        half_dur:              0.6989458964552759
        time_offset_sec:       None
        standardized:          False
        obs_processed:         False
        syn_processed:         False
        nwin   [windows]:      None
        misfit [adjsrcs]:      None


The load function is also used to load previously saved data from an
ASDFDataSet. See the `Saving data with ASDF <storage.html>`__ doc page for
more information.


Providing Data
``````````````

The simplest method to provide the
:class:`Manager <pyatoa.core.manager.Manager>` with data is to set it's
attributes. Data are provided and stored as ObsPy objects.

At a minimum, Manager expects two waveforms, observed (``st_obs``) and synthetics
(``st_syn``). Despite the labels, these can be any types of waveforms (i.e.,
two synthetics; two sets of observed waveforms).

.. code:: python

    from obspy import read

    st_obs = read()
    st_syn = read()

    mgmt = Manager(st_obs=st_obs, st_syn=st_syn)


To unlock the full potential of the Manager, metadata should also be provided.
These include station metadata, including response (``inv``) and event metadata
(``event``).

.. code:: python

    from obspy import read_events, read_inventory

    event = read_events("some_example_catalog.xml")[0]
    inv = read_inventory("some_example_stationxml.xml")

    mgmt.inv = inv
    mgmt.event = event


.. warning::

    If metadata are not provided, some check criteria during the windowing and
    preprocessing will be skipped. Similarly, the Manager will not be able to
    plot a source-receiver map.

Accessing Data
``````````````

Accessing data is done by accessing the Manager's attributes. Data are stored
as ObsPy objects. Use the `print` command to determine the names of the
relevant attributes.

.. code:: python

    >>> from pyatoa import Manager
    >>> mgmt = Manager()
    >>> mgmt.load()
    >>> print(mgmt)
    Manager Data
        dataset   [ds]:        None
        quakeml   [event]:     smi:nz.org.geonet/2018p130600
        station   [inv]:       NZ.BFZ
        observed  [st_obs]:    3
        synthetic [st_syn]:    3
    Stats & Status
        half_dur:              0.6989458964552759
        time_offset_sec:       None
        standardized:          False
        obs_processed:         False
        syn_processed:         False
        nwin   [windows]:      None
        misfit [adjsrcs]:      None
    >>> print(mgmt.event)
    Event:	2018-02-18T07:43:48.127644Z | -39.949, +176.300 | 5.16 M  | manual

                          resource_id: ResourceIdentifier(id="smi:nz.org.geonet/2018p130600")
                           event_type: 'earthquake'
                        creation_info: CreationInfo(agency_id='WEL(GNS_Primary)', author='scevent@kseqp01.geonet.org.nz', creation_time=UTCDateTime(2018, 2, 18, 7, 44, 9, 156454))
                  preferred_origin_id: ResourceIdentifier(id="smi:nz.org.geonet/Origin#20180226021110.13419.62761")
               preferred_magnitude_id: ResourceIdentifier(id="smi:nz.org.geonet/Origin#20180226021110.13419.62761#netMag.M")
         preferred_focal_mechanism_id: ResourceIdentifier(id="smi:local/ad83e11b-cc91-4de7-9cd0-5c51f99e1062")
                                 ---------
                   event_descriptions: 1 Elements
                     focal_mechanisms: 1 Elements
                              origins: 1 Elements
                           magnitudes: 3 Elements
    >>> print(mgmt.inv)
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
    >>> print(mgmt.st_obs)
    3 Trace(s) in Stream:
    NZ.BFZ.10.HHE | 2018-02-18T07:43:28.128394Z - 2018-02-18T07:49:38.128394Z | 100.0 Hz, 37001 samples
    NZ.BFZ.10.HHN | 2018-02-18T07:43:28.128394Z - 2018-02-18T07:49:38.128394Z | 100.0 Hz, 37001 samples
    NZ.BFZ.10.HHZ | 2018-02-18T07:43:28.128394Z - 2018-02-18T07:49:38.128394Z | 100.0 Hz, 37001 samples
    >>> print(mgmt.st_syn)
    3 Trace(s) in Stream:
    NZ.BFZ..BXE | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples
    NZ.BFZ..BXN | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples
    NZ.BFZ..BXZ | 2018-02-18T07:43:28.127644Z - 2018-02-18T07:48:28.097644Z | 33.3 Hz, 10000 samples


Processing Functions
````````````````````

The :class:`Manager <pyatoa.core.manager.Manager>` has four main processing
functions that it applies on data and synthetics.

- :meth:`standardize <pyatoa.core.manager.Manager.standardize>`: match the time series of the data and synthetics
- :meth:`preprocess <pyatoa.core.manager.Manager.preprocess>`: remove response, detrend and filter data
- :meth:`window <pyatoa.core.manager.Manager.window>`: generate misfit windows based on preprocessed data
- :meth:`measure <pyatoa.core.manager.Manager.measure>`: calculate misfit and generate adjoint sources for given windows

Standardize Time Series
++++++++++++++++++++++++

Oftentimes, observed and synthetic waveforms will differ in sampling rate,
start and end time. The
:meth:`standardize <pyatoa.core.manager.Manager.standardize>`
function matches time series for the two waveforms: `st_obs` and `st_syn`.

.. code:: python

    mgmt.standardize(standardize_to="syn")

.. note::

    By default, Manager will standardize both time series' to the synthetic
    trace, as it is assumed that the adjoint source resulting from the
    processing will require the same time array as the synthetics.

Preprocess Waveforms
+++++++++++++++++++++

The :meth:`preprocess <pyatoa.core.manager.Manager.preprocess>` function
involves detrending and filtering, with additional instrument response removal
for observed waveforms.

.. code:: python

    mgmt.preprocess(which="both")

.. note::

    By default, Manager will preprocess both `st_obs` and `st_syn`. Users can
    choose selectively with the `which` parameter.

Custom Preprocessing Scripts
.............................

Pyatoa has a default preprocessing script which it applies to both observed and
synthetic data. Some users may wish to use their own preprocessing function.
This can be achieved using the ``overwrite`` command.

.. code:: python

    def custom_preprocessing(st, choice, inv, unit_output="disp", **kwargs):
        """
        This function performs a custom preprocessing for the Manager class.

        :type st: obspy.core.stream.Stream
        :param st: Stream object to preprocess
        :type choice: str
        :param choice: choice of output, either "obs" or "syn"
        :rtype: obspy.core.stream.Stream
        :return: A preprocessed ObsPy Stream object
        """
        # The `choice` argument for different processing of `st_obs`, `st_syn`
        if choice == "obs":
            st.remove_response(inventory=inv, output=unit_output)

            # Here we add a random action to scale data
            for tr in st:
                tr.data *= 2

        # Access to Config parameters is still possible
        st.filter("bandpass", freqmin=1/max_period, freqmax=1/min_period)

        # MUST output a Stream
        return st

    mgmt.preprocess(overwrite=custom_preprocessing)


Generate Misfit Windows
++++++++++++++++++++++++

Pyatoa uses Pyflex to window observed and synthetic waveforms. Windowing
parameters are stored in ``Config.pyflex_config`` and is set internally via
the function :meth:`pyatoa.core.config.set_pyflex_config`

Under the hood, the :meth:`window <pyatoa.core.manager.Manager.window>` calls
the Pylex package to generate misfit windows for the two waveforms ``st_obs``
and ``st_syn``.


.. code:: python

    mgmt.window()

To access created misfit windows, check the `windows` attribute.

Have a look at the `Saving data with ASDF <storage.html>`__ doc page to
see how misfit windows are stored in ASDFDataSets.

.. code:: python

    >>> mgmt.windows
    {'E': [Window(left=990, right=3187, center=2088, channel_id=NZ.BFZ.10.HHE, max_cc_value=0.8899832380628487, cc_shift=31, dlnA=-0.6397761404459611)],
     'N': [Window(left=941, right=3006, center=1973, channel_id=NZ.BFZ.10.HHN, max_cc_value=0.9753605590906922, cc_shift=63, dlnA=-0.8370414140149721)]}


The total number of collected windows is stored in the `stats` attribute:

.. code:: python

    >>> mgmt.stats.nwin
    2


Rejected time windows, useful for plotting or to aid in fine-tuning of the
windowing algorithm can be accessed in the `rejwins` attribute.

.. code:: python

    >>> mgmt.rejwins
    {'E': {'water_level': [Window(left=990, right=3787, center=1499, channel_id=NZ.BFZ.10.HHE, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=990, right=4062, center=1499, channel_id=NZ.BFZ.10.HHE, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=990, right=4558, center=1499, channel_id=NZ.BFZ.10.HHE, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=990, right=4741, center=1499, channel_id=NZ.BFZ.10.HHE, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=990, right=4902, center=1499, channel_id=NZ.BFZ.10.HHE, max_cc_value=None, cc_shift=None, dlnA=None)]},
     'N': {'water_level': [Window(left=941, right=3387, center=1445, channel_id=NZ.BFZ.10.HHN, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=941, right=3789, center=1445, channel_id=NZ.BFZ.10.HHN, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=941, right=4170, center=1445, channel_id=NZ.BFZ.10.HHN, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=941, right=4902, center=1445, channel_id=NZ.BFZ.10.HHN, max_cc_value=None, cc_shift=None, dlnA=None)]},
     'Z': {'min_length': [Window(left=909, right=1377, center=1291, channel_id=NZ.BFZ.10.HHZ, max_cc_value=None, cc_shift=None, dlnA=None)],
      'water_level': [Window(left=909, right=4902, center=1291, channel_id=NZ.BFZ.10.HHZ, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=909, right=4902, center=1716, channel_id=NZ.BFZ.10.HHZ, max_cc_value=None, cc_shift=None, dlnA=None),
       Window(left=1377, right=4902, center=1716, channel_id=NZ.BFZ.10.HHZ, max_cc_value=None, cc_shift=None, dlnA=None)],
      'dlna': [Window(left=909, right=3800, center=2354, channel_id=NZ.BFZ.10.HHZ, max_cc_value=0.9101699760469617, cc_shift=85, dlnA=-1.3118583378421544),
       Window(left=1377, right=3800, center=2588, channel_id=NZ.BFZ.10.HHZ, max_cc_value=0.9267908457090609, cc_shift=86, dlnA=-1.3247299121081957)]}}

Fixed Time Windows
...................

Users can use a previously generated set of time windows to evaluate
misfit on new waveforms. Rather than select new windows, the Manager can load
a previous set of windows from an ASDFDataSet.

The :class:`Config <pyatoa.core.config.Config>` parameters ``iteration`` and
``step_count`` are important here, as they are used to tag saved windows and
load them at a later time.

.. code:: python

    from pyasdf import ASDFDataSet as asdf
    from pyatoa import Config, Manager

    # Load in dataset that has saved misfit windows
    ds = ASDFDataSet("test_dataset.h5")

    mgmt = Manager(ds=ds, config=cfg)
    mgmt.load()  # some example data, this could be any data

    mgmt = Manager(ds=ds)
    mgmt.standardize().preprocess()  # it is possible to chain functions

    # Load in previously saved windows
    mgmt.window(fix_windows=True, iteration="i01", step_count="s00")

Generate Adjoint Sources
+++++++++++++++++++++++++

Manager uses Pyadjoint to measure misfit within time windows, and generate
adjoint sources for a seismic inversion. The type of adjoint source is defined
by ``Config.adj_src_type`` and parameters are set internally with the function
:meth:`pyatoa.core.config.set_pyadjoint_config`.


The :meth:`measure <pyatoa.core.manager.Manager.measure>` function calls
Pyadjoint under the hood to generate an adjoint source within the time windows
selected by the :meth:`window <pyatoa.core.manager.Manager.window>` function.

.. note::

    If no windows are provided or calculated, the Manager will calcualte misfit
    along the entire time series

.. code:: python

    mgmt.measure()

To access the generated adjoint sources, check the `adjsrcs` attribute:

.. code:: python

    >>> mgmt.adjsrcs
    {'E': <pyadjoint.adjoint_source.AdjointSource at 0x104999a10>,
     'N': <pyadjoint.adjoint_source.AdjointSource at 0x132c354d0>}
    >>> vars(mgmt.adjsrcs["E"])
    {'adj_src_type': 'cc_traveltime_misfit',
     'adj_src_name': 'Cross Correlation Traveltime Misfit',
     'misfit': 0.30411925696681014,
     'dt': 0.03,
     'min_period': 10,
     'max_period': 100,
     'component': 'BXE',
     'network': 'NZ',
     'station': 'BFZ',
     'location': '10',
     'starttime': 2018-02-18T07:43:28.127644Z,
     'adjoint_source': array([0., 0., 0., ..., 0., 0., 0.])}
    >>> mgmt.adjsrcs["E"].adjoint_source
    array([0., 0., 0., ..., 0., 0., 0.])

Misfit information is stored in the `stats` attribute:

.. code:: python

    >>> mgmt.stats.misfit
    2.09016925696681


Plotting
+++++++++

The Manager has built-in plotting functions to plot waveforms, misfit windows
adjoint sources and a source receiver map.

To plot waveforms and map in the same figure (done by default),

.. code:: python

    mgmt.plot(choice="both")

Otherwise Users can plot the waveforms on their own

.. code:: python

    mgmt.plot(choice="wav")

Or the map on its own

.. code:: python

    mgmt.plot(choice="map")


Flow Function
++++++++++++++

The :meth:`flow <pyatoa.core.manager.Manager.flow>` function simply chains all
the preprocessing steps together. It is equivalent to running standardize,
preprocess, window and measure one after another.

.. code:: python

    mgmt.flow()




