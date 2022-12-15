Misfit Quantification
=====================

Pyatoa is comprised of a two core classes which take care of the main package
functionality. These are the Config, Manager classes. The following
page explains each class and their functionalities.

Config
~~~~~~

The `Config` class is how Users set parameters in Pyatoa. Config parameters
determine how waveforms are gathered, processed, windowed and measured.

Internal bookkeeping parameters ensure that all data is maintained to the same
standard during an inversion.

.. code:: python

    from pyatoa import Config

    cfg = Config(min_period=10., max_period=30.)


Printing the Config instance will provide a useful display of Config
parameters.

.. code:: python

    print(cfg)


Reading/Writing Config
``````````````````````

A set of Config parameters can be read to and written from YAML files and
ASDFDataSets. To write the Config to a YAML file:

.. code:: python

    cfg.write(write_to="config.yaml", fmt="yaml")


See the `'Storage with ASDFDataSets' <storage.html>`__ doc page to see how
the Config object is written to ASDFDataSets.


File Naming Convention
``````````````````````
The Config object includes parameters that are used to keep track of files
during an inversion.

Iteration and Step Count
++++++++++++++++++++++++

The `iteration` and line search `step_count` parameters are used to tag
synthetic waveform data and output figures.

.. code:: python

    cfg = Config(iteration=1, step_count=0)


Users can access the string representations used to tag files with:

.. code:: python

    cfg.iter_tag  # i01
    cfg.step_tag  # s00


Synthetic Tag
+++++++++++++

The `synthetic_tag` parameter is used to save synthetic waveforms to
ASDFDataSets by distinguishing which model they were created with. This tag is
derived directly from the iteration and step count tags.

.. code:: python

    cfg.synthetic_tag  # synthetic_i01s00


See the `standards <standards.html>`__ docs page for more information on
the standards that Pyatoa uses for internal and external file naming.

Windowing and Measurement Parameters
````````````````````````````````````

Under the hood, Config controls the
`Pyflex Config <http://adjtomo.github.io/pyflex/#config-object>`__ and
`Pyadjoint Config
<https://github.com/krischer/pyadjoint/blob/master/src/pyadjoint/config.py>`__
objects.

Valid parameters of those Config objects can be passed directly to Config.
The `pyflex_preset` and `adj_src_type` parameter lets the User define the
misfit function.

Click for available `pyflex_preset <https://github.com/adjtomo/pyatoa/blob/master/pyatoa/plugins/pyflex_presets.py>`__
Click for available `adj_src_types <http://adjtomo.github.io/pyadjoint/adjoint_sources/index.html>`__

.. code:: python

    from pyatoa import Config

    cfg = Config(pyflex_preset="default",
                 adj_src_type="cc_traveltime_misfit",
                 tshift_acceptance_level=8.0,  # Pyflex parameter,
                 min_cycle_in_window=1.0       # Pyadjoint parameter
                 )


The underlying Pyflex and Pyadjoint Configs can be accessed as attributes:

.. code:: python

    cfg.pyflex_config
    cfg.pyadjoint_config


Manager
~~~~~~~

The `Manager` is the main workhorse of Pyatoa. Its job is to group waveforms
and metadata, process misfit, and output misfit windows and adjoint sources.
The Manager takes `Config` as input, which controls internal processing.

If no Config object is provided, the Manager will instantiate its own with
default parameters.

.. code:: python

    from pyatoa import Config, Manager

    cfg = Config()
    mgmt = Manager(config=cfg)


Loading Example Data
````````````````````

To load some example data and play around with Manager, you can use the load
function.

.. code:: python

    mgmt.load()

The load function is also used to load previously saved data from an
ASDFDataSet. See the `'Storage with ASDFDataSets' <storage.html>`__ doc page for
more information.

Providing Data
``````````````

The simplest method to provide the Manager with data is to set it's attributes.
Data must be provided as ObsPy objects.

At a minimum, Manager expects two waveforms, observed (`st_obs`) and synthetics
(`st_syn`). Despite the labels, these can be any types of waveforms (i.e.,
two synthetics; two sets of observed waeveforms).

.. code:: python

    from obspy import read

    st_obs = read("some_example_waveform_data.mseed")
    st_syn = read("some_example_synthetic_data.mseed")

    mgmt = Manager(st_obs=st_obs, st_syn=st_syn)


To unlock the full potential of the Manager, metadata should also be provided.
These include station metadata, inlcuding response ('inv') and event metadata
('event')

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

Processing Functions
````````````````````

The Manager has four main processing functions that it applies on data and
synthetics.

- standardize: match the time series of the data and synthetics
- preprocess: remove response, detrend and filter data
- window: generate misfit windows based on preprocessed data
- measure: calculate misfit and generate adjoint sources for given windows

Standardize
+++++++++++

Oftentimes, observed and synthetic waveforms will differ in sampling rate,
start and end time. Standardize matches time series for `st_obs` and `st_syn`.

.. code:: python

    mgmt.standardize(standardize_to="syn")


.. note::

    By default, Manager will standardize both time series' to the synthetic
    trace, as it is assumed that the adjoint source resulting from the
    processing will require the same time array as the synthetics.

Preprocess
++++++++++

Preprocessing involves detrending and filtering, with additional instrument
response removal for observed waveforms.

.. code:: python

    mgmt.preprocess(which="both")


.. note::

    By default, Manager will preprocess both `st_obs` and `st_syn`. Users can
    choose selectively with the `which` parameter.

Custom Preprocessing Scripts
.............................

Pyatoa has a default preprocessing script which it applies to both observed and
synthetic data. Some users may wish to use their own preprocessing function.
This can be achieved using the `overwrite` command.

.. code:: python

    def custom_preprocessing(mgmt, choice):
        """
        This function performs a custom preprocessing for the Manager class.

        :type mgmt: pyatoa.core.manager.Manager
        :param mgmt: the Manager class, which contains standardized data
        :type choice: str
        :param choice: choice of output, either "obs" or "syn"
        :rtype: obspy.core.stream.Stream
        :return: A preprocessed ObsPy Stream object
        """
        if choice == "obs":
            st = mgmt.st_obs
        elif choice == "syn":
            st = mgmt.st_syn

        # The `choice` argument allows different preprocessing for `obs` and `syn`
        if choice == "obs":
            st.remove_response(inventory=mgmt.inv,
                               output=mgmt.config.unit_output)

            # Here we add a random action to scale data
            for tr in st:
                tr.data *= 2

        # Access to Config parameters is still possible
        st.filter("bandpass", freqmin=1/mgmt.config.max_period,
                  freqmax=1/mgmt.config.min_period)

        # MUST output a Stream
        return st

    mgmt.preprocess(overwrite=custom_preprocessing)


Window
++++++

Pyatoa uses Pyflex to window observed and synthetic waveforms. Windowing
parameters are stored in `Config.pyflex_config`.

.. code:: python

    mgmt.window()



Fixed Time Windows
...................

Pyatoa has the ability to use a previous set of time windows to evaluate
misfit. That is, rather than select new windows, the Manager will load
a previous set of windows from an ASDFDataSet.

The Config parameters `iteration` and `step_count` are important here, as they '
are used to tag saved windows and load them up at a later time.

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


To access created misfit windows, check the `windows` attribute

.. code:: python

    mgmt.windows


The total number of collected windows is stored in the `stats` attribute:

.. code:: python

    mgmt.stats.nwin


Rejected time windows, useful for plotting or to aid in fine-tuning of the
windowing algorithm can be accessed in the `rejwins` attribute

.. code:: python

    mgmt.rejwins


Measure
+++++++

Manager uses Pyadjoint to measure misfit within time windows, and generate
adjoint sources for a seismic inversion. The type of adjoint source is defined
by `Config.adj_src_type`.

.. note::

    If no windows (Manager.windows) are provided or calculated, Manager will
    calcualte misfit along the entire time series

.. code:: python

    mgmt.measure()


To access the generated adjoint sources, check the `adjsrcs` attribute:

.. code:: python

    mgmt.adjsrcs


Misfit information is stored in the `stats` attribute:

.. code:: python

    mgmt.stats.misfit


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

The Flow function simply chains all the preprocessing steps together. It is
equivalent to running standardize, preprocess, window and measure one after
another.

.. code:: python

    mgmt.flow()




