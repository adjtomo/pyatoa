Config
======

The :doc:`Config </modules/core.config>` class controls the internal workflow and structure of Pyatoa, and is accessed by almost all of the core classes. Configuration parameters are used to determine how waveforms are gathered, processed, windowed and measured. There are additional bookkeeping parameters to ensure that all data is maintained to the same standard throughout an inversion. Configs can be saved to text files, or into ASDFDataSets, as a form of provenance.

--------------

Initialization
--------------

An empty configuration class comes with some preset values that are
acceptable for a long-period regional seismic inversion. These
parameters will probably not satisfy use-cases outside of this scenario,
but provide a template starting point for future adjustments.

.. code:: ipython3

    from pyatoa import Config
    cfg = Config()

.. code:: ipython3

    cfg




.. parsed-literal::

    CONFIG
        iteration:               None
        step_count:              None
        event_id:                None
    GATHER
        client:                  None
        start_pad:               20
        end_pad:                 500
        save_to_ds:              True
    PROCESS
        min_period:              10.0
        max_period:              30.0
        filter_corners:          2.0
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
        pyflex_config:           <pyflex.config.Config object at 0x7f8c007d6250>
        pyadjoint_config:        <pyadjoint.config.Config object at 0x7f8c007d64d0>



File naming convention
~~~~~~~~~~~~~~~~~~~~~~

The ``iteration`` and ``step_count`` parameters are used for internal
naming. They can be set using integer values or strings. Formatted tags
can be accessed using the ``iter_tag`` and ``step_tag`` parameters.

   **NOTE**: The formatted convetion for iterations is *i??*. For step
   counts it is *s??*. (``?`` takes the place of a single integer from
   0-9). Iterations start from 1, step counts start from 0.

The parameter ``synthetic_tag``, used to save synthetic waveforms,
automatically reflects changes to the ``iteration`` and ``step_count``
variables.

.. code:: ipython3

    cfg = Config(iteration=1, step_count=0)
    print(f"{cfg.iter_tag} == {cfg.iteration}")
    print(f"{cfg.step_tag} == {cfg.step_count}")
    print(cfg.synthetic_tag)


.. parsed-literal::

    i01 == 1
    s00 == 0
    synthetic_i01s00


.. code:: ipython3

    cfg.iteration = 2
    cfg.step_count = 3
    
    print(f"{cfg.iter_tag} == {cfg.iteration}")
    print(f"{cfg.step_tag} == {cfg.step_count}")
    print(cfg.synthetic_tag)


.. parsed-literal::

    i02 == 2
    s03 == 3
    synthetic_i02s03


--------------

External Configurations
-----------------------

The ``Config`` class also contains `Pyflex
Config <http://krischer.github.io/pyflex/#config-object>`__ and
`Pyadjoint
Config <https://github.com/krischer/pyadjoint/blob/master/src/pyadjoint/config.py>`__
objects. Preset parameters can be defined using the ``pyflex_preset``
and ``adj_src_type`` parameters, which take map names and converts them
into a set of parameters. Alternatively Pyflex Config and Pyadjoint
Config keyword arguments can be passed directly to the ``Pyatoa.Config``
class.

For specific arguments of the Pyflex and Pyadjoint Config parameters,
see their respective documentation pages. For available choices of
``pyflex_preset`` and ``adj_src_type``, see the following API.

.. code:: ipython3

    # A few randomly chosen arguments to check
    example_pyflex_kwargs = ["s2n_limit", "c_0", "max_time_before_first_arrival"]
    
    # List the arguments for the 'example' preset
    cfg = Config(pyflex_preset="example")
    print("PYFLEX CONFIG")
    for ex in example_pyflex_kwargs:
        print(f"\t{ex}: {getattr(cfg.pyflex_config, ex)}")
    
    # Modify the arguments of the 'example' preset
    cfg = Config(pyflex_preset="example", s2n_limit=2.0, c_0=1.0, max_time_before_first_arrival=25.0)
    print("\nMODIFIED PYFLEX CONFIG")
    for ex in example_pyflex_kwargs:
        print(f"\t{ex}: {getattr(cfg.pyflex_config, ex)}")


.. parsed-literal::

    PYFLEX CONFIG
    	s2n_limit: 1.5
    	c_0: 0.7
    	max_time_before_first_arrival: 50.0
    
    MODIFIED PYFLEX CONFIG
    	s2n_limit: 2.0
    	c_0: 1.0
    	max_time_before_first_arrival: 25.0


--------------

Reading / Writing
-----------------

The ``Config`` class can be read to and written from YAML files and
ASDFDataSets. This is accomplished using the ``read`` and ``write``
functions. This is handy if a specific suite of configuration parameters
will need to be accessed in the future, as may happen in an inversion
workflow. Here we show this capability using a YAML file. ASDFDataSet
capabilites are showcased in the ``storage`` documentation page.

.. code:: ipython3

    # Reading and writing from a yaml file
    cfg = Config(min_period=12.345)
    cfg.write(write_to="../tests/test_data/docs_data/test_config", fmt="yaml")
    cfg_check = Config(yaml_fid="../tests/test_data/docs_data/test_config.yaml")
    print(cfg_check.min_period)


.. parsed-literal::

    12.345

