Windowing Parameters
====================

`Pyflex <adjtomo.github.io/pyflex>`__ is a Python port of the FLEXWIN
algorithm for automatically selecting windows between data and synthetic 
waveforms for the purposes of seismic tomography.

Below are a number of 'presets' for Pyflex which have been used for various
tomography studies. They are provided as reference and starting points for
determining the optimal set of windowing parameters for your own study.

A small description is provided for each Pyflex preset. Only values which
differ from the default values are provided.

Descriptions of a few commonly used parameters that are not self explanatory:

    1. Short Term Average Long Term Average
        - ``stalta_waterlevel``: reject windows where sta/lta waveform dips
          below this threshold value. Values between 0 and 1
    2. Water level rejection
        - ``c_0`` reject if window.stalta.min() < c_0 * stalta_waterlevel
        - ``c_1`` min_acceptable_window_length = c_1 * T_min
    3. Prominence rejection
        - ``c_2`` reject if window.stalta.min() < c_2 * window.stalta.max()
    4. Separation height in phase separation
        - ``c_3a`` `d_stlta > c_3a * d_stalta_center * f_time`
          where *d_stalta* is the current max height above min and
          *d_stalta_center* is the central max height above min and
          *f_time* is the time decay function
        - ``c_3b`` d_time = separation between center of window and internal
          maxima if `d_time > c_3b` then *f_time* is a time decay function else
          it is 1 if `c_3b` goes down
    5. Emergent start/stops and coda wave curtailing
        - ``c_4a``: *time_decay_left = T_min * c_4a / dt*
        - ``c_4b``: *time_decay_right: T_min * c_4b / dt*

Default
~~~~~~~
The default windowing parameters set in ``pyflex.core.config.Config()``. All the
remaining presets override the parameter values set here.

.. code:: yaml

    stalta_waterlevel: 0.07
    tshift_acceptance_level: 10.0
    tshift_reference: 0.0
    dlna_acceptance_level: 1.3
    dlna_reference: 0.0
    cc_acceptance_level: 0.7
    s2n_limit: 1.5
    earth_model: "ak135"
    min_surface_wave_velocity: 3.0
    max_time_before_first_arrival: 50.0
    c_0: 1.0
    c_1: 1.5
    c_2: 0.0
    c_3a: 4.0
    c_3b: 2.5
    c_4a: 2.0
    c_4b: 6.0
    check_global_data_quality: False
    snr_integrate_base: 3.5
    snr_max_base: 3.0
    noise_start_index: 0
    noise_end_index: None
    signal_start_index: None
    signal_end_index: -1
    window_weight_fct: None
    window_signal_to_noise_type: "amplitude"
    resolution_strategy: "interval_scheduling"

-----------------------

Homogeneous Halfspace
~~~~~~~~~~~~~~~~~~~~~
Very relaxed windowing parameters that should work with simple synthetics
generated using a homogeneous halfspace 

.. code:: yaml

    stalta_waterlevel: 0.05
    tshift_acceptance_level: 15.0
    dlna_acceptance_level: 2.0
    s2n_limit: 3.
    min_surface_wave_velocity: 1.5
    c_0: 0.7
    c_1: 2.
    c_3a: 3.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0

Maggi et al. (2009)
~~~~~~~~~~~~~~~~~~~
The following windowing parameters are defined in Table 3 of Maggi et al. (2009)
for various regions globally.

Global (50--150s)
`````````````````

.. code:: yaml

    s2n_limit: 2.5
    stalta_waterlevel: 0.08
    cc_acceptance_level: 0.85
    tshift_acceptance_level: 15.0
    dlna_acceptance_level: 1.
    tshift_reference: 0.
    dlna_reference: 0.
    c_0: 0.7
    c_1: 4.
    c_2: 0.3
    c_3a: 1.0
    c_3b: 2.0
    c_4a: 3.
    c_4b: 10.0

Japan (6--30s)
`````````````````

.. code:: yaml

    s2n_limit: 3.
    stalta_waterlevel: 0.12
    cc_acceptance_level: 0.73
    tshift_acceptance_level: 3.0
    dlna_acceptance_level: 1.5
    tshift_reference: 0.
    dlna_reference: 0.
    c_0: 0.7
    c_1: 3.
    c_2: 0.6
    c_3a: 1.0
    c_3b: 2.0
    c_4a: 3.
    c_4b: 12.0

Southern California (6--30s)
``````````````````````````````
    
.. code:: yaml

    s2n_limit: 3.
    stalta_waterlevel: 0.18
    cc_acceptance_level: 0.71
    tshift_acceptance_level: 8.0
    dlna_acceptance_level: 1.5
    tshift_reference: 4.
    dlna_reference: 0.
    c_0: 0.7
    c_1: 2.
    c_2: 0.
    c_3a: 3.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0

Southern California (3--30s)
````````````````````````````

.. code:: yaml

    s2n_limit: 4.
    stalta_waterlevel: 0.11
    cc_acceptance_level: 0.8
    tshift_acceptance_level: 4.0
    dlna_acceptance_level: 1.
    tshift_reference: 2.
    dlna_reference: 0.
    c_0: 1.3
    c_1: 4.
    c_2: 0.
    c_3a: 4.0
    c_3b: 2.5
    c_4a: 2.
    c_4b: 6.0

Southern California (2--30s)
````````````````````````````

.. code:: yaml

    s2n_limit: 4.
    stalta_waterlevel: 0.07
    cc_acceptance_level: 0.85
    tshift_acceptance_level: 3.0
    dlna_acceptance_level: 1.
    tshift_reference: 1.
    dlna_reference: 0.
    c_0: 1.
    c_1: 5.
    c_2: 0.
    c_3a: 4.0
    c_3b: 2.5
    c_4a: 2.
    c_4b: 6.0


NZATOM_NORTH
~~~~~~~~~~~~
The following parameter sets were used to derive NZATOM_NORTH an adjoint
tomography model for the North Island of New Zealand. The results of this study
are published in Chow et al. (2020) and Chow etal. (2022a).

The parameter set is split into various inversion legs which tackle different
passbands of the dataset.

LEG 1 (15--30s)
```````````````````````
.. code:: yaml

    stalta_waterlevel: 0.08
    tshift_acceptance_level: 12.0
    dlna_acceptance_level: 2.5
    cc_acceptance_level: 0.7
    s2n_limit: 2.5
    max_time_before_first_arrival: 10.
    min_surface_wave_velocity: 1.2 
    check_global_data_quality: True
    c_0: 0.7
    c_1: 2.0
    c_3a: 1.0
    c_3b: 2.0
    c_4a: 3.0
    c_4b: 10.0

LEG 2 (10--30s)
```````````````````````
.. code:: yaml

    stalta_waterlevel: 0.10
    tshift_acceptance_level: 8.0  # based on sign-flip
    dlna_acceptance_level: 2.0
    cc_acceptance_level: 0.7
    s2n_limit: 3.
    max_time_before_first_arrival: 5.
    min_surface_wave_velocity: 1.2
    check_global_data_quality: True
    c_0: 0.7
    c_1: 2.0
    c_3a: 3.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0

LEG 3 (8--30s)
``````````````````````
.. code:: yaml

    stalta_waterlevel: 0.10 
    tshift_acceptance_level: 8.0
    dlna_acceptance_level: 1.5
    cc_acceptance_level: 0.7
    s2n_limit: 3.
    max_time_before_first_arrival: 5.
    min_surface_wave_velocity: 1.1
    check_global_data_quality: True
    c_0: 0.7
    c_1: 2.0  # min window = c1 * tmin = 16s 
    c_3a: 4.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0

LEG 4 (6--30s)
``````````````````````
.. code:: yaml

    stalta_waterlevel: 0.08
    tshift_acceptance_level: 8.  
    dlna_acceptance_level: 1.5
    cc_acceptance_level: 0.60
    s2n_limit: 3.
    max_time_before_first_arrival: 5. 
    min_surface_wave_velocity: 1.05
    check_global_data_quality: True
    snr_integrate_base: 3.5  # exclude noisy data
    c_0: 0.8     # reject if win.stalta.min < c_0 * stalta_wl
    c_1: 2.0     # min window = c1 * tmin = 12s
    c_3a: 3.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0

LEG 5 (4--30s)
``````````````````````
.. code:: yaml

    stalta_waterlevel: 0.075
    tshift_acceptance_level: 6.
    dlna_acceptance_level: 1.5
    cc_acceptance_level: 0.65
    s2n_limit: 4.
    max_time_before_first_arrival: 5.
    min_surface_wave_velocity: 1.0
    check_global_data_quality: True
    snr_integrate_base: 3.5  # exclude noisy data
    c_0: 0.9     # reject if win.stalta.min < c_0 * stalta_wl
    c_1: 3.
    c_3a: 3.5
    c_3b: 2.25
    c_4a: 2.25
    c_4b: 9.0

POSTHOC (6--30s) 
````````````````````````
This was used for posthoc evaluation of the final model using events not
inverted for during the inversion.

.. code:: yaml

    stalta_waterlevel: 0.08
    tshift_acceptance_level: 12.  
    dlna_acceptance_level: 1.5
    cc_acceptance_level: 0.60
    s2n_limit: 3.
    max_time_before_first_arrival: 5. 
    min_surface_wave_velocity: 1.05
    check_global_data_quality: True
    snr_integrate_base: 3.5  # exclude noisy data
    c_0: 0.8     # reject if win.stalta.min < c_0 * stalta_wl
    c_1: 2.0     # min window = c1 * tmin = 12s
    c_3a: 3.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0

RISTAU 1D (10--30s) 
```````````````````````````
Used for windowing synthetic waveforms generated using the 1D model of the
North Island of New Zealand from Ristau (2008)

.. code:: yaml

    stalta_waterlevel: 0.10
    tshift_acceptance_level: 120  # based on sign-flip
    dlna_acceptance_level: 20
    cc_acceptance_level: 0675
    s2n_limit: 3
    max_time_before_first_arrival: 5
    min_surface_wave_velocity: 16
    check_global_data_quality: True
    c_0: 07
    c_1: 20
    c_3a: 30
    c_3b: 20
    c_4a: 25
    c_4b: 120

RISTAU 1D (8--30s) 
```````````````````````````

.. code:: yaml

    stalta_waterlevel: 0.08 
    tshift_acceptance_level: 10.0
    dlna_acceptance_level: 2.0
    cc_acceptance_level: 0.675
    s2n_limit: 3.
    max_time_before_first_arrival: 5.
    min_surface_wave_velocity: 1.4
    check_global_data_quality: True
    c_0: 0.7
    c_1: 2.5 
    c_3a: 3.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0

ALASKA
~~~~~~
Values from group at University of Alaska Fairbanks performing regional studies
of Alaska

General
```````
.. code:: yaml

    stalta_waterlevel: 0.18
    tshift_acceptance_level: 4.0
    dlna_acceptance_level: 1.5
    cc_acceptance_level: 0.71
    c_0: 0.7
    c_1: 2.0
    c_3a: 3.0
    c_3b: 2.0
    c_4a: 2.5
    c_4b: 12.0