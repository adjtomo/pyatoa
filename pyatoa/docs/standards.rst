Standards
=========

Pyatoa enforces a few standards for book keeping. Many of these
standards are encapsulated in the functions found in
:module:``~pyatoa.utils.form``.

--------------

Receiver codes
--------------

Receivers can be identified multiple ways in the workflow. Most of the
time when Pyatoa asks the user to identify a receiver, it will be
through its **code**:

-  **CODE**: The full SEED standard code **NN.SSS.LL.CCC** (where
   **N**\ =Network, **S**\ =Station, **L**\ =Location, **C**\ =Channel);
   e.g. \ `NZ.BFZ.10.HHZ <https://www.geonet.org.nz/data/network/sensor/BFZ>`__

..

   **> wildcards**: Wildcards are accepted and encouraged when
   specifying full codes. For example, in data gathering, selecting all
   available components is accomplished by specifying code
   NZ.BFZ.10.HH?, where the ? character is a single-valued wildcard.

-  **NET (network)**: Two character network code **NN**; e.g. NZ
-  **STA (station)**: Multi (3 or 4) character station code **SSS**:
   e.g. BFZ
-  **LOC (location)**: Two digit location code, for receiver sites that
   have multiple instruments (e.g. borehole and surface broadbands).
   e.g. 10

..

   **> unknown location code**: In my experience, location naming is not
   standard across data services. It is usually acceptable to provide
   wildcards to avoid having to determine the specific location code,
   e.g. NZ.BFZ.??.HH?

-  **CHA (channel)**: Three character channel code **BIO** (where
   **B**\ =Band code, **I**\ =Instrument code, **O**\ = Orientation
   code), e.g. HHZ

..

   **> seed convention**: When saving new data, band code is enforced
   via the data sampling rate, following the format of the `SEED
   convention <https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/>`__.

-  **COMP (component)**: One character channel component, a.k.a the
   orientation code in **CHANNEL**

..

   **> expected components**: In seismology **COMPONENT** is usually one
   of the following N: North, E: East, Z: Vertical, R: Radial, T:
   Transverse

--------------

Event ID
--------

The natural separation point for data objects is per event. Said
differently, processing in Pyatoa happens independently, one event at a
time. Unique event ID’s are **critical** for book keeping functionality
within Pyatoa. Since Pyatoa was developed in New Zealand, the standards
reflects the GeoNet Event ID formatting.

The following **standard** is enforced throughout Pyatoa:

-  Event IDs are generated from the resource_id of the ObsPy Event
   objects that are used to define event parameters.
-  Event IDs must be strings and all event ids within an inversion must
   be unique.
-  Each ASDFDataSet must correspond to a single event, and be saved
   using the associated event id.

Currently Pyatoa has been tested with the following source locations:

-  `GeoNet <https://quakesearch.geonet.org.nz/>`__
-  IRIS
-  .ndk files from SPUD, GCMT
-  USGS Events, and associated ANSS ComCat moment tensors

..

   **Note**: If unique event identifiers are required, you can monkey
   patch the function :func:``~pyatoa.utils.form.format_event_name``

--------------

Iteration
---------

Iteration number is used to keep track of the number of full iterations
(one iteration more or less contains a forward simulation, adjoint
simulation, and line search). The :class:``~pyatoa.core.config.Config``
class keeps track of the iteration as an internal attribute, which is
used for saving waveforms, auxiliary data, figures, etc.

The following **standard** is enforced throughout Pyatoa:

-  Iteration must be an integer that starts at 1. Iterations can not be
   negative.
-  The ``iteration tag`` is a string formatted version of the iteration
   used for saving data, and is formatted e.g. \ **i01** for iteration
   1.

--------------

Step Count
----------

Within a single iteration, multiple fuction evaluations take place,
relating to the intial model assessment (step 0), and each subsequent
line search function evaluation, each which try to reduce the objective
function via some step length of model perturbations (step count > 0).
As with iterations, the Config class keeps track of step count
information for internal tagging.

The following **standard** is enforced throughout Pyatoa:

-  Step count must be an integer that starts from 0. Step counts can not
   be negative.
-  The ``step tag`` is a string formatted version of the iteration used
   for saving data, and is formatted e.g. \ **s00** for step count 0.

--------------

Observation Waveforms
---------------------

Observation data that will not be gathered from online web-services,
e.g. temporary network data that is locally stored, must be available
using a specific directory and file naming format that mimics data
center convention. This is reflected in the function
:func:``~pyatoa.core.gatherer.Fetcher.fetch_obs_by_dir``.

The following **standard** is enforced throughout Pyatoa:

-  Default Directory Template: path/to/observed/YYYY/NN/SSS/CCC/
-  Default File ID Template: NN.SSS.LL.CCC.YYYY.DDD

Where:

-  YYYY: The year with the century (e.g., 1987)
-  NN: The network code (e.g. NZ)
-  SSS: The station code (e.g. BFZ)
-  LL: The location code (e.g. 10)
-  CCC: The channel code (e.g. HHZ.D)
-  DDD: The julian day of the year (January 1 is 001)

Example directory for station NZ.BFZ, for the day 2018-02-18:
**path/to/observed/2018/NZ/BFZ/HHZ/NZ.BFZ.10.HHZ.D.2018.049**

--------------

Station Response
----------------

As with observation data, response files that are stored locally on disk
(a.k.a StationXML files, dataless files) must be saved to a specific
directory and with specific file naming. These will be searched for by
the function :func:``pyatoa.core.gatherer.Fetcher.fetch_resp_by_dir``

The following **standard** is enforced throughout Pyatoa:

-  Default Directory Template: path/to/responses/\ **SSS.NN**
-  Default File ID Template: **RESP.NN.SSS.LL.CCC**

Where:

-  NN: The network code (e.g. NZ)
-  SSS: The station code (e.g. BFZ)
-  LL: The location code (e.g. 10)
-  CCC: The channel code (e.g. HHZ.D)

Example directory for station NZ.BFZ:
**path/to/response/BFZ.NZ/RESP.NZ.BFZ.10.HHZ**

--------------

Synthetic Waveforms
-------------------

Currently Pyatoa is written to work with ASCII synthetic seismograms
outputted by SPECFEM3D. Synthetic waveforms therefore must follow the
given naming convention set by SPECFEM3D.

The following **standard** is enforced throughout Pyatoa:

-  Synthetic waveforms must saved in the form:
   \__NN.SSS.BIO.__sem\ **U**

Where:

-  **N** = Network
-  **S** = Station
-  **B** = Band code
-  **I** = Instrument code (Must always be **X** for synthetics)
-  **O** = Orientation code
-  **U** = Unit code

Unit code **U** is dictated by the chosen output units in the SPECFEM3D
Par_file, where: 

-  **d** = displacement 
-  **v** = velocity 
-  **a** = acceleration

Example for displacement synthetic waveforms for the vertical component
of New Zealand station BFZ: **NZ.BFZ.BXZ**
