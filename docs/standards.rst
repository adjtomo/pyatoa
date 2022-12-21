Naming Standards
=================

Pyatoa enforces a few standards for book keeping which are listed here.

Receiver codes
--------------

Receivers codes define station names and components for a given station.
Wildcards are typically acceptable when specifying station codes in Pyatoa.

-  **CODE**: Standard code are expected in format **NN.SSS.LL.CCC** (where
   **N**\ =Network, **S**\ =Station, **L**\ =Location, **C**\ =Channel);
   e.g. `NZ.BFZ.10.HHZ <https://www.geonet.org.nz/data/network/sensor/BFZ>`__
-  **NET (network)**: Two character network code **NN** (e.g. NZ)
-  **STA (station)**: Multi character station code **SSS** (e.g.BFZ)
-  **LOC (location)**: Two digit location code (e.g., 10).
-  **CHA (channel)**: Three character channel code **BIO** (where
   **B**\ =Band code, **I**\ =Instrument code, **O**\ = Orientation
   code); e.g.HHZ
-  **COMP (component)**: One character channel component, a.k.a the
   orientation code in **CHANNEL**

.. note::

    When saving new data, band code is enforced via the data sampling rate,
    following the format of the `SEED convention
    <https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/>`__.

Event ID
--------

Unique event ID’s are critical for book keeping. The following **standard** is
enforced throughout Pyatoa:

-  Event IDs are generated from: ``Obspy.Event.resource_id.id``
-  Event IDs must be strings
-  All Event IDs must be unique for a given inversion
-  ASDFDataSets must correspond to a single event, saved as: ``{event_id}.h5``

Iteration
---------

Iteration number is used to keep track of the number of full iterations in an
inversion.

-  Iteration must be an integer [1, inf).
-  ``iteration tag``: string formatted iteration. E.g., **i01** for iteration 1.


Step Count
----------

During a line search, multiple forward simulations may take place. Step Count
is used to track each forward simulation in a line search.

-  Step count must be an integer [0, inf).
-  ``step tag``: string formatted step count. E.g., **s00** for step count 0.

Evaluation
----------

An evaluation is a single step within a single inversion.

- Evaluation is a string comprised of ``iteration_tag`` and ``step_tag``
- The first evaluation of an inversion is: `i01s00`
- Whenever iterations are incremented, step counts reset to 0


Observation Waveforms
---------------------

Related to the `Data Discovery page <discovery.html>`__, Observation waveforms
stored locally (on disk) should be saved in a specific format:

-  Default Directory Template: ``path/to/observed/YYYY/NN/SSS/CCC/``
-  Default File ID Template: ``NN.SSS.LL.CCC.YYYY.DDD``

Where:

-  YYYY: Year with century (e.g., 1987)
-  NN: Network code (e.g.NZ)
-  SSS: Station code (e.g.BFZ)
-  LL: Location code (e.g.10)
-  CCC: Channel code (e.g.HHZ.D)
-  DDD: Julian day of the year (January 1 is 001)

Example directory for station NZ.BFZ, for the day 2018-02-18:
**path/to/observed/2018/NZ/BFZ/HHZ/NZ.BFZ.10.HHZ.D.2018.049**


Station Response
----------------

Related to the `Data Discovery page <discovery.html>`__, station metadata or
response information stored locally (on disk) should be saved in a specific
format:

-  Default Directory Template: path/to/responses/**SSS.NN**
-  Default File ID Template: **RESP.NN.SSS.LL.CCC**

Where:

-  NN: Network code (e.g. NZ)
-  SSS: Station code (e.g. BFZ)
-  LL: Location code (e.g. 10)
-  CCC: Channel code (e.g. HHZ.D)

Example directory for station NZ.BFZ:
**path/to/response/BFZ.NZ/RESP.NZ.BFZ.10.HHZ**

Synthetic Waveforms
-------------------
Synthetic waveforms therefore must follow the naming convention dictated by
SPECFEM.

-   NN.SSS.BIO.sem**U**

Where:

-  **N** = Network
-  **S** = Station
-  **B** = Band code
-  **I** = Instrument code (**Must always be **X** for synthetics**)
-  **O** = Orientation code
-  **U** = Unit code

Unit code **U** is dictated by the chosen output units in SPECFEM, e.g.:

-  **d** = displacement 
-  **v** = velocity 
-  **a** = acceleration

Example for displacement synthetic waveforms for the vertical component
of New Zealand station BFZ: **NZ.BFZ.BXZ.sem\***
