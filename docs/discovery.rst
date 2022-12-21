Data Discovery
==============

The :class:`Manager <pyatoa.core.manager.Manager>` class has internal
routines in the :meth:`gather <pyatoa.core.manager.Manager.gather>` function
that automatically discover/gather waveforms and metadata, given a set
of search paths and rules. These routines are defined in the source code by the
:class:`Gatherer <pyatoa.core.gatherer.Gatherer>` class.

This functionality replaces the direct passing of ObsPy objects to the Manager
as described in the `Misfit Quantification <misfit.html>`__ documentation. It
is meant to facilitate automated data discovery within larger workflow tools
like `SeisFlows <https://github.com/adjtomo/seisflows>`__.

Data can also be discovered via ASDFDataSets, which is explained in the
`storage docs page <storage.html>`__.

.. note::

    *In a nutshell*: Users provide path and tag information (event id,
    station code) and the Gatherer class builds a set of wild card search paths
    in which to look for and read corresponding data/metadata.

Event Metadata
~~~~~~~~~~~~~~

Event metadata can be gathered locally via a search path, using an event ID.

.. code:: python

    from pyatoa import Config, Manager

    event_paths = ["path/to/events", "path/to/other_events"]

    cfg = Config(paths={"events": event_paths})
    mgmt = Manager(config=cfg)
    mgmt.gather_events(event_id="001", prefix="CMTSOLUTION_", suffix="test")

The above code block will attempt to look for events under a specific path
structure:

.. code:: bash

    {path}/{prefix}{event_id}{suffix}

Which in this example would be:

- `path/to/events/CMTSOLUTION_001test`
- `path/to/other_events/CMTSOLUTION_001test`

The optional ``prefix`` and ``suffix`` allow a User to wrap text around
the ``event_id``.

The Manager will loop through available paths until a matching file is found.
If no file is found, the Manager will throw a
:class:`pyatoa.core.gatherer.GathererNoDataException`.

Station Metadata
~~~~~~~~~~~~~~~~

Station response information can be gathered via a local search path using the
station and network codes:

.. code:: python

    from pyatoa import Config, Manager

    response_paths = ["path/to/responses", "path/to/other_responses"]

    cfg = Config(paths={"responses": response_paths})
    mgmt = Manager(config=cfg)
    mgmt.gather_stations(code="NZ.BFZ..HH?", loc="", cha="*")

The above code block will attempt to look for StationXML files under the
following path structure:

.. code:: bash

    {path}/{net}.{sta}/RESP.{net}.{sta}.{loc}.{cha}

Which in this example would be be:

- path/to/responses/NZ.BFZ/RESP.NZ.BFZ..*
- path/to/other_responses/NZ.BFZ/RESP.NZ.BFZ..*

The optional `loc` and `cha` parameters can be set as specific location and
channel codes, or as wildcards to make the search more general.

.. note::

    Users who use `PySEP <https://github.com/adjtomo/pysep>`__ to gather their
    data can automatically output station metadata in this required directory
    format.

Custom Path Structure
`````````````````````

Users who want to input their own custom path and filename structure for
gathering station metadata can do so. The path can make use of formatting codes
`net`, `sta`, `loc`, and `cha` (network, station, location, channel).

.. code:: python

    cfg = Config(paths={"responses": ["./"]})
    mgmt.gather_stations(code="NZ.BFZ..HH?", resp_dir_template="",
                         resp_fid_template="{net}_{sta}_RESPONSE")

The above code block will attempt to look for StationXML files under the
following path structure:

.. code:: bash

    {path}//{net}_{sta}_RESPONSE

Which in this example would be be:

- `.//NZ_BFZ_RESPONSE`


Observed Waveforms
~~~~~~~~~~~~~~~~~~

Observed waveforms can be searched for via local paths under specific path
and filenaming structure. Observed waveforms require an event origin time
to search for specific waveform start and end times.

.. code:: python

    from obspy import UTCDateTime
    from pyatoa import Config, Manager

    waveform_paths = ["path/to/waveforms", "path/to/other_waveforms"]

    cfg = Config(paths={"observed": waveform_paths})

    mgmt = Manager(config=cfg, start_pad=50, end_pad=200)
    mgmt.gatherer.origintime = UTCDateTime("2000-01-01T00:00:00")  # example
    mgmt.gather_observed(code="NZ.BFZ..*")

The above code block will attempt to look for observed waveform data under
the following path structure:

.. code:: bash

    {path}/{year}/{network}/{station}/{channel}*/{net}.{sta}.{loc}.{cha}.{year}.{jday:0>3}

Which in this example would be be:

- `path/to/waveforms/NZ/BFZ/**/NZ.BFZ..*.2000.001`
- `path/to/waveforms/NZ/BFZ/**/NZ.BFZ..*.1999.365`
- `path/to/other_waveforms/NZ/BFZ/**/NZ.BFZ..*.2000.001`
- `path/to/other_waveforms/NZ/BFZ/**/NZ.BFZ..*.1999.365`


The gathering routine uses the
:class:`Config <pyatoa.core.config.Config>` attributes
``start_pad`` and ``end_pad`` to define how the search period over time. In this
example, the example origin time is at midnight (00:00:00), so
the Manager knows to look for data on both 2000-01-01 and 1999-12-31.

.. note::

    Users who use `PySEP <https://github.com/adjtomo/pysep>`__ to gather their
    data can automatically output observed waveforms in this format.


Custom Path Structure
`````````````````````

Users who want to input their own custom path and filename structure for
gathering station metadata can do so. The path can make use of formatting codes
`net`, `sta`, `loc`, and `cha`, `year` and `jday` (network, station, location,
channel, year, julian day).

.. code:: python

    cfg = Config(paths={"observed": "./"})
    mgmt = Manager(config=cfg, start_pad=50, end_pad=200)
    mgmt.gatherer.origintime = UTCDateTime("2000-01-01T00:00:00")  # example
    mgmt.gather_observed(code="NZ.BFZ..*", obs_dir_template="",
                         obs_fid_template="{net}_{sta}_{year}.{jday}.ms")

The above code block will attempt to look for observed waveforms under the
following path structure:

.. code:: bash

    {path}/{net}_{sta}_{year}.{jday}.ms

Which in this example would be be:

- `./NZ_BFZ_2000.001.ms`
- `./NZ_BFZ_1999.365.ms`


Synthetic Waveforms
~~~~~~~~~~~~~~~~~~~

Synthetic waveforms can be discovered via local path given a specific event ID.
The default filenaming structure is meant to match synthetics output by
SPECFEM2D/3D/3D_GLOBE.

.. code:: python

    from pyatoa import Config, Manager

    synthetic_paths = ["path/to/synthetics", "path/to/other_synthetics"]

    cfg = Config(paths={"synthetics": synthetic_paths})
    mgmt = Manager(config=cfg)
    mgmt.gather_synthetic(code="NZ.BFZ..HH?", syn_unit="?")

The above code block will attempt to look for synthetic waveforms under the
following path structure:

.. code:: bash

    {path}/{net}.{sta}.*{cmp}.sem{dva}*

Which in this example would be be:

- `path/to/synthetics/NZ.BFZ.*?.sem?*`
- `path/to/other_synthetics/NZ.BFZ.*?.sem?*`


Custom Path Structure
`````````````````````

Users who want to input their own custom path and filename structure for
gathering synthetics can do so.

.. code:: python

    cfg = Config(paths={"synthetics": "./"})
    mgmt = Manager(config=cfg)
    mgmt.gather_synthetic(code="NZ.BFZ..HH?", syn_unit="?",
                          syn_dir_template="synthetics",
                          syn_fid_template="{net}_{sta}_{cha}")

The above code block will attempt to look for synthetic waveforms under the
following path structure:

.. code:: bash

    {path}/{syn_dir_template}/{syn_fid_template}

Which in this example would be be:

- `./synthetics/NZ_BFZ_HH?`


Combined Gathering Call
~~~~~~~~~~~~~~~~~~~~~~~

Users can chain together all of the above gathering into a single call by
setting paths all together and calling the
:meth:`gather <pyatoa.core.manager.Manager.gather>` function.

.. code:: python

        cfg = Config(paths={"event": [path_to_events],
                            "responses": [path_to_responses],
                            "observed": [path_to_observed],
                            "synthetics": [path_to_synthetics]
                            })
        mgmt.gather(event_id="CMTSOlUTION_001", code="NZ.BFZ..HH?",
                    choice=["event", "inv", "st_obs", "st_syn"])

This function will run all gathering functions one after another. If any part
of the gathering call fails, the Manager will throw a
:class:`pyatoa.core.gatherer.GathererNoDataException`.
