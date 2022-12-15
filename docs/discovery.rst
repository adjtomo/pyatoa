Data Discovery
==============

The Pyatoa Manager class has internal routines to automatically discover
observed and synthetic waveforms, and metadata, given an ASDFDataSet or a set
of search paths and rules.

This functionality replaces the direct passing of ObsPy objects to the Manager.
All search paths are provided to the ``Config.paths`` attribute.

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

- `{path}/{prefix}{event_id}{suffix}`
- `path/to/events/CMTSOLUTION_001test`
- `path/to/other_events/CMTSOLUTION_001test`

The optional `prefix` and `suffix` allow a User to add additional text around
the `event_id`. However the entire file name can be passed the `event_id`
directly.

The Manager will loop through available paths until a matching file is found.
If no file is found, the Manager will throw a ``GathererNoDataException``.

Station Metadata
~~~~~~~~~~~~~~~~

Station response information can be gathered via a local search path using
station and network code.

.. note::

    Response information is searched for with a specific path and filenaming
    scheme that needs to be adhered to.

.. code:: python

    from pyatoa import Config, Manager

    response_paths = ["path/to/responses", "path/to/other_responses"]

    cfg = Config(paths={"responses": response_paths})
    mgmt = Manager(config=cfg)
    mgmt.gather_stations(code="NZ.BFZ..HH?", loc="", cha="*")

The above code block will attempt to look for StationXML files under the
following path structure:

- `{path}/{net}.{sta}/RESP.{net}.{sta}.{loc}.{cha}`
- path/to/responses/NZ.BFZ/RESP.NZ.BFZ..*
- path/to/other_responses/NZ.BFZ/RESP.NZ.BFZ..*

The optional `loc` and `cha` parameters can be set as specific location and
channel codes, or as wildcards to make the search more general.

.. note::

    Users who use `PySEP <https://github.com/adjtomo/pysep>`__ to gather their
    data can automatically output station metadata in the required directory
    format.

Custom Path Structure
`````````````````````

Users who want to input their own custom path and filename structure for
gathering station metadata can do so. The path can make use of formatting codes
`net`, `sta`, `loc`, and `cha` (networkk, station, location, channel).

.. code:: python

    cfg = Config(paths={"responses": ["./"]})
    mgmt.gather_stations(code="NZ.BFZ..HH?", resp_dir_template="",
                         resp_fid_template="{net}_{sta}_RESPONSE")

The above code block will attempt to look for StationXML files under the
following path structure:

- `{path}//{net}_{sta}_RESPONSE`
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
the following path structure, taking into account the Config parameters
`start_pad` and `end_pad`, which define how much data to gather before and
after the origin time.

In this case, since the example origin time is at midnight (00:00:00),
the Manager knows to look for data on both 2000-01-01 and 1999-12-31.

.. note::

    This default path structure and filename is meant to adhere to the
    typical data storage in data centers.

.. note::

    Users who use `PySEP <https://github.com/adjtomo/pysep>`__ to gather their
    data can automatically output observed waveforms in this format.

- `{path}/{year}/{network}/{station}/{channel}*/{net}.{sta}.{loc}.{cha}.{year}.{jday:0>3}
- `path/to/waveforms/NZ/BFZ/**/NZ.BFZ..*.2000.001`
- `path/to/waveforms/NZ/BFZ/**/NZ.BFZ..*.1999.365`
- `path/to/other_waveforms/NZ/BFZ/**/NZ.BFZ..*.2000.001`
- `path/to/other_waveforms/NZ/BFZ/**/NZ.BFZ..*.1999.365`

Custom Path Structure
`````````````````````

Users who want to input their own custom path and filename structure for
gathering station metadata can do so. The path can make use of formatting codes
`net`, `sta`, `loc`, and `cha`, `year` and `jday` (networkk, station, location,
channel, year, julian day).

.. code:: python

    cfg = Config(paths={"observed": "./"})
    mgmt = Manager(config=cfg, start_pad=50, end_pad=200)
    mgmt.gatherer.origintime = UTCDateTime("2000-01-01T00:00:00")  # example
    mgmt.gather_observed(code="NZ.BFZ..*", obs_dir_template="",
                         obs_fid_template="{net}_{sta}_{year}.{jday}.ms")

The above code block will attempt to look for observed waveforms under the
following path structure:

- `{path}/{net}_{sta}_{year}.{jday}.ms"
- `./NZ_BFZ_2000.001.ms`
- `./NZ_BFZ_1999.365.ms`


Synthetic Waveforms
~~~~~~~~~~~~~~~~~~~

Synthetic waveforms can be discovered via local path, corresponding to a
given event ID.

.. code:: python

    from pyatoa import Config, Manager

    synthetic_paths = ["path/to/synthetics", "path/to/other_synthetics"]

    cfg = Config(paths={"synthetics": synthetic_paths})
    mgmt = Manager(config=cfg)
    mgmt.gather_synthetic(code="NZ.BFZ..HH?", syn_unit="?")

The above code block will attempt to look for synthetic waveforms under the
following path structure:

- `{path}/{net}.{sta}.*{cmp}.sem{dva}*`
- `path/to/synthetics/NZ.BFZ.*?.sem?*`
- `path/to/other_synthetics/NZ.BFZ.*?.sem?*`

.. note::

    The default filenaming structure is meant to match synthetics output by
    SPECFEM2D/3D/3D_GLOBE.

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

- `{path}/{syn_dir_template}/{syn_fid_template}`
- `./synthetics/NZ_BFZ_HH?`


Combined Gathering Call
~~~~~~~~~~~~~~~~~~~~~~~

Users can chain together all of the above gathering into a single call by
setting paths all together and calling the ``Manager.gather()`` function.

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
``GathererNoDataException``.
