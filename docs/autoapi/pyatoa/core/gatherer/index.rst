:py:mod:`pyatoa.core.gatherer`
==============================

.. py:module:: pyatoa.core.gatherer

.. autoapi-nested-parse::

   Mid and Low level data gathering classes to retrieve data from local filesystems
   either via an ASDFDataSet or through a pre-defined data directory structure.

   Gatherer directly called by the Manager class and shouldn't need to be called
   by the User unless for bespoke data gathering functionality.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   pyatoa.core.gatherer.Gatherer




.. py:exception:: GathererNoDataException

   Bases: :py:obj:`Exception`

   Custom exception to be thrown generally to Manager class in the case that
   gathering of any data fails.


.. py:class:: Gatherer(config, ds=None, origintime=None)

   A mid-level data gathering class used to get data internally and externally.
   All saving to ASDFDataSet taken care of by the Gatherer class.

   .. py:method:: gather_event(event_id=None, **kwargs)

      Gather an ObsPy Event object by searching disk
      Event info need only be retrieved once per Pyatoa workflow.

      :type event_id: str
      :param event_id: a unique event idenfitier to search and tag event info
      :rtype: obspy.core.event.Event
      :return: event retrieved either via internal or external methods
      :raises GathererNoDataException: if no event information is found.


   .. py:method:: gather_station(code, **kwargs)

      Gather StationXML information from disk

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.inventory.Inventory
      :return: inventory containing relevant network and stations


   .. py:method:: gather_observed(code, **kwargs)

      Gather observed waveforms from disk as ObsPy streams

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.stream.Stream
      :return: stream object containing relevant waveforms


   .. py:method:: gather_synthetic(code, **kwargs)

      Gather synthetic waveforms as ObsPy streams.

      Only possible to check ASDFDataSet and local filesystem, cannot gather
      synthetics from webservice.

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.stream.Stream
      :return: stream object containing relevant waveforms
      :raises GathererNoDataException: if no synthetic data is found


   .. py:method:: fetch_event_from_dataset()

      Return Event information from ASDFDataSet.

      .. note::
          Assumes that the ASDF Dataset will only contain one event, which is
          dictated by the structure of Pyatoa.

      :rtype event: obspy.core.event.Event
      :return event: event object
      :raises AttributeError: if no event attribute found in ASDFDataSet
      :raises IndexError: if event attribute found but no events


   .. py:method:: fetch_inv_from_dataset(code)

      Return StationXML from ASDFDataSet based on station code.

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :rtype: obspy.core.inventory.network.Network
      :return: network containing relevant station information
      :raises KeyError: if no matching StationXML found


   .. py:method:: fetch_waveform_from_dataset(code, tag)

      Return waveforms as Stream objects from ASDFDataSet.

      .. note:
          * Allows for wildcard selection of component (? or *)
          * Selects by component because synthetic channel naming may differ
          from observation channels.
          * Component is assumed to be the last index in the channel,
          following SEED convention.

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type tag: str
      :param tag: internal asdf tag labelling waveforms
      :rtype: obspy.core.stream.Stream
      :return: waveform contained in a stream, or None if no matching value


   .. py:method:: fetch_event_by_dir(event_id, prefix='', suffix='', format_=None, **kwargs)

      Fetch event information via directory structure on disk. Developed to
      parse CMTSOLUTION and QUAKEML files, but theoretically accepts any
      format that the ObsPy read_events() function will accept.

      Will search through all paths given until a matching source file found.

      .. note::
          This function will search for the following path
          /path/to/event_dir/{prefix}{event_id}{suffix}

          so, if e.g., searching for a CMTSOLUTION file in the current dir:
          ./CMTSOLUTION_{event_id}

          Wildcards are okay but the function will return the first match

      :type event_id: str
      :param event_id: Unique event identifier to search source file by.
          e.g., a New Zealand earthquake ID '2018p130600'. A prefix or suffix
          will be tacked onto this
      :rtype event: obspy.core.event.Event or None
      :return event: event object if found, else None.
      :type prefix: str
      :param prefix Prefix to prepend to event id for file name searching.
          Wildcards are okay.
      :type suffix: str
      :param suffix: Suffix to append to event id for file name searching.
          Wildcards are okay.
      :type format_: str or NoneType
      :param format_: Expected format of the file to read, e.g., 'QUAKEML',
          passed to ObsPy read_events. NoneType means read_events() will guess


   .. py:method:: fetch_inv_by_dir(code, resp_dir_template='{sta}.{net}', resp_fid_template='RESP.{net}.{sta}.{loc}.{cha}', **kwargs)

      Fetch station dataless via directory structure on disk.
      Will search through all paths given until StationXML found.

      .. note::
          Default path naming follows SEED convention, that is:
          path/to/dataless/{NET}.{STA}/RESP.{NET}.{STA}.{LOC}.{CHA}
          e.g. path/to/dataless/NZ.BFZ/RESP.NZ.BFZ.10.HHZ

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type resp_dir_template: str
      :param resp_dir_template: Directory structure template to search for
          response files. By default follows the SEED convention:
          'path/to/RESPONSE/{sta}.{net}/'
      :type resp_fid_template: str
      :param resp_fid_template: Response file naming template to search for
          station dataless. By default, follows the SEED convention:
          'RESP.{net}.{sta}.{loc}.{cha}'
      :rtype inv: obspy.core.inventory.Inventory or None
      :return inv: inventory containing relevant network and stations


   .. py:method:: fetch_observed_by_dir(code, obs_dir_template='{year}/{net}/{sta}/{cha}', obs_fid_template='{net}.{sta}.{loc}.{cha}.{year}.{jday:0>3}', **kwargs)

      Fetch observation waveforms via directory structure on disk.

      .. note::
          Default waveform directory structure assumed to follow SEED
          convention. That is:
          path/to/data/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/{FID}
          e.g. path/to/data/2017/NZ/OPRZ/HHZ.D/NZ.OPRZ.10.HHZ.D

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type obs_dir_template: str
      :param obs_dir_template: directory structure to search for observation
          data. Follows the SEED convention:
          'path/to/obs_data/{year}/{net}/{sta}/{cha}'
      :type obs_fid_template: str
      :param obs_fid_template: File naming template to search for observation
          data. Follows the SEED convention:
          '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
      :rtype stream: obspy.core.stream.Stream or None
      :return stream: stream object containing relevant waveforms, else None


   .. py:method:: fetch_synthetic_by_dir(code, syn_cfgpath='synthetics', syn_unit='?', syn_dir_template='', syn_fid_template='{net}.{sta}.*{cmp}.sem{dva}*', **kwargs)

      Fetch synthetic waveforms from Specfem3D via directory structure on
      disk, if necessary convert native ASCII format to Stream object.

      .. note::
          By default, synthetics will be searched for with the following path
          config.paths[syn_cfgpath]/syn_dir_template/syn_fid_template.format()

      :type code: str
      :param code: Station code following SEED naming convention.
          This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
          L=location, C=channel). Allows for wildcard naming. By default
          the pyatoa workflow wants three orthogonal components in the N/E/Z
          coordinate system. Example station code: NZ.OPRZ.10.HH?
      :type syn_cfgpath: str
      :param syn_cfgpath: Config.paths key to search for synthetic data.
          Defaults to 'synthetics', but for the may need to be set to
          'waveforms' in certain use-cases.
      :type syn_unit: str
      :param syn_unit: Optional argument to specify the letter used to
          identify the units of the synthetic data: For Specfem3D:
          ["d", "v", "a", "?"] 'd' for displacement, 'v' for velocity,
          'a' for acceleration. Wildcards okay. Defaults to '?'
      :type syn_dir_template: str
      :param syn_dir_template: Directory structure template to search for
          synthetic waveforms. Defaults to empty string
      :type syn_fid_template: str
      :param syn_fid_template: The naming template of synthetic waveforms
          defaults to "{net}.{sta}.*{cmp}.sem{syn_unit}"
      :rtype stream: obspy.core.stream.Stream or None
      :return stream: stream object containing relevant waveforms


   .. py:method:: save_waveforms_to_dataset(st, tag)

      Save waveformsm to the ASDFDataSet with a simple check for existence
      of dataset and save parameter. Passes if waveforms already exist while
      ignoring the PyASDF warning that gets thrown if waveforms exist.

      :type st: obspy.core.stream.Stream
      :param st: Stream object to be saved into the dataset
      :type tag: str
      :param tag: unique identifier to save the waveforms under



