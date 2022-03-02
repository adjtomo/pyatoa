Inversion Prep
==============

Here we provide a general template for setting up a tomographic
inversion using ObsPy + Pyatoa. This will involve agathering an event
catalog with moment tensors, collecting observation waveforms and
response files, organizing data into the optimal directory structure,
and generating ASDFDataSets that can be used in a SeisFlows or
standalone inversion.

Alaska Event Catalog
--------------------

Alaska is a good region for an example problem, let’s work there. First
we’ll use ObsPy to gather our initial catalog of events from the past
decade in a box bounding Anchorage and Fairbanks.

.. code:: ipython3

    from obspy import UTCDateTime, Catalog
    from obspy.clients.fdsn import Client

.. code:: ipython3

    c = Client("USGS")
    cat = c.get_events(starttime=UTCDateTime("2010-01-01T00:00:00"), 
                       endtime=UTCDateTime("2020-01-01T00:00:00"), 
                       maxdepth=60.0,
                       minmagnitude=5.0,
                       maxmagnitude=6.0, 
                       minlatitude=59.75, 
                       maxlatitude=65.50, 
                       minlongitude=-154.5, 
                       maxlongitude=-143.789
                      )
    cat




.. parsed-literal::

    15 Event(s) in Catalog:
    2019-01-13T16:45:55.437000Z | +61.299, -150.065 | 5.0 ml | manual
    2019-01-06T03:45:34.525000Z | +65.407, -153.280 | 5.1 ml | manual
    ...
    2011-06-16T19:06:05.214000Z | +60.765, -151.076 | 5.1 mw | manual
    2011-01-23T02:50:04.629000Z | +63.542, -150.865 | 5.2 mw | manual
    To see all events call 'print(CatalogObject.__str__(print_all=True))'



Lets have a look at the Event IDs of our events. If we knew these
apriori, we could have gathered our catalog based on event ids

.. code:: ipython3

    from pyatoa.utils.form import format_event_name
    
    event_ids = []
    for event in cat:
        event_ids.append(format_event_name(event))
    event_ids




.. parsed-literal::

    ['ak019lrs7iu',
     'ak0199za3yf',
     'ak0191pccr7',
     'ak018fe5jk85',
     'ak20421672',
     'ak018fcpk9xi',
     'us1000hyge',
     'ak018fcntv5m',
     'ak018dsf3btv',
     'ak017f7s3c06',
     'ak014dlss56k',
     'ak014b5xf1in',
     'ak013ae2ycca',
     'ak0117oi3hnt',
     'ak011122ukq6']



Getting moment tensors
----------------------

Great, we have an event catalog now, but to do waveform simulations we
need moment tensors.

Unfortunately it’s not straight forward to grab moment tensor
information directly from USGS as they do not directly provide XML
files. It would be possible to manually generate moment tensor objects
from each `individual event
pages <https://earthquake.usgs.gov/earthquakes/eventpage/ak019lrs7iu/moment-tensor>`__,
but that seems tedious for a tutorial.

Instead we’ll use `GCMT <https://www.globalcmt.org/CMTsearch.html>`__.
Pyatoa has a function to read .ndk files hosted online with GCMT,
finding events based on origintime and magnitude.

.. code:: ipython3

    from pyatoa.core.gatherer import get_gcmt_moment_tensor
    
    events = []
    for event in cat:
        origintime = event.preferred_origin().time
        magnitude = event.preferred_magnitude().mag
        try:
            events.append(get_gcmt_moment_tensor(origintime, magnitude))
        except FileNotFoundError:
            print(f"No GCMT event found for: {format_event_name(event)}")
            continue
        
    gcmt_catalog = Catalog(events)
    print(f"\n{len(gcmt_catalog)}/{len(cat)} events with GCMT solutions found")


.. parsed-literal::

    No GCMT event found for: ak018fcpk9xi
    No GCMT event found for: us1000hyge
    No GCMT event found for: ak018fcntv5m
    No GCMT event found for: ak013ae2ycca
    
    11/15 events with GCMT solutions found


Great, 11 out of 15 isn’t bad, we’ll go ahead with and use the GCMT
catalog that we just collected. However if we wanted to retain the
(probably more accurate) origin information from the USGS catalog, we
would need to move the moment tensor objects from the GCMT catalog over
to the USGS catalog, an exercise left for the reader…

Gathering Observation Data
--------------------------

Now we need seismic waveform data for all the events in our catalog. We
can use the multithreaded data gathering functioality of Pyatoa’s
Gatherer class. First we need to determine the available broadband
stations in the area, using ObsPy.

Some pieces of relevant information that help motivate our search: \*
The Alaska Earthquake Center (AEC) operates stations under the network
code “AK”. \* The SEED standard seismometer instrument code is “H” \*
The SEED standard for broadband instruments is “B” or “H”

.. code:: ipython3

    c = Client("IRIS")
    inv = c.get_stations(network="AK", 
                         station="*", 
                         location="*",
                         channel="BH?",
                         starttime=UTCDateTime("2010-01-01T00:00:00"), 
                         endtime=UTCDateTime("2020-01-01T00:00:00"), 
                         minlatitude=59.75,                    
                         maxlatitude=65.50, 
                         minlongitude=-154.5, 
                         maxlongitude=-143.789,
                         level="channel"
                        )
    inv




.. parsed-literal::

    Inventory created at 2022-03-02T23:08:39.770000Z
    	Created by: IRIS WEB SERVICE: fdsnws-station | version: 1.1.48
    		    http://service.iris.edu/fdsnws/station/1/query?starttime=2010-01-01...
    	Sending institution: IRIS-DMC (IRIS-DMC)
    	Contains:
    		Networks (1):
    			AK
    		Stations (76):
    			AK.BMR (Bremner River, AK, USA)
    			AK.BPAW (Bear Paw Mountain, AK, USA)
    			AK.BRLK (Bradley Lake, AK, USA)
    			AK.BWN (Browne, AK, USA)
    			AK.CAPN (Captain Cook Nikiski, AK, USA)
    			AK.CAST (Castle Rocks, AK, USA)
    			AK.CCB (Clear Creek Butte, AK, USA)
    			AK.CHUM (Lake Minchumina, AK, USA)
    			AK.CUT (Chulitna, AK, USA)
    			AK.DDM (Donnely Dome, AK, USA)
    			AK.DHY (Denali Highway, AK, USA)
    			AK.DIV (Divide Microwave, AK, USA)
    			AK.DOT (Dot Lake, AK, USA)
    			AK.EYAK (Cordova Ski Area, AK, USA)
    			AK.FIB (Fire Island, AK, USA)
    			AK.FID (Fidalgo, AK, USA)
    			AK.FIRE (Fire Island, AK, USA)
    			AK.GHO (Gloryhole, AK, USA)
    			AK.GLB (Gilahina Butte, AK, USA)
    			AK.GLI (Glacier Island, AK, USA)
    			AK.GLM (Gilmore Dome, AK, USA)
    			AK.GOAT (Goat Mountain, AK, USA)
    			AK.HDA (Harding Lake, AK, USA) (2x)
    			AK.HIN (Hinchinbrook, AK, USA)
    			AK.HMT (Hamilton, AK, USA)
    			AK.I21K (Tanana, AK, USA)
    			AK.I23K (Minto, Yukon-Koyukuk, AK, USA)
    			AK.J20K (Nowitna River, AK, USA)
    			AK.J25K (Salcha River, AK, USA)
    			AK.K20K (Telida, AK, USA)
    			AK.K24K (Donnelly Dome, AK, USA)
    			AK.KAI (Kayak Island, AK, USA)
    			AK.KLU (Klutina Pass, AK, USA)
    			AK.KNK (Knik Glacier, AK, USA)
    			AK.KTH (Kantishna Hills, AK, USA)
    			AK.L20K (Farewell, AK, USA)
    			AK.L22K (Petersville, AK, USA)
    			AK.M19K (Big River Lodge, Big River, AK, USA)
    			AK.M20K (Styx River, AK, USA)
    			AK.MCK (McKinley Park, AK, USA)
    			AK.MDM (Murphy Dome, AK, USA)
    			AK.MLY (Manley Hot Springs, AK, USA)
    			AK.N19K (Bonanza Creek NPS repeater, AK, USA)
    			AK.NEA (Nenana, AK, USA)
    			AK.NEA2 (Nenana, AK, USA)
    			AK.NICH (Nichawak Mountain, AK, USA)
    			AK.NKA (Nikiski, AK, USA)
    			AK.O19K (Port Alsworth, AK, USA)
    			AK.O20K (Slope Mountain, AK, USA)
    			AK.P23K (Montague Island, AK, USA)
    			AK.PAX (Paxson, AK, USA)
    			AK.PPLA (Purkeypile, AK, USA)
    			AK.PWL (Port Wells, AK, USA)
    			AK.RAG (Ragged Mountain, AK, USA)
    			AK.RC01 (Rabbit Creek, AK, USA)
    			AK.RIDG (Independent Ridge, AK, USA)
    			AK.RND (Reindeer, AK, USA)
    			AK.SAW (Sawmill, AK, USA)
    			AK.SCM (Sheep Mountain, AK, USA)
    			AK.SCRK (Sand Creek, AK, USA)
    			AK.SGA (Sherman Glacier, AK, USA)
    			AK.SKN (Skwentna, AK, USA)
    			AK.SLK (Skilak Lake, AK, USA)
    			AK.SSN (Susitna, AK, USA)
    			AK.SWD (Seward, AK, USA)
    			AK.TRF (Thorofare Mountian, AK, USA) (2x)
    			AK.WAT1 (Susitna Watana 1, AK, USA)
    			AK.WAT2 (Susitna Watana 2, AK, USA)
    			AK.WAT3 (Susitna Watana 3, AK, USA)
    			AK.WAT4 (Susitna Watana 4, AK, USA)
    			AK.WAT5 (Susitna Watana 5, AK, USA)
    			AK.WAT6 (Susitna Watana 6, AK, USA)
    			AK.WAT7 (Susitna Watana 7, AK, USA)
    			AK.WRH (Wood River Hill, AK, USA)
    		Channels (543):
    			AK.BMR..BHZ (3x), AK.BMR..BHN (3x), AK.BMR..BHE (3x), 
    			AK.BPAW..BHZ (3x), AK.BPAW..BHN (3x), AK.BPAW..BHE (3x), 
    			AK.BRLK..BHZ (2x), AK.BRLK..BHN (2x), AK.BRLK..BHE (2x), 
    			AK.BWN..BHZ (3x), AK.BWN..BHN (3x), AK.BWN..BHE (3x), 
    			AK.CAPN..BHZ (2x), AK.CAPN..BHN (2x), AK.CAPN..BHE (2x), 
    			AK.CAST..BHZ (4x), AK.CAST..BHN (4x), AK.CAST..BHE (4x), 
    			AK.CCB..BHZ (3x), AK.CCB..BHN (3x), AK.CCB..BHE (3x), 
    			AK.CHUM..BHZ (2x), AK.CHUM..BHN (2x), AK.CHUM..BHE (2x), 
    			AK.CUT..BHZ (2x), AK.CUT..BHN (2x), AK.CUT..BHE (2x), 
    			AK.DDM..BHZ (2x), AK.DDM..BHN (2x), AK.DDM..BHE (2x), 
    			AK.DHY..BHZ (4x), AK.DHY..BHN (4x), AK.DHY..BHE (4x), 
    			AK.DIV..BHZ (4x), AK.DIV..BHN (4x), AK.DIV..BHE (4x), 
    			AK.DOT..BHZ (5x), AK.DOT..BHN (5x), AK.DOT..BHE (5x), 
    			AK.EYAK..BHZ (4x), AK.EYAK..BHN (4x), AK.EYAK..BHE (4x), 
    			AK.FIB..BHZ (2x), AK.FIB..BHN (2x), AK.FIB..BHE (2x), 
    			AK.FID..BHZ (4x), AK.FID..BHN (4x), AK.FID..BHE (4x), AK.FIRE..BHZ
    			AK.FIRE..BHN, AK.FIRE..BHE, AK.GHO..BHZ, AK.GHO..BHN, AK.GHO..BHE
    			AK.GLB..BHZ, AK.GLB..BHN, AK.GLB..BHE, AK.GLI..BHZ (3x), 
    			AK.GLI..BHN (3x), AK.GLI..BHE (3x), AK.GLM..BHZ, AK.GLM..BHN, 
    			AK.GLM..BHE, AK.GOAT..BHZ, AK.GOAT..BHN, AK.GOAT..BHE, 
    			AK.HDA..BHZ (2x), AK.HDA..BHN (2x), AK.HDA..BHE (2x), 
    			AK.HIN..BHZ (3x), AK.HIN..BHN (3x), AK.HIN..BHE (3x), 
    			AK.HMT..BHZ (2x), AK.HMT..BHN (2x), AK.HMT..BHE (2x), AK.I21K..BHZ
    			AK.I21K..BHN, AK.I21K..BHE, AK.I23K..BHZ, AK.I23K..BHN, 
    			AK.I23K..BHE, AK.J20K..BHZ, AK.J20K..BHN, AK.J20K..BHE, 
    			AK.J25K..BHZ, AK.J25K..BHN, AK.J25K..BHE, AK.K20K..BHZ, 
    			AK.K20K..BHN, AK.K20K..BHE, AK.K24K..BHZ, AK.K24K..BHN, 
    			AK.K24K..BHE, AK.KAI..BHZ (2x), AK.KAI..BHN (2x), AK.KAI..BHE (2x)
    			AK.KLU..BHZ (2x), AK.KLU..BHN (2x), AK.KLU..BHE (2x), 
    			AK.KNK..BHZ (2x), AK.KNK..BHN (2x), AK.KNK..BHE (2x), 
    			AK.KTH..BHZ (2x), AK.KTH..BHN (2x), AK.KTH..BHE (2x), AK.L20K..BHZ
    			AK.L20K..BHN, AK.L20K..BHE, AK.L22K..BHZ, AK.L22K..BHN, 
    			AK.L22K..BHE, AK.M19K..BHZ, AK.M19K..BHN, AK.M19K..BHE, 
    			AK.M20K..BHZ, AK.M20K..BHN, AK.M20K..BHE, AK.MCK..BHZ (3x), 
    			AK.MCK..BHN (3x), AK.MCK..BHE (3x), AK.MDM..BHZ (8x), 
    			AK.MDM..BHN (8x), AK.MDM..BHE (8x), AK.MLY..BHZ (4x), 
    			AK.MLY..BHN (4x), AK.MLY..BHE (4x), AK.N19K..BHZ, AK.N19K..BHN, 
    			AK.N19K..BHE, AK.NEA..BHZ, AK.NEA..BHN, AK.NEA..BHE, AK.NEA2..BHZ, 
    			AK.NEA2..BHN, AK.NEA2..BHE, AK.NICH..BHZ (5x), AK.NICH..BHN (5x), 
    			AK.NICH..BHE (5x), AK.NKA..BHZ, AK.NKA..BHN, AK.NKA..BHE, 
    			AK.O19K..BHZ, AK.O19K..BHN, AK.O19K..BHE, AK.O20K..BHZ, 
    			AK.O20K..BHN, AK.O20K..BHE, AK.P23K..BHZ, AK.P23K..BHN, 
    			AK.P23K..BHE, AK.PAX..BHZ (4x), AK.PAX..BHN (4x), AK.PAX..BHE (4x)
    			AK.PPLA..BHZ (2x), AK.PPLA..BHN (2x), AK.PPLA..BHE (2x), 
    			AK.PWL..BHZ (3x), AK.PWL..BHN (3x), AK.PWL..BHE (3x), 
    			AK.RAG..BHZ (3x), AK.RAG..BHN (3x), AK.RAG..BHE (3x), 
    			AK.RC01..BHZ (2x), AK.RC01..BHN (2x), AK.RC01..BHE (2x), 
    			AK.RIDG..BHZ (3x), AK.RIDG..BHN (3x), AK.RIDG..BHE (3x), 
    			AK.RND..BHZ (3x), AK.RND..BHN (3x), AK.RND..BHE (3x), 
    			AK.SAW..BHZ (4x), AK.SAW..BHN (4x), AK.SAW..BHE (4x), 
    			AK.SCM..BHZ (2x), AK.SCM..BHN (2x), AK.SCM..BHE (2x), 
    			AK.SCRK..BHZ (2x), AK.SCRK..BHN (2x), AK.SCRK..BHE (2x), 
    			AK.SGA..BHZ (4x), AK.SGA..BHN (4x), AK.SGA..BHE (4x), 
    			AK.SKN..BHZ (4x), AK.SKN..BHN (4x), AK.SKN..BHE (4x), AK.SLK..BHZ, 
    			AK.SLK..BHN, AK.SLK..BHE, AK.SSN..BHZ (5x), AK.SSN..BHN (5x), 
    			AK.SSN..BHE (5x), AK.SWD..BHZ (3x), AK.SWD..BHN (3x), 
    			AK.SWD..BHE (3x), AK.TRF..BHZ (4x), AK.TRF..BHN (4x), 
    			AK.TRF..BHE (4x), AK.WAT1..BHZ (3x), AK.WAT1..BHN (3x), 
    			AK.WAT1..BHE (3x), AK.WAT2..BHZ (3x), AK.WAT2..BHN (3x), 
    			AK.WAT2..BHE (3x), AK.WAT3..BHZ (4x), AK.WAT3..BHN (4x), 
    			AK.WAT3..BHE (4x), AK.WAT4..BHZ (4x), AK.WAT4..BHN (4x), 
    			AK.WAT4..BHE (4x), AK.WAT5..BHZ (2x), AK.WAT5..BHN (2x), 
    			AK.WAT5..BHE (2x), AK.WAT6..BHZ (2x), AK.WAT6..BHN (2x), 
    			AK.WAT6..BHE (2x), AK.WAT7..BHZ (3x), AK.WAT7..BHN (3x), 
    			AK.WAT7..BHE (3x), AK.WRH..BHZ (2x), AK.WRH..BHN (2x), 
    			AK.WRH..BHE (2x)



.. code:: ipython3

    # We'll need to create a list of station ids for data gathering
    station_codes = []
    for net in inv:
        for sta in net:
            station_codes.append(f"{net.code}.{sta.code}.*.BH?")
            
    # Let's just take a look at the first 10 as an example
    station_codes[:10]




.. parsed-literal::

    ['AK.BMR.*.BH?',
     'AK.BPAW.*.BH?',
     'AK.BRLK.*.BH?',
     'AK.BWN.*.BH?',
     'AK.CAPN.*.BH?',
     'AK.CAST.*.BH?',
     'AK.CCB.*.BH?',
     'AK.CHUM.*.BH?',
     'AK.CUT.*.BH?',
     'AK.DDM.*.BH?']



If we look at the inventory we see that there are 76 available stations
in our domain, quite a lot! Lets see how many have waveform data for the
events in question. We will do this by creating an ASDFDataSet for a
single event, and trying to fill it with all available data.

.. code:: ipython3

    from pyasdf import ASDFDataSet
    from pyatoa import Gatherer, Config

.. code:: ipython3

    # Here we are just using the first event in our catalog
    event = gcmt_catalog[0]
    event_id = format_event_name(event)

.. code:: ipython3

    # The gatherer needs to know where to look (Client) and when to look (origintime)
    cfg = Config(client="IRIS")
    origintime = event.preferred_origin().time

.. code:: ipython3

    # Now we initate the Gatherer and use its multithreading capabilities to gather waveform and metadata
    # Here the 'return_count' argument means we only want to save stations that return data including 
    # metadata (1) + 3 waveforms (3) = 4 
    
    # First make sure were writing to an empty dataset
    ds_fid = f"../tests/test_data/docs_data/{event_id}.h5"
    if os.path.exists(ds_fid):
        os.remove(ds_fid)
        
    with ASDFDataSet(ds_fid) as ds:
        ds.add_quakeml(event)
        gthr = Gatherer(config=cfg, ds=ds, origintime=origintime)
        gthr.gather_obs_multithread(codes=station_codes, return_count=4, print_exception=True)


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    /tmp/ipykernel_77870/1443690712.py in <module>
          5 # First make sure were writing to an empty dataset
          6 ds_fid = f"../tests/test_data/docs_data/{event_id}.h5"
    ----> 7 if os.path.exists(ds_fid):
          8     os.remove(ds_fid)
          9 


    NameError: name 'os' is not defined


.. code:: ipython3

    with ASDFDataSet(f"../tests/test_data/docs_data/{event_id}.h5") as ds:
        print(ds.waveforms.list())
        print(f"\n{len(ds.waveforms.list())} stations collected")


.. parsed-literal::

    ['AK.BMR', 'AK.BPAW', 'AK.BRLK', 'AK.BWN', 'AK.CAPN', 'AK.CAST', 'AK.CCB', 'AK.CHUM', 'AK.CUT', 'AK.DIV', 'AK.DOT', 'AK.EYAK', 'AK.FIRE', 'AK.GHO', 'AK.GLB', 'AK.GOAT', 'AK.HDA', 'AK.HIN', 'AK.HMT', 'AK.KAI', 'AK.KLU', 'AK.KNK', 'AK.KTH', 'AK.MCK', 'AK.NEA2', 'AK.NICH', 'AK.PAX', 'AK.PPLA', 'AK.PWL', 'AK.RAG', 'AK.RC01', 'AK.RIDG', 'AK.RND', 'AK.SAW', 'AK.SCM', 'AK.SCRK', 'AK.SKN', 'AK.SLK', 'AK.SWD', 'AK.TRF', 'AK.WRH']
    
    41 stations collected


Great! Looks like we’ve got data for 41 stations just for this one
event. Some stations did not return any data, as expected, but many of
them returned a StationXML plus three component waveforms (as explained
by data_count == 4).

--------------

Next Steps
~~~~~~~~~~

Now you can repeat the above data gathering steps for the remainder of
the events in your catalog. Each event should get it’s own ASDFDataSet
to keep data organized nicely. Take a look at the Storage tutorial to
get an idea of how to navigate and manipulate the ASDFDataSets. Also
have a look at the Pyaflowa tutorial in order to figure out how to
process the data you’ve just collected, either in a standalone manner
using Pyatao + SPECFEM3D, or with an automated workflow tool like
SeisFlows.
