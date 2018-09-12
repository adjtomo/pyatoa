"""
make sure the fdsn getter works
"""
import os
import unittest

from obspy import read, read_events, read_inventory, UTCDateTime

from pyatoa.config import Config
from pyatoa.utils.data_gather.external_getter import Getter
from pyatoa.utils.data_gather.internal_fetcher import Fetcher

class TestDataGatherMethods(unittest.TestCase):
    def setUpClass(cls):
        """
        set up the classls
        :return:
        """
        cls.test_data_path  = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'test_data', '')

    def test_external_getter(self):
        """
        make sure external getter returns what you want
        :return:
        """
        config = Config(event_id='2014p240655')
        getter = Getter(config)
        event_got, inv_got, st_got = getter.get_all('NZ.BFZ.10.HH?')

        event_check = read_events(
            os.path.join(self.test_data_path),'test_event.xml')
        inv_check = read_inventory(
            os.path.join(self.test_data_path),'test_inv.xml')
        stream_check = read(
            os.path.join(self.test_path_data),'test_stream.mseed')

        # assert event_got == event_check
        # assert inv_got == inv_check
        # assert st_got == stream_check

    def test_internal_fetcher():
        from obspy import UTCDateTime
        from pyatoa.config import Config
        from pyatoa.utils.gathering.internal_fetcher import Fetcher
        config = Config(event_id='2018p130600',
            paths_to_waveforms=['/Users/chowbr/Documents/subduction/seismic',
            '/seis/prj/fwi/bchow/seismic','/geonet/seismic'],
            paths_to_responses=['/seis/prj/fwi/bchow/seed/RESPONSE'],
            path_to_faults=['']
            model_number=0
                        )
        origintime = UTCDateTime('2018-02-18T07:43:48.127644Z')
        fetcher = Fetcher(config, origintime=origintime)
        inv = fetcher.fetch_response('NZ.MRZ.*.HH?')
        inv = fetcher.fetch_response('XX.RD01.*.HH?')

        st = fetcher.fetch_by_event('NZ.MRZ.*.HH?')
        st = fetcher.fetch_by_event('XX.RD01.*.HH?')
        st = fetcher.fetch_by_directory('XX.RD01.*.HH?')

    def test_data_gatherer():
        import logging
        from importlib import reload
        from pyatoa.core.config import Config
        from pyatoa.utils.gathering.data_gatherer import Gatherer
        logger = logging.getLogger("pyatoa")
        logger.setLevel(logging.DEBUG)
        config = Config(event_id='2013p613824',
            paths_to_waveforms=['/Users/chowbr/Documents/subduction/seismic',
            '/seis/prj/fwi/bchow/seismic','/geonet/seismic'],
            paths_to_responses=['/seis/prj/fwi/bchow/seed/RESPONSE'],
            model_number=0
                        )
        gatherer = Gatherer(config)
        gatherer.gather_moment_tensor()
        gatherer.gather_all('NZ.MRZ.*.HH?')

    def test_workflow():
        import logging
        import pyasdf
        from pyatoa import Config, Processor
        logger = logging.getLogger("pyatoa")
        logger.setLevel(logging.DEBUG)
        ds = pyasdf.ASDFDataSet("testasdf.h5")
        config = Config(event_id='2018p130600',
            min_period=9,
            paths_to_waveforms=['/Users/chowbr/Documents/subduction/seismic',
            '/seis/prj/fwi/bchow/seismic','/geonet/seismic'],
            paths_to_responses=['/seis/prj/fwi/bchow/seed/RESPONSE'],
            model_number=0
                        )
        proc = Processor(config,ds=ds)
        sta = 'NZ.BKZ.10.HH?'
        proc.populate(sta)

        f = window_maker(st_obs=proc.crate.st_obs, st_syn=proc.crate.st_syn,
                         windows=proc.crate.windows, staltas=proc.crate.staltas,
                         adj_srcs=proc.crate.adj_srcs,
                         stalta_wl=proc.config.pyflex_config[0],
                         unit_output=proc.config.unit_output,
                         config=proc.config, figsize=(11.69, 8.27), dpi=100,
                         show=show
                         )
        m = generate_map(config=proc.config, event=proc.crate.event,
                         inv=proc.crate.inv, show_faults=True,
                         show=show, figsize=(11.69, 8.27), dpi=100
                         )

if __name__ == "__main__":
    test_external_getter()