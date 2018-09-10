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
        from pyatoa.core.config import Config
        from pyatoa.utils.gathering.data_gatherer import Gatherer
        logger = logging.getLogger("pyatoa")
        logger.setLevel(logging.DEBUG)
        config = Config(event_id='2018p130600',
            paths_to_waveforms=['/Users/chowbr/Documents/subduction/seismic',
            '/seis/prj/fwi/bchow/seismic','/geonet/seismic'],
            paths_to_responses=['/seis/prj/fwi/bchow/seed/RESPONSE'],
            model_number=0
                        )
        gatherer = Gatherer(config)
        gatherer.gather_all('NZ.MRZ.*.HH?')


if __name__ == "__main__":
    test_external_getter()