"""
make sure the fdsn getter works
"""
import os
import unittest

import obspy
import pyasdf
from pyatoa.config import Config
from pyatoa.utils.gathering.external_getter import Getter
from pyatoa.utils.gathering.internal_fetcher import Fetcher


class TestDataGatherMethods(unittest.TestCase):
    def setUpClass(cls):
        """
        set up the classls
        :return:
        """
        cls.test_data_path  = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'test_data', '')
        cls.ds = pyasdf.ASDFDataSet(
            os.path.join(cls.test_data_path, 'test_dataset.h5')
        )
        cls.st_obs = obspy.read(
            os.path.join(cls.test_data_path, 'test_obs_data.ms')
        )
        cls.st_syn = obspy.read(
            os.path.join(cls.test_data_path, 'test_syn_m00_data.ms')
        )
        cls.catalog = obspy.read_events(
            os.path.join(cls.test_data_path, 'test_event.xml')
        )
        cls.event = cls.catalog[0]
        cls.config = Config(
            model_number="m00",
            event_id=cls.event.resource_id.resource_id.split('/')[1]
        )

    def test_external_getter(self):
        """
        make sure external getter returns what you want
        :return:
        """
        comp = "Z"
        getter = Getter(self.config)
        station_name = self.st_obs[0].get_id()
        station_name = station_name[:-1] + "?"

        # Check to see Getter can collect the relevant information
        event, inv, st = getter.get_all(station_name)

        # Check random bits of information to check if Got items are correct
        self.assertEqual(event.resource_id, self.event.resource_id)
        self.assertEqual(inv[0][0].latitude, self.inv[0][0].latitude)
        self.assertEqual(st.select(component=comp)[0].data.max(),
                         self.st_obs.select(component=comp).data.max()
                         )

    def test_internal_fetcher(self):
        """

        :return:
        """

    def test_data_gatherer():
        """

        :return:
        """


if __name__ == "__main__":
    test_external_getter()