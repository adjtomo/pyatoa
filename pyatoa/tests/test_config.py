"""
Make sure the Manager class works as advertised
"""
import unittest

import os
import obspy
from pyatoa import Config


class TestConfig(unittest.TestCase):
    def setUp(self):
        """
        set up the class
        :return:
        """
        self.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'data', '')
        self.inv = obspy.read_inventory(
            os.path.join(self.test_data_path, 'test_inv.xml')
        )
        self.st_obs = obspy.read(os.path.join(
            self.test_data_path, 'test_obs_data.ms'))
        self.st_syn = obspy.read(
            os.path.join(self.test_data_path, 'test_syn_m00_data.ms')
        )
        self.catalog = obspy.read_events(
            os.path.join(self.test_data_path, 'test_event.xml')
        )
        self.event = self.catalog[0]
        self.config = Config(
            model_number="m00",
            event_id=self.event.resource_id.resource_id.split('/')[1]
        )
        self.empty_ds_fid = os.path.join(
            self.test_data_path, 'test_empty_ds.h5'
        )

    def test_init(self):
        """
        Make sure config class can be created
        :return:
        """
        model_number_check = "m00"
        for model_number in [0, '0', 'm00']:
            config = Config(model_number=model_number)
            self.assertTrue(config.model_number, model_number_check)

    def test_incorrect_parameters(self):
        """

        :return:
        """
        with self.assertRaises(AssertionError):
            # Check period range
            Config(min_period=100, max_period=10)
            # Check map corner dictionary
            Config(map_corners={'test': None})
            # Check unit output
            Config(unit_output='test')
            # Check synthetic output
            Config(synthetic_ouput='test')
            # Check pyflex output
            Config(pyflex_config='incorrect_config')









if __name__ == "__main__":
    unittest.main()








