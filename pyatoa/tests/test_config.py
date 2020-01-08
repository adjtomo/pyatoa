"""
Test the Config class, initialization, checks, reading and writing
"""
import os
import pyasdf
import unittest
from pyatoa import Config


class TestConfig(unittest.TestCase):
    def setUp(self):
        """
        Set up the class
        """
        self.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'data', '')

        self.empty_ds_fid = os.path.join(
            self.test_data_path, 'test_empty_ds.h5'
        )

    def test_print(self):
        """
        Testing the string representation
        """
        cfg = Config()
        print(cfg)

    def test_init(self):
        """
        Make sure config class can be initiated with various model numbers
        """
        model_number_check = "m00"
        for model_number in [0, '0', 'm00']:
            config = Config(model_number=model_number)
            self.assertTrue(config.model_number, model_number_check)

    def test_check(self):
        """
        Test the check function by providing incorrectly formated parameters
        """
        with self.assertRaises(AssertionError):
            # Check period range
            Config(min_period=100, max_period=10)
            # Check map corner dictionary
            Config(map_corners={'test': None})
            # Check unit output
            Config(unit_output='test')
            # Check synthetic output
            Config(cfgpaths={'test': 'test'})
            # Check window amplitude ratio
            Config(window_amplitude_ratio=-1)
            Config(window_amplitude_ratio=2)
            # External configs
            Config(pyflex_map='test')
            Config(adj_src_type='test')

    def test_read_write(self):
        """
        Ensure that the config can be written in each format
        """
        cfg = Config()
        for fid in ["test_config.yaml", "test_config.txt"]:
            cfg.write(fid)
            cfg.read(fid)
            if "yaml" in fid:
                # Make sure you can initiate using a yaml file
                Config(yaml_fid="test_config.yaml")
            os.remove(fid)

        # reading and writing from ASDF
        ds = pyasdf.ASDFDataSet(self.empty_ds_fid)
        cfg = Config()
        cfg.write(ds)
        cfg.read(ds, path="default")
        os.remove(self.empty_ds_fid)


if __name__ == "__main__":
    unittest.main()





