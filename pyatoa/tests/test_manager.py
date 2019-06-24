"""
Make sure the Manager class works as advertised
"""
import unittest

import os
import obspy
import pyasdf
from pyatoa import Config, Manager
from pyatoa.utils.asdf.deletions import clean_ds


class TestManager(unittest.TestCase):
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

    def test_manager_init(self):
        """
        Make sure empty manager can be initiated
        :return:
        """
        with self.assertRaises(TypeError):
            mgmt = Manager()
        mgmt = Manager(config=self.config, empty=True)
        self.assertIsNone(mgmt.event)
        with self.assertWarns(UserWarning):
            mgmt.preprocess()
            mgmt.run_pyflex()
            mgmt.run_pyadjoint()
            mgmt.plot_wav()

    def test_preprocess(self):
        """
        Test that the preprocess function works
        :return:
        """
        mgmt = Manager(config=self.config, empty=True, event=self.event,
                       st_obs=self.st_obs, st_syn=self.st_syn, inv=self.inv)
        mgmt.preprocess()

        # Check that the proper processing has occurred
        for tr in mgmt.st_obs:
            self.assertTrue(hasattr(tr.stats, 'processing'))
            self.assertEqual(len(tr.stats.processing), 13)

        for tr in mgmt.st_syn:
            self.assertTrue(hasattr(tr.stats, 'processing'))
            self.assertEqual(len(tr.stats.processing), 9)

    def test_preprocess_unit(self):
        """

        :return:
        """
        config = self.config
        config.output_unit = "DISP"
        config.synthetic_unit = "VEL"

        comp = "Z"
        mgmt = Manager(config=self.config, empty=True, event=self.event,
                       st_obs=self.st_obs, st_syn=self.st_syn, inv=self.inv)
        mgmt.preprocess()

        # Make sure preprocess recognized different units and adjusted syn's
        self.assertNotEqual(mgmt.st_syn.select(component=comp)[0].max(),
                            self.st_syn.select(component=comp)[0].max())

    def test_preprocess_extras(self):
        """

        :return:
        """
        config = self.config

        config.zero_pad = 50
        expected_length = 42500
        mgmt = Manager(config=self.config, empty=True, event=self.event,
                       st_obs=self.st_obs, st_syn=self.st_syn, inv=self.inv)
        mgmt.preprocess()

        # Make sure preprocess recognized different units and adjusted syn's
        self.assertEqual(len(mgmt.st_syn[0].data), expected_length)

    def test_run_pyflex_pyadjoint(self):
        """

        :return:
        """
        assert(self.config.pyflex_config[0] == "default")
        assert(self.config.pyadjoint_config[0] == "multitaper_misfit")

        mgmt = Manager(config=self.config, empty=True, event=self.event,
                       st_obs=self.st_obs, st_syn=self.st_syn, inv=self.inv)
        mgmt.preprocess()

        # Simple check to make sure Pyflex found windows using default params
        mgmt.run_pyflex()
        self.assertEqual(len(mgmt.windows), 2)
        self.assertEqual(len(mgmt.staltas), 2)
        self.assertEqual(mgmt.num_windows, 2)

        # Check that Pyadjoint multitaper retrieves the correct misfit value
        misfit_check = 0.06395629407523007
        mgmt.run_pyadjoint()
        self.assertEqual(len(mgmt.adj_srcs), 2)
        self.assertAlmostEqual(mgmt.total_misfit, misfit_check)

    def test_fill_dataset(self):
        """
        Test pyatoa's ability to fill a PyASDF ASDFDataSet
        Will query FDSN
        :return:
        """
        with pyasdf.ASDFDataSet(self.empty_ds_fid) as test_ds:
            clean_ds(test_ds)
            mgmt = Manager(config=self.config, ds=test_ds, empty=True,
                           event=self.event, st_obs=self.st_obs,
                           st_syn=self.st_syn, inv=self.inv
                           )
            mgmt.populate_dataset()
            mgmt.preprocess()
            mgmt.run_pyflex()

            # Check that the dataset contains the correct auxiliary data
            check_model = self.config.model_number
            check_window = 'NZ_BFZ_E_0'
            check_param = ("left_index", 4111)

            self.assertTrue(hasattr(mgmt.ds.auxiliary_data, 'MisfitWindows'))
            self.assertTrue(hasattr(mgmt.ds.auxiliary_data.MisfitWindows,
                                    check_model))
            self.assertTrue(hasattr(
                mgmt.ds.auxiliary_data.MisfitWindows[check_model], check_window)
            )

            # Shorten call for cleanliness
            wd = mgmt.ds.auxiliary_data.MisfitWindows[check_model][check_window]
            self.assertTrue(hasattr(wd, 'parameters'))
            self.assertEqual(wd.parameters[check_param[0]], check_param[1])

            # Remove the contents of the test database for future testing
            clean_ds(test_ds)

        del test_ds
        os.remove(self.empty_ds_fid)


if __name__ == "__main__":
    unittest.main()








