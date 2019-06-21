"""
Make sure the Manager class works as advertised
"""
import os
import obspy
import pyasdf
import unittest

from pyatoa import Config, Manager


class TestManager(unittest.TestCase):
    def setUpClass(cls):
        """
        set up the class
        :return:
        """
        cls.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'test_data', '')
        cls.ds = pyasdf.ASDFDataSet(
            os.path.join(cls.test_data_path, 'test_dataset.h5')
        )
        cls.st_obs = obspy.read(os.path.join(cls.test_data_path, 'test_obs_data.ms'))
        cls.st_syn = obspy.read(
            os.path.join(cls.test_data_path, 'test_syn_m00_data.ms')
        )
        cls.event = obspy.read_events(
            os.path.join(cls.test_data_path, 'test_event.xml')
        )
        cls.config = Config(
            model_number="m00",
            event_id=cls.event.resource_id.resource_id.split('/')[1]
        )


    def test_manager_init(self):
        """
        Make sure setup plot works for multiple ranges of windows and twin ax
        turned on and off
        :return:
        """
        mgmt = Manager(config=self.config, empty=False)

    def test_get_data(self):
        """
        Check that the data gathering process works
        Requires a working internet connection for FDSN query
        :return:
        """
        comp = "Z"
        station_code = self.st_obs.select(component=comp)[0].get_id()
        station_code = station_code[:-1] + "?"

        mgmt = Manager(config=self.config, empty=False)
        mgmt.gather_data(station_code=station_code)

        # Make sure gatherer downloaded the right data
        self.assertEqual(type(mgmt.st_obs), obspy.core.stream.Stream)
        self.assertEqual(len(mgmt.st_obs), 3)
        self.assertEqual(type(mgmt.inv), obspy.core.inventory.Inventory)
        self.assertEqual(type(mgmt.event), obspy.core.event.Event)

        # Check waveforms gathered are the same by checking their max values
        self.assertEqual(mgmt.st_obs.select(component=comp).data.max(),
                         self.st_obs.select(component=comp).data.max())

        # Check inventory captured the correct network
        self.assertEqual(mgmt.inv[0].code, self.inv[0].code)

        # Check event has the right origin time
        self.assertEqual(mgmt.event.preferred_origin().time,
                         self.event.preferred_origin().time)

    def test_preprocess(self):
        """
        Test that the preprocess function works
        :return:
        """
        mgmt = Manager(config=self.config, empty=True, event=self.event,
                       st_obs=self.st_obs, st_syn=self.st_syn, inv=self.inv)
        mgmt.preprocess()

        # Check that the proper processing has occurred
        for st in mgmt.st_obs:
            for tr in st:
                self.assertTrue(hasattr(tr.stats, 'processing'))
                self.assertEqual(tr.stats.processing, 13)

        for st in mgmt.st_syn:
            for tr in st:
                self.assertTrue(hasattr(tr.stats, 'processing'))
                self.assertEqual(tr.stats.processing, 9)

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
        mgmt = Manager(config=self.config, empty=True, event=self.event,
                       st_obs=self.st_obs, st_syn=self.st_syn, inv=self.inv)
        mgmt.preprocess()

        # Make sure preprocess recognized different units and adjusted syn's
        self.assertEqual(len(mgmt.st_obs[0].data),
                         len(self.st_obs[0].data) + 2 * config.zero_pad)

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
        misfit_check = 0.07336431876981384
        mgmt.run_pyadjoint()
        self.assertEqual(len(self.adj_srcs), 2)
        self.assertAlmostEqual(mgmt.total_misfit, misfit_check)


    def test_fill_dataset(self):
        """
        Test pyatoa's ability to fill a PyASDF ASDFDataSet
        :return:
        """
        self.empty_ds = pyasdf.ASDFDataSet(
            os.path.join(cls.test_data_path, 'test_empty_dataset.h5')
        )




