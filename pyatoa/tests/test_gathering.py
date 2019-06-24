"""
Make sure all the data gathering methods, including internal fetching
and external getting, are working properly
"""
import os
import unittest

import obspy
import pyasdf
from pyatoa.core.config import Config
from pyatoa.core.manager import Manager
from pyatoa.core.gatherer import Gatherer

from pyatoa.utils.gathering.external_getter import Getter
from pyatoa.utils.gathering.internal_fetcher import Fetcher
from pyatoa.utils.gathering import grab_auxiliaries


class TestDataGatherMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set up class with test data to check against gathered data
        :return:
        """
        # Set up data to use
        cls.test_data_path  = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'data')
        cls.ds = pyasdf.ASDFDataSet(
            os.path.join(cls.test_data_path, 'test_dataset.h5')
        )
        cls.inv = obspy.read_inventory(
            os.path.join(cls.test_data_path, 'test_inv.xml')
        )

        cls.st_obs = obspy.read(
            os.path.join(cls.test_data_path, 'test_obs_data.ascii')
        )
        cls.st_syn = obspy.read(
            os.path.join(cls.test_data_path, 'test_syn_m00_data.ascii')
        )

        # For internal file searching
        cls.station_name = cls.st_obs[0].get_id()
        cls.station_code = cls.station_name[:-1] + '?'

        # Create event information
        cls.catalog = obspy.read_events(
            os.path.join(cls.test_data_path, 'test_event.xml')
        )
        cls.event = cls.catalog[0]
        cls.origintime = cls.event.preferred_origin().time

        cls.test_comp = "Z"

        cls.test_dir = os.path.join(cls.test_data_path, 'test_directories')
        cls.config = Config(
            model_number="m00",
            event_id=cls.event.resource_id.resource_id.split('/')[1],
            paths={"waveforms": os.path.join(cls.test_dir, "waveforms"),
                   "synthetics": os.path.join(cls.test_dir, "synthetics"),
                   "responses": os.path.join(cls.test_dir, "responses")
                   }
        )

    def test_external_getter_geonet(self):
        """
        make sure external getter returns what you want
        :return:
        """
        client = "GEONET"
        getter = Getter(config=self.config, client=client)

        # Check that the FDSN queries work
        event = getter.event_get()
        self.assertEqual(event.resource_id, self.event.resource_id)

        inv = getter.station_get(station_code=self.station_code)
        self.assertEqual(inv[0][0].latitude, self.inv[0][0].latitude)

        st = getter.waveform_get(station_code=self.station_code)
        self.assertEqual(
            st.select(component=self.test_comp)[0].data.max(),
            self.st_obs.select(component=self.test_comp)[0].data.max()
        )

    def test_internal_fetcher_by_dir(self):
        """
        Test fetching by internal directories
        :return:
        """
        fetcher = Fetcher(config=self.config, origintime=self.origintime)

        # Make sure wildcard naming works
        st_obs = fetcher.fetch_obs_by_dir(station_code=self.station_code)
        self.assertEqual(len(st_obs), 3)
        self.assertEqual(
            st_obs.select(component=self.test_comp)[0].data.max(),
            self.st_obs.select(component=self.test_comp)[0].data.max()
        )
        import ipdb;ipdb.set_trace()
        st_syn = fetcher.fetch_syn_by_dir(station_code=self.station_code)
        self.assertEqual(len(st_syn), 3)
        self.assertEqual(
            st_syn.select(component=self.test_comp)[0].data.max(),
            self.st_syn.select(component=self.test_comp)[0].data.max()
        )

        # Make sure specific station naming works
        st_obs_one = fetcher.fetch_obs_by_dir(station_code=self.station_name)
        self.assertEqual(len(st_obs_one), 1)

        st_syn_one = fetcher.fetch_syn_by_dir(station_code=self.station_name)
        self.assertEqual(len(st_syn_one), 1)

        # Make sure station fetching works
        inv = fetcher.fetch_resp_by_dir(station_code=self.station_code)
        self.assertEqual(len(inv.networks), 1)
        self.assertEqual(len(inv.networks[0].stations), 1)
        self.assertEqual(len(inv.networks[0].stations[0]), 18)

        self.assertEqual(inv[0][0].latitude, self.inv[0][0].latitude)

        # Make sure specific station naming works
        inv = fetcher.fetch_resp_by_dir(station_code=self.station_name)
        self.assertEqual(len(inv.networks), 1)
        self.assertEqual(len(inv.networks[0].stations), 1)
        self.assertEqual(len(inv.networks[0].stations[0]), 6)

    def test_internal_fetcher_by_asdf(self):
        """
        Test fetching by an ASDF dataset
        :return:
        """
        # Check that an empty fetcher throws the correct errors when
        # an empty dataset is given
        fetcher = Fetcher(config=self.config, origintime=self.origintime)


        fetcher = Fetcher(config=self.config, ds=self.ds,
                          origintime=self.origintime)

        event = fetcher.asdf_event_fetch()
        self.assertEqual(event.resource_id.resource_id.split('/')[1],
                         self.config.event_id)

        inv = fetcher.asdf_station_fetch(station_code=self.station_code)
        self.assertEqual(inv[0][0].latitude, self.inv[0][0].latitude)

        st_obs = fetcher.asdf_waveform_fetch(station_code=self.station_code,
                                             tag=self.config.observed_tag)
        self.assertEqual(
            st_obs.select(component=self.test_comp)[0].data.max(),
            self.st_obs.select(component=self.test_comp)[0].data.max()
        )
        st_syn = fetcher.asdf_waveform_fetch(station_code=self.station_code,
                                             tag=self.config.synthetic_tag)
        self.assertEqual(
            st_syn.select(component=self.test_comp)[0].data.max(),
            self.st_syn.select(component=self.test_comp)[0].data.max()
        )

    def test_internal_fetcher_main(self):
        """
        Make sure the internal fetcher makes the correct choices
        :return:
        """
        fetcher = Fetcher(config=self.config, ds=self.ds)
        st_obs = fetcher.obs_waveform_fetch(station_code=self.station_code)
        self.assertEqual(
            st_obs.select(component=self.test_comp)[0].data.max(),
            self.st_obs.select(component=self.test_comp)[0].data.max()
        )

        st_syn = fetcher.syn_waveform_fetch(station_code=self.station_code)
        self.assertEqual(
            st_syn.select(component=self.test_comp)[0].data.max(),
            self.st_syn.select(component=self.test_comp)[0].data.max()
        )

        inv = fetcher.station_fetch(station_code=self.station_code)
        self.assertEqual(inv[0][0].latitude, self.inv[0][0].latitude)

    def test_grab_auxiliaries(self):
        """
        Test the auxiliary data gatherer, which gathers non-essential
        information such as moment tensors
        :return:
        """
        expected_mxx = 1314.88
        expected_dip1 = 55.0
        expected_mw = 4.9
        expected_date = obspy.UTCDateTime('2018-02-18T07:43:00.000000Z')

        # Grab moment tensor from github repo query
        mom_ten = grab_auxiliaries.grab_geonet_moment_tensor(
            event_id=self.config.event_id)

        self.assertEqual(mom_ten['Mxx'], expected_mxx)
        self.assertEqual(mom_ten['dip1'], expected_dip1)
        self.assertEqual(mom_ten['Mw'], expected_mw)
        self.assertEqual(mom_ten['Date'], expected_date)

        # Grab moment tensor from GCMT catalog. The New Zealand event used
        # for testing is not available so use 2011 Tohoku as reference
        tohoku_origin = obspy.UTCDateTime("2011-03-11T5:47:32.8")
        tohoku_mw = 9.1
        tohoku_scalar_moment = 5.312E22

        event = grab_auxiliaries.grab_gcmt_moment_tensor(datetime=tohoku_origin,
                                                         magnitude=tohoku_mw)
        self.assertTrue(
            event.preferred_focal_mechanism().moment_tensor.scalar_moment,
            tohoku_scalar_moment
        )

    def test_data_gatherer(self):
        """
        Data gatherer is the upper level class that Manager interacts with,
        make sure its functionalities work. This is somewhat redundant, as
        Fetcher and Getter are tested first.

        TO DO: write some more test cases where data gatherer is forced to look
        in various places. Similar to what is tested in Fetcher, but want to
        make sure that gathering looks properly
        :return:
        """
        # Test the Gatherer empty
        expected_double_couple = 0.77
        gatherer = Gatherer(config=self.config)

        event = gatherer.gather_event(append_focal_mechanism=False)
        self.assertEqual(event.resource_id.resource_id.split('/')[1],
                         self.config.event_id)

        gatherer.append_focal_mechanism()
        self.assertEqual(gatherer.event.preferred_focal_mechanism(
            ).moment_tensor.double_couple,
                         expected_double_couple
                         )

        inv = gatherer.gather_station(station_code=self.station_code)
        self.assertEqual(inv[0][0].latitude, self.inv[0][0].latitude)

        st_obs = gatherer.gather_observed(station_code=self.station_code)
        self.assertEqual(
            st_obs.select(component=self.test_comp)[0].data.max(),
            self.st_obs.select(component=self.test_comp)[0].data.max()
        )

        st_syn = gatherer.gather_synthetic(station_code=self.station_code)
        self.assertEqual(
            st_syn.select(component=self.test_comp)[0].data.max(),
            self.st_syn.select(component=self.test_comp)[0].data.max()
        )

    def test_manager_gather_data(self):
        """
        Check that the data gathering process works
        Requires a working internet connection for FDSN query
        :return:
        """
        comp = "Z"
        station_code = self.st_obs.select(component=comp)[0].get_id()
        station_code = station_code[:-1] + "?"

        mgmt = Manager(config=self.config, event=self.event)
        mgmt.gather_data(station_code=station_code)

        # Make sure gatherer downloaded the right data
        self.assertIsInstance(mgmt.st_obs, obspy.core.stream.Stream)
        self.assertIs(mgmt.st_syn, None)
        self.assertEqual(len(mgmt.st_obs), 3)
        self.assertIsInstance(mgmt.inv, obspy.core.inventory.Inventory)
        self.assertIsInstance(mgmt.event, obspy.core.event.Event)

        # Check waveforms gathered are the same by checking their max values
        self.assertEqual(mgmt.st_obs.select(component=comp)[0].data.max(),
                         self.st_obs.select(component=comp)[0].data.max())

        # Check inventory captured the correct network
        self.assertEqual(mgmt.inv[0].code, self.inv[0].code)

        # Check event has the right origin time
        self.assertEqual(mgmt.event.preferred_origin().time,
                         self.event.preferred_origin().time)


if __name__ == "__main__":
    unittest.main()