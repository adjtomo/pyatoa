"""
UNFINISHED 4-7-19
"""
import os
import pyasdf
import unittest
import matplotlib as mpl

import obspy
from pyatoa import Config, Manager
from pyatoa.utils.visuals.waveforms import setup_plot, window_maker


class TestWaveformPlot(unittest.TestCase):
    def setUpClass(cls):
        """
        set up the class
        :return:
        """
        cls.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'test_data', '')

        cls.st_obs = obspy.read(
            os.path.join(cls.test_data_path, 'test_obs_data.ascii')
        )
        cls.st_syn = obspy.read(
            os.path.join(cls.test_data_path, 'test_syn_m00_data.ascii')
        )
        cls.catalog = obspy.read_events(
            os.path.join(cls.test_data_path, 'test_cat.xml')
        )
        cls.event = cls.catalog[0]
        cls.config = Config(
            model_number="m00",
            event_id=cls.event.resource_id.resource_id.split('/')[1],
        )

        cls.mgmt = Manager(cls.config, st_obs=cls.st_obs, st_syn=cls.st_syn,
                           event=cls.event
                           )

    def test_setup_plot(self):
        """
        Make sure setup plot works for multiple ranges of windows and twin ax
        turned on and off
        :return:
        """
        # Test with twin ax
        twax_bool = True
        for number_of in range(0, 4):
            axes, twaxes = setup_plot(number_of, twax_bool)
            self.assertEqual(len(axes), number_of)
            self.assertEqual(len(twaxes), number_of)
            for ax, tw in zip(axes, twaxes):
                self.assertEqual(type(ax), mpl.axes._subplots.Subplot)
                self.assertEqual(type(tw), mpl.axes._subplots.Subplot)

        # Test with not twin ax
        twax_bool = False
        for number_of in range(0, 4):
            axes, twaxes = setup_plot(number_of, twax_bool)
            self.assertEqual(len(axes), number_of)
            self.assertEqual(len(twaxes), 0)
            for ax in axes:
                self.assertEqual(type(ax), mpl.axes._subplots.Subplot)

    def test_window_maker_barebones(self):
        """
        Check that window maker works just as a waveform plotter, without
        the need to plot misfit windows or adjoint sources
        :return:
        """
        # Test data has uneven lengths, so plotting will throw a value error
        with self.assertRaises(ValueError):
            f = window_maker(st_obs=self.mgmt.st_obs, st_syn=self.mgmt.st_syn,
                             config=self.config)
        # Preprocess data and plot
        self.mgmt.preprocess()
        f = window_maker(st_obs=self.mgmt.st_obs, st_syn=self.mgmt.st_syn,
                         config=self.config)

        # Run pyflex and plot staltas and windows
        self.mgmt.run_pyflex()
        f = window_maker(st_obs=self.mgmt.st_obs, st_syn=self.mgmt.st_syn,
                         config=self.config, staltas=self.mgmt.staltas,
                         windows=self.mgmt.windows)

        self.mgmt.run_pyadjoint()
        f = window_maker(st_obs=self.mgmt.st_obs, st_syn=self.mgmt.st_syn,
                         config=self.config, staltas=self.mgmt.staltas,
                         windows=self.mgmt.windows)

if __name__ == "__main__":
    test_external_getter()
