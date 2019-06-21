"""
make sure the fdsn getter works
"""
import os
import pyasdf
import unittest
import matplotlib as mpl

from obspy import read, read_events, read_inventory, UTCDateTime

from pyatoa import Config, Manager
from pyatoa.visuals.waveforms import setup_plot, window_maker


class TestWaveformPlot(unittest.TestCase):
    def setUpClass(cls):
        """
        set up the class
        :return:
        """
        cls.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'test_data', '')
        cls.ds = pyasdf.ASDFDataSet(
            os.path.join(cls.test_data_path, 'test_dataset.h5'))
        cls.st_obs = read(os.path.join(cls.test_data_path, 'test_obs_data.ms'))
        cls.st_syn = read(
            os.path.join(cls.test_data_path, 'test_syn_m00_data.ms'))
        cls.event = read_events(
            os.path.join(cls.test_data_path, 'test_event.xml'))
        cls.config = Config(
            model_number="m00",
            event_id=cls.event.resource_id.resource_id.split('/')[1],
        )
        cls.mgmt = Manager(cls.config, empty=True)
        cls.mgmt.launch(set_event=cls.event)

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
            f = window_maker(
                st_obs=self.st_obs, st_syn=self.st_syn, config=self.config,
            )
        return f

    def test_window_maker_only_staltas(self):
        """
        Check that window maker works plotting waveforms with STA/LTA info only
        """
        f = window_maker(
            st_obs=self.st_obs, st_syn=self.st_syn, config=self.config)
        )







if __name__ == "__main__":
    test_external_getter()