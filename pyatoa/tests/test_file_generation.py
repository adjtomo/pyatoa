"""
Make sure the Manager class works as advertised
"""
import unittest

import os
import pyasdf
import shutil

from pyatoa.core.config import Config
import pyatoa.utils.tools.io as pyatoa_io


class TestFileGeneration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        set up the class
        :return:
        """
        cls.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'data', '')
        cls.ds_path = os.path.join(cls.test_data_path, 'test_dataset.h5')
        cls.config = Config('m00', '2018p130600')

        # Keep all the test files generated in a scratch directory
        cls.scratch_dir = os.path.join(cls.test_data_path, "scratch")
        if not os.path.exists(cls.scratch_dir):
            os.mkdir(cls.scratch_dir)
        cls.scratch_file = os.path.join(cls.scratch_dir, "scratchfile")

    @classmethod
    def tearDownClass(cls):
        """
        Testing finalizations
        :return:
        """
        # Remove the empty dataset
        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

    def setUp(self):
        """open dataset for each new test"""
        self.ds = pyasdf.ASDFDataSet(self.ds_path)

    def tearDown(self):
        """delete dataset after each test"""
        del self.ds
        if os.path.exists(self.scratch_file):
            os.remove(self.scratch_file)

    def test_write_misfit_json(self):
        """
        Test that the JSON misfit writer works and writes out the proper
        misfit. Only check one line
        :return:
        """
        io.write_misfit_json(
            ds=self.ds, model=self.config.model_number, step_count=0,
            fidout=self.scratch_file
        )
        check_line = 5
        check_statement = '"misfit":0.4827'
        with open(self.scratch_file, 'r') as f:
            lines = f.readlines()
        self.assertTrue(check_statement in lines[check_line])

    def test_write_misfit_stats(self):
        """
        Check that the simple misfit statistics writer works
        :return:
        """
        pyatoa_io.write_misfit_stats(
            ds=self.ds, model=self.config.model_number,
            fidout=self.scratch_file
        )

        check = '3.447875e-01'
        with open(self.scratch_file, 'r') as f:
            line = f.read()
        self.assertEqual(line.strip(), check)

    def test_create_srcrcv_vtk_single(self):
        """
        Source receiver vtk file generator
        :return:
        """
        pyatoa_io.create_srcrcv_vtk_single(
            ds=self.ds, model=self.config.model_number,
            pathout=self.scratch_dir, event_separate=True,
            utm_zone=60
        )

        # Check the srcrcv VTK file has been properly written
        filename_check = os.path.join(
            self.scratch_dir, "{}_{}.vtk".format(
                os.path.basename(self.ds.filename).split('.')[0],
                self.config.model_number)
        )
        self.assertTrue(os.path.exists(filename_check))
        check_line = 6
        check = "5.689377E+05"
        with open(filename_check, 'r') as f:
            lines = f.readlines()
        self.assertEqual(lines[check_line].split()[0], check)

        # Check the event VTK file has been properly written
        event_filename_check = os.path.join(
            self.scratch_dir, "{}_{}_event.vtk".format(
                os.path.basename(self.ds.filename).split('.')[0],
                self.config.model_number)
        )
        self.assertTrue(os.path.exists(event_filename_check))
        check_line = 5
        check = "4.401620E+05"
        with open(event_filename_check, 'r') as f:
            lines = f.readlines()
        self.assertEqual(lines[check_line].split()[0], check)

    def test_create_srcrcv_vtk_multiple(self):
        """
        Same as single test except with different input
        :return:
        """
        pyatoa_io.create_srcrcv_vtk_multiple(
            pathin=self.test_data_path, pathout=self.scratch_dir,
            model=self.config.model_number, utm_zone=60
        )
        # Check the srcrcv VTK file has been properly written
        filename_check = os.path.join(
            self.scratch_dir, "srcrcv_{}_1.vtk".format(
                self.config.model_number)
        )
        self.assertTrue(os.path.exists(filename_check))
        check_line = 6
        check = "5.689377E+05"
        with open(filename_check, 'r') as f:
            lines = f.readlines()
        self.assertEqual(lines[check_line].split()[0], check)

        # Check the event VTK file has been properly written
        event_filename_check = os.path.join(
            self.scratch_dir, "events_{}_1.vtk".format(
                self.config.model_number)
        )

        self.assertTrue(os.path.exists(event_filename_check))
        check_line = 5
        check = "4.401620E+05"
        with open(event_filename_check, 'r') as f:
            lines = f.readlines()
        self.assertEqual(lines[check_line].split()[0], check)

    def test_create_stations_adjoint(self):
        """
        Make sure the adjoint station writer works
        :return:
        """
        station_file = os.path.join(self.test_data_path, 'test_STATIONS')

        pyatoa_io.create_stations_adjoint(
            ds=self.ds, model=self.config.model_number,
            specfem_station_file=station_file, pathout=self.scratch_dir
        )
        check_filename = os.path.join(self.scratch_dir, 'STATIONS_ADJOINT')
        self.assertTrue(os.path.exists(check_filename))
        with open(check_filename, 'r') as f:
            lines = f.readlines()
        check_line = 0
        check = "HAZ"
        self.assertEqual(lines[check_line].split()[0], check)

    def test_write_adj_src_to_ascii(self):
        """

        :return:
        """
        import numpy as np

        pyatoa_io.write_adj_src_to_ascii(
            ds=self.ds, model=self.config.model_number, pathout=self.scratch_dir
        )
        # Check that adjoint source writing works
        check_adj = os.path.join(self.scratch_dir, "NZ.HAZ.HXE.adj")
        check_max = 21029.258174999999
        self.assertTrue(os.path.exists(check_adj))
        adj = np.loadtxt(check_adj)
        self.assertAlmostEqual(adj.max(), check_max)

        # Check that blank adjoint source generation works
        check_blank = os.path.join(self.scratch_dir, "NZ.WAZ.HXN.adj")
        self.assertTrue(os.path.exists(check_blank))
        blank = np.loadtxt(check_blank)
        self.assertFalse(np.any(blank[:, 1]))


if __name__ == "__main__":
    unittest.main()








