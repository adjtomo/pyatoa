"""
Make sure the Manager class works as advertised
"""
import unittest

import os
import json
import glob
import shutil
import pyasdf

from pyatoa.core.config import Config
import pyatoa.plugins.seisflows.sfprocess as sfprocess
import pyatoa.plugins.seisflows.sfconfig as sfconfig


class TestSeisflowsPlugin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        set up the class
        :return:
        """
        cls.config = Config('m00', '2018p130600')

        cls.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'data', '')

        # Keep all the test files generated in a scratch directory
        cls.scratch_dir = os.path.join(cls.test_data_path, "scratch")
        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

        # Create the necessary directory structure
        cls.working_dir = os.path.join(cls.scratch_dir, "working_dir")
        cls.current_dir = os.path.join(cls.scratch_dir, "current_dir")
        cls.pyatoa_dir = os.path.join(cls.scratch_dir, "pyatoa_dir")
        cls.output_dir = os.path.join(cls.working_dir, "output")
        os.makedirs(cls.scratch_dir)
        os.makedirs(cls.working_dir)
        os.makedirs(cls.current_dir)
        os.makedirs(cls.output_dir)

        # Make a fake Json file in the output directory
        jsonfid = os.path.join(cls.output_dir, "seisflows_paths.json")
        jsondata = {"PYATOA_IO": cls.pyatoa_dir}
        with open(jsonfid, 'w') as f:
            json.dump(jsondata, f)

        # Copy the STATION file to the current dir for 'process'
        sem_data_dir = os.path.join(cls.current_dir, "DATA")
        os.makedirs(sem_data_dir)
        shutil.copy(src=os.path.join(cls.test_data_path, "test_STATIONS"),
                    dst=os.path.join(sem_data_dir, "STATIONS")
                    )

        # Copy the test dataset into the pyatoa data dir
        # for 'process' and 'finalize'
        data_dir = os.path.join(cls.pyatoa_dir, "data")
        os.makedirs(data_dir)
        shutil.copy(src=os.path.join(cls.test_data_path, "test_dataset.h5"),
                    dst=os.path.join(data_dir, "2018p130600.h5")
                    )

        # Create the config file
        sfconfig_fid = os.path.join(cls.pyatoa_dir, "sfconfig.json")
        sfconfig.sfconfig(fidout=sfconfig_fid)

        # Copy the synthetic data into the current dir
        syn_dir = os.path.join(cls.current_dir, 'traces', 'syn')
        os.makedirs(syn_dir)
        synthetics = glob.glob(
            os.path.join(cls.test_data_path,
                         'test_directories', 'synthetics', '*'))
        for syn in synthetics:
            shutil.copy(src=syn,
                        dst=os.path.join(syn_dir, os.path.basename(syn))
                        )

        # Initialize a parser object using the sfprocess function
        cls.parser = sfprocess.initialize_parser([
            '--mode', 'initialize',
            '--working_dir', cls.working_dir,
            '--event_id', cls.config.event_id,
            '--model_number', cls.config.model_number,
            '--step_count', 0,
            '--current_dir', cls.current_dir,
            '--suffix', 'try'
        ])

    @classmethod
    def tearDownClass(cls):
        """
        Remove all the scratch data that was created
        :return:
        """
        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

    def test_initialize_parser(self):
        """
        Check that argparser works. Also adds the parser attribute to
        the test class so it can be used by following functions
        :return:
        """
        self.assertTrue(self.parser.mode is "initialize")

    def test_assemble_paths(self):
        """
        Check that the
        :return:
        """
        paths, usrcfg = sfprocess.assemble_paths(self.parser)

        self.assertTrue("PYATOA_DATA" in paths.keys())
        self.assertTrue("adj_src_type" in usrcfg.keys())

    def test_initialize(self):
        """

        :return:
        """
        paths, usrcfg = sfprocess.assemble_paths(self.parser)
        sfprocess.initialize(self.parser)

        self.assertTrue(os.path.exists(paths['PYATOA_FIGURES']))
        self.assertFalse(os.path.exists(paths['MISFIT_FILE']))

    def test_finalize(self):
        """
        Check that the initialize mode works
        :return:
        """
        sfprocess.initialize(self.parser)
        sfprocess.finalize(self.parser)
        # Check that the vtk generation worked
        check_file = os.path.join(
            self.pyatoa_dir, 'figures', 'vtks', 'srcrcv_m00_1.vtk')
        self.assertTrue(os.path.exists(check_file))

        # Check that snapshot works
        self.assertTrue(os.path.exists(
            self.pyatoa_dir, 'data', 'snapshot', '2018p130600.h5')
            )

    def test_process(self):
        """
        Check that the main processing script works
        :return:
        """
        sfprocess.initialize(self.parser)

        # sfconfig will be read from pyatoa.plugins.seisflows.sfconfig
        sfprocess.process(self.parser)

        # Quick check to make sure all expected output files were created
        for check_path in [
            os.path.join(self.current_dir, 'DATA', 'STATIONS_ADJOINT'),
            os.path.join(self.current_dir, 'traces', 'adj', 'NZ.BFZ.HXE.adj'),
            os.path.join(self.pyatoa_dir, 'data', 'misfits', '2018p130600'),
            os.path.join(
                self.pyatoa_dir, 'data', 'm00', '2018p130600', 'wav_BFZ.png'),
            os.path.join(
                self.pyatoa_dir, 'data', 'm00', '2018p130600', 'map_BFZ.png'),
            os.path.join(self.pyatoa_dir, 'data', 'vtks', 'srcrcv_m00_1.vtk'),
            os.path.join(self.pyatoa_dir, 'data', 'vtks', 'events_m00_1.vtk'),
            os.path.join(self.pyatoa_dir, 'misfits.json')
        ]:
            self.assertTrue(check_path)


class TestSeisflowsPluginSynOnly(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        set up the class
        :return:
        """
        # Set the config to a synthetics only case
        cls.config = Config('m00', '2018p130600')

        cls.test_data_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'data', '')

        # Keep all the test files generated in a scratch directory
        cls.scratch_dir = os.path.join(cls.test_data_path, "scratch")
        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

        # Create the necessary directory structure
        cls.working_dir = os.path.join(cls.scratch_dir, "working_dir")
        cls.current_dir = os.path.join(cls.scratch_dir, "current_dir")
        cls.pyatoa_dir = os.path.join(cls.scratch_dir, "pyatoa_dir")
        cls.output_dir = os.path.join(cls.working_dir, "output")
        os.makedirs(cls.scratch_dir)
        os.makedirs(cls.working_dir)
        os.makedirs(cls.current_dir)
        os.makedirs(cls.output_dir)

        # Make a fake Json file in the output directory
        jsonfid = os.path.join(cls.output_dir, "seisflows_paths.json")
        jsondata = {"PYATOA_IO": cls.pyatoa_dir}
        with open(jsonfid, 'w') as f:
            json.dump(jsondata, f)

        # Copy the STATION file to the current dir for 'process'
        sem_data_dir = os.path.join(cls.current_dir, "DATA")
        os.makedirs(sem_data_dir)
        shutil.copy(src=os.path.join(cls.test_data_path, "test_STATIONS"),
                    dst=os.path.join(sem_data_dir, "STATIONS")
                    )

        # Create an empty test dataset to be filled
        data_dir = os.path.join(cls.pyatoa_dir, "data")
        os.makedirs(data_dir)
        cls.ds = pyasdf.ASDFDataSet(os.path.join(data_dir, "2018p130600.h5"))

        # Copy the synthetic data into the current dir and to the obs directory
        # so that Pyatoa searches for data there
        syn_dir = os.path.join(cls.current_dir, 'traces', 'syn')
        os.makedirs(syn_dir)

        obs_dir = os.path.join(cls.current_dir, 'traces', 'obs')
        os.makedirs(obs_dir)

        # Copy the synthetics into test obs and syn directories
        synthetics = glob.glob(
            os.path.join(cls.test_data_path,
                         'test_directories', 'synthetics', '*'))
        for syn in synthetics:
            shutil.copy(src=syn,
                        dst=os.path.join(syn_dir, os.path.basename(syn))
                        )
            shutil.copy(src=syn,
                        dst=os.path.join(obs_dir, os.path.basename(syn))
                        )

        # Create the config file and assign the proper paths
        sfconfig_fid = os.path.join(cls.pyatoa_dir, "sfconfig.json")
        sfconfig.sfconfig(fidout=sfconfig_fid, set_logging="debug",
                          synthetics_only=True
                          )
        # Initialize a parser object using the sfprocess function
        cls.parser = sfprocess.initialize_parser([
            '--mode', 'initialize',
            '--working_dir', cls.working_dir,
            '--event_id', cls.config.event_id,
            '--model_number', cls.config.model_number,
            '--step_count', 0,
            '--current_dir', cls.current_dir,
            '--suffix', 'try'
        ])

    @classmethod
    def tearDownClass(cls):
        """
        Remove all the scratch data that was created
        :return:
        """
        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

    def test_process(self):
        """
        Test that synthetics only case works
        :return:
        """
        sfprocess.initialize(self.parser)

        sfprocess.process(self.parser)



if __name__ == "__main__":
    unittest.main()








