"""
Test the pyatoa
"""
import os
import pytest
from pyasdf import ASDFDataSet
from pyatoa.core.config import Config


def test_read_config_from_yaml():
    """
    Test that reading from an external YAML file works
    """
    cfg = Config(yaml_fid="./test_data/seisflows/parameters.yaml")
    assert(cfg.client == "TEST_CLIENT")  # Check a random variable


def test_read_write_from_ASDFDataSet(tmpdir):
    """
    Initiate a dataset, write the Config in with unique parameters, then read
    it back from the same dataset. Also tests _check_io_format()
    """
    cfg = Config(client="TEST_CLIENT", min_period=1, max_period=2)
    with ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5")) as ds:
        cfg.write(write_to=ds)
        cfg_check = Config(ds=ds, path="default")
        assert cfg == cfg_check



def test_incorrect_parameter_check():
    """
    Check that incorrect values passed to Config will set off the correct errors 
    """
    cfg = Config()

    # unacceptable period range
    with pytest.raises(AssertionError):
        setattr(cfg, "min_period", 100)
        setattr(cfg, "max_period", 10)
        cfg._check()

    # unaccetapable parameter intputs
    incorrect_data = {"unit_output": "DISPLACEMENT",
                      "synthetic_unit": "ACCELERATION",
                      "cfgpaths": [],
                      "cfgpaths": {"wave"},
                      "win_amp_ratio": 1.5
                      }
    with pytest.raises(AssertionError):
        for key, value in incorrect_data:
            cfg = Config()
            setattr(cfg, key, value)
            cfg._check()

    # unused key word arguments should result in ValueError
    with pytest.raises(ValueError):
        cfg = Config(ununused_kwarg="I dont belong :(")




