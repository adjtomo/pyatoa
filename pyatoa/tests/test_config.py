"""
Test the Config object from Pyatoa
"""
import pytest
from pyatoa.core.config import Config


def set_and_check(config, attribute, value):
    """
    Convenience function to set an attribute and run the config check
    """
    setattr(config, attribute, value)
    config._check()


def test_print():
    """
    Test the print statement of the config object
    """
    assert Config().__str__()


def test_read_yaml():
    """
    Test reading the Config object from a yaml file
    """
    for yaml_fid in ["./test_data/sfconfig.yaml",
                     "./test_data/seisflows/parameters.yaml"]:
        cfg = Config(yaml_fid)
        assert cfg


def test_check():
    """
    Check that incorrect values will create errors in the Config check
    """
    cfg = Config()

    with pytest.raises(AssertionError):
        set_and_check(cfg, "min_period", 100)
        set_and_check(cfg, "max_period", 10)
        set_and_check(cfg, "map_corners", {"test": 5})
        set_and_check(cfg, "unit_output", "test")
        set_and_check(cfg, "synthetic_unit", "test")
        set_and_check(cfg, "cfgpaths", [])
        set_and_check(cfg, "window_amplitude_ratio", 1E3)


