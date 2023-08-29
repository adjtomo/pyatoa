"""
Test the I/O and sanity checks of the Config object
"""
import os
import pytest
import yaml
from pyasdf import ASDFDataSet
from pyatoa.core.config import Config


def test_io_asdf(tmpdir):
    """
    Saving a Config object using an ASDFDataSet will flatten and distribute
    the dictionaries and Pyflex and Pyadjoint Config objects. Ensure that when
    this functionality is called, that it doesn't lose any information in the
    conversion process.
    """
    test_value = 9999.

    # Set some ridiculous values to check upon re-reading
    cfg = Config(pyflex_parameters={'c_0': test_value},
                 pyadjoint_parameters={'dt_sigma_min': test_value}
                 )
 
    ds = ASDFDataSet(os.path.join(tmpdir, "test_dataset.h5"))
    cfg.write(write_to=ds)
    cfg_check = Config(ds=ds, path="default")
    assert cfg_check.pyflex_config.c_0 == test_value
    assert cfg_check.pyadjoint_config.dt_sigma_min == test_value

    del ds


def test_io_yaml(tmpdir):
    """
    Ensure that reading and writing to a Yaml file retains Config parameters
    """
    test_path = "this/is/a/test/path"
    test_value = 9999.

    # Set some ridiculous values to check upon re-reading
    cfg = Config(pyflex_parameters={'c_0': test_value},
                 pyadjoint_parameters={'dt_sigma_min': test_value}
                 )
    cfg.write(os.path.join(tmpdir, "test_config.yaml"))

    # Ensure ridiculous values read back in
    cfg_check = Config(yaml_fid=os.path.join(tmpdir, "test_config.yaml"))
    assert cfg_check.pyflex_config.c_0 == test_value
    assert cfg_check.pyadjoint_config.dt_sigma_min == test_value


def test_incorrect_io_yaml(tmpdir):
    """
    Ensure that wrong parameters cannot be passed in through yaml files
    """
    fid = os.path.join(tmpdir, "test_config.yaml")

    cfg = Config()
    cfg.write(fid)

    # Quickly edit the yaml file with an unrecognized parameter
    attrs = yaml.load(open(fid), Loader=yaml.Loader)
    attrs["test_parameter"] = 9999.
    with open(fid, "w") as f:
        yaml.dump(attrs, f, default_flow_style=False, sort_keys=False)

    with pytest.raises(ValueError):
        Config(yaml_fid=fid)


def test_incorrect_parameter_check():
    """
    Check that incorrect values passed to Config will set off the correct errors 
    """
    cfg = Config()

    # Try to set iteration and step incorrectly
    with pytest.raises(AssertionError):
        Config(iteration=0)
        Config(step_count=-1)

    # Try to use an unacceptable period range
    with pytest.raises(AssertionError):
        setattr(cfg, "min_period", 100)
        setattr(cfg, "max_period", 10)
        cfg._check()

    # Try to set unaccetapable parameter intputs
    incorrect_data = {"unit_output": "DISPLACEMENT",
                      "synthetic_unit": "ACCELERATION",
                      }
    with pytest.raises(AssertionError):
        for key, value in incorrect_data.items():
            cfg = Config()
            setattr(cfg, key, value)
            cfg._check()
