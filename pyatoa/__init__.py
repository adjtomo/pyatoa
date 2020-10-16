#!/usr/bin/env python
import logging

# Set up logging
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.WARNING)  # Default level
logger.propagate = 0  # Prevent propagating to higher loggers
ch = logging.StreamHandler()  # Console log handler
FORMAT = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT, datefmt="%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)

# Actually this is annoying, don't do that
# Redirect warning messages into the logger
# logging.captureWarnings(True)
# warnings_logger = logging.getLogger("py.warnings")
# warnings_logger.addHandler(logger)

from pyatoa.core.config import Config # NOQA
from pyatoa.core.manager import Manager, ManagerError # NOQA
from pyatoa.core.gatherer import (Gatherer, append_focal_mechanism,
                                  get_gcmt_moment_tensor)  # NOQA
from pyatoa.core.inspector import Inspector  # NOQA
from pyatoa.core.pyaflowa import Pyaflowa  # NOQA
from pyatoa.utils.read import read_sem, read_stations # NOQA
from pyatoa.utils.write import write_sem, write_stations # NOQA

# Dont include this so that Mayavi is not a default requirement
#from pyatoa.visuals.vtk_modeler import VTKModeler  # NOQA

