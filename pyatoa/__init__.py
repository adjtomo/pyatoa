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

from pyatoa.core.config import Config # NOQA
from pyatoa.core.manager import Manager, ManagerError # NOQA
from pyatoa.core.executive import Executive # NOQA
from pyatoa.core.gatherer import Gatherer  # NOQA
from pyatoa.core.inspector import Inspector  # NOQA
from pyatoa.utils.read import read_stations # NOQA
from pyatoa.utils.write import write_sem # NOQA

