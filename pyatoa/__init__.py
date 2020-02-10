#!/usr/bin/env python
import logging

# set up logging (copied from pyflex)
logger = logging.getLogger("pyatoa")
logger.setLevel(logging.WARNING)
# Prevent propagating to higher loggers.
logger.propagate = 0
# Console log handler.
ch = logging.StreamHandler()
# Add formatter
FORMAT = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT)
ch.setFormatter(formatter)
logger.addHandler(ch)


from pyatoa.core.config import Config # NOQA
from pyatoa.core.manager import Manager # NOQA
from pyatoa.core.gatherer import Gatherer # NOQA
from pyatoa.plugins.seisflows.pyaflowa import Pyaflowa  # NOQA
from pyatoa.plugins.seisflows.inspector import Inspector  # NOQA
