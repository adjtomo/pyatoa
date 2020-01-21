#!/usr/bin/env python3
"""
Pyaflowa

The Seisflows plugin class that allows easy scripting of Pyatoa
functionality into a Seisflows workflow. Pre-written functionalities simplify
calls made in Seisflows to Pyatoa, to reduce clutter inside the workflow
"""
import os

import sys
import json
import time
import glob
import pyasdf
import pyatoa
import shutil
import logging
import argparse
import warnings
import traceback
import numpy as np

from pyatoa.utils.asdf.deletions import clean_ds
from pyatoa.utils.asdf.additions import write_stats_to_asdf
from pyatoa.utils.asdf.extractions import windows_from_ds
from pyatoa.utils.tools.io import create_stations_adjoint, write_misfit_json, \
    write_adj_src_to_ascii, write_misfit_stats, tile_combine_imgs


class Pyaflowa:
    """
    The plugin object that is created to exist within Seisflows, keep track of
    the Seisflows workflow, and create the necessary outputs when requested
    """
    def __init__(self, working_dir):
        """
        Pyaflowa only needs to know where the main Seisflows working directory
        is located. With this information it can create the internal directory
        structure that it uses to navigate around Seisflows

        :type working_dir: str
        :param working_dir: location of the main Seisflows working directory
        """
        

