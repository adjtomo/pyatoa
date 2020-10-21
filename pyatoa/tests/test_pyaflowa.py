"""
Test the workflow management class Pyaflowa and its supporting classes
PathStructure and IO
"""
import pytest
from pyatoa.core.pyaflowa import IO, PathStructure, Pyaflowa

def test_path_structure_standalone(tmpdir):
    """
    Test that the path structure standalone creates the expected
    directory structure. Just make sure things work.
    """
    path_structure = PathStructure(workdir=tmpdir)
    path_structure.format(source_name="TEST_EVENT_0001")

def test_path_structure_seisflows(tmpdir):
    """
    Test that seisflows directory structure matches what's expected
    """








