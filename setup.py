import os
from setuptools import setup, find_packages


# Get the list of required dependencies from requirements.txt
path_pyatoa = os.path.dirname(os.path.realpath(__file__))
req_file = os.path.join(path_pyatoa, "requirements.txt")
if os.path.exists(req_file):
    with open(req_file, "r") as f:
        install_requires = list(f.read().splitlines())
else:
    install_requires = []

setup(name='pyatoa',
      version='0.0.1',
      description="Python's Adjoint Tomography Operations Assitance",
      url='http://github.com/bch0w/pyatoa',
      author='Bryant Chow',
      author_email='bryant.chow@vuw.ac.nz',
      license='GPL',
      packages=find_packages(),
      install_requires=[],
      zip_safe=False)
