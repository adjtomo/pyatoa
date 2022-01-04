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
      version='0.1.0',
      description="Python's Adjoint Tomography Operations Assistant",
      url='http://github.com/bch0w/pyatoa',
      download_url='https://github.com/bch0w/pyatoa/archive/refs/tags/v0.1.0.tar.gz',
      author='Bryant Chow',
      author_email='bhchow@alaska.edu',
      license='GPL-3.0',
      packages=find_packages(),
      install_requires=[],
      zip_safe=False)
