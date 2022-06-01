import os
from setuptools import setup, find_packages


# Check if we're on the RTD, don't require dependencies if so because we aren't
# actually running the code
READ_THE_DOCS = os.environ.get("READTHEDOCS", None) == "True"


# RTD doesn't require dependencies
if READ_THE_DOCS:
    install_requires = []
else:
    install_requires=[
            # Install requires has been shifted to 'requirements.txt' to avoid
            # Using Pip for package installation. Conda install preferred
            "obspy>=1.2.2",
            "pyasdf>=0.7.2",
            "pandas>=1.1.0",
            "pypdf2>=1.26.0",
            "pyyaml>=5.4",
            "matplotlib>=3.0.3",  
            "cartopy>=0.20.2",
            "pillow>=8.4.0",
            "pyflex @ git+https://github.com/bch0w/pyflex",
            "pyadjoint @ git+https://github.com/bch0w/pyadjoint"
            # "pyflex @ git+https://github.com/adjtomo/pyflex",
            # "pyadjoint @ git+https://github.com/adjtomo/pyadjoint"
            ]

setup(name='pyatoa',
      version='0.1.0',
      description="Python's Adjoint Tomography Operations Assistant",
      url='http://github.com/bch0w/pyatoa',
      download_url='https://github.com/bch0w/pyatoa/archive/refs/tags/v0.1.0.tar.gz',
      author='Bryant Chow',
      author_email='bhchow@alaska.edu',
      license='GPL-3.0',
      python_requires=">=3.7",
      packages=find_packages(),
      install_requires=install_requires,
      zip_safe=False
      )
