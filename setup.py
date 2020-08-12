from setuptools import setup

setup(name='pyatoa',
      version='0.0.1',
      description="Python's Adjoint Tomography Operational Assitance",
      url='http://github.com/bch0w/pyatoa',
      author='Bryant Chow',
      author_email='bryant.chow@vuw.ac.nz',
      license='GPL',
      packages=['pyatoa', 'pyatoa.core', 'pyatoa.utils', 'pyatoa.visuals'],
      zip_safe=False)
