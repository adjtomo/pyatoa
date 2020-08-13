from setuptools import setup

setup(name='pyatoa',
      version='0.0.1',
      description="Python's Adjoint Tomography Operational Assitance",
      url='http://github.com/bch0w/pyatoa',
      author='Bryant Chow',
      author_email='bryant.chow@vuw.ac.nz',
      license='GPL',
      packages=['pyatoa', 'pyatoa.core', 'pyatoa.utils', 'pyatoa.visuals',
                'pyatoa.plugins', 'pyatoa.utils.asdf'],
      # install_requires=[
      #     # 'obspy==1.2.2',
      #     # 'pyasdf==0.7.2',
      #     # 'pandas==1.1.0',
      #     # 'pyyaml==5.3.1',
      #     'repo @ http://github.com/bch0w/pyflex/tarball/return_reject_windows',
      #     'repo @ http://github.com/krischer/pyadjoint/tarball/master'
      #     ],
      # dependency_link = [
      #     'http://github.com/bch0w/pyflex/tarball/return_reject_windows',
      #     'http://github.com/krischer/pyadjoint/tarball/master'
      #     ],
      zip_safe=False)
