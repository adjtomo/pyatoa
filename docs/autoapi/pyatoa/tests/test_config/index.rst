:py:mod:`pyatoa.tests.test_config`
==================================

.. py:module:: pyatoa.tests.test_config

.. autoapi-nested-parse::

   Test the I/O and sanity checks of the Config object



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.tests.test_config.test_io_asdf
   pyatoa.tests.test_config.test_io_yaml
   pyatoa.tests.test_config.test_incorrect_io_yaml
   pyatoa.tests.test_config.test_incorrect_parameter_check



.. py:function:: test_io_asdf(tmpdir)

   Saving a Config object using an ASDFDataSet will flatten and distribute
   the dictionaries and Pyflex and Pyadjoint Config objects. Ensure that when
   this functionality is called, that it doesn't lose any information in the
   conversion process.


.. py:function:: test_io_yaml(tmpdir)

   Ensure that reading and writing to a Yaml file retains Config parameters


.. py:function:: test_incorrect_io_yaml(tmpdir)

   Ensure that wrong parameters cannot be passed in through yaml files


.. py:function:: test_incorrect_parameter_check()

   Check that incorrect values passed to Config will set off the correct errors


