:py:mod:`pyatoa.scripts.process_data_w_mpi`
===========================================

.. py:module:: pyatoa.scripts.process_data_w_mpi

.. autoapi-nested-parse::

   Processes (preprocess, window and generate adjoint source) a large amount of
   synthetic waveform data in parallel using Pyatoa and MPI (with mpi4py).

   Writes out a text file with misfit values for each source-receiver pair,
   saves figures of waveform misfit for each pair, and generates adjoint sources
   for each of the source receiver pairs.

   .. note::

       Requires `mpi4py` which is not a dependency of Pyatoa. To install into your
       conda environment you can run:

       $ conda install mpi4py

   .. rubric:: usage

       # To run this on 4 cores
       mpirun -n 4 process_data_w_mpi.py



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.scripts.process_data_w_mpi.get_rank_start_stop



Attributes
~~~~~~~~~~

.. autoapisummary::

   pyatoa.scripts.process_data_w_mpi.comm


.. py:function:: get_rank_start_stop(ntask, nproc)

   Determine the share of work that each processor gets based on the total
   number of tasks `ntask`, and the total number of processors `nproc`

   https://stackoverflow.com/questions/15658145/how-to-share-work-roughly-                          evenly-between-processes-in-mpi-despite-the-array-size

   :type ntask: int
   :param ntask: total number of tasks to divy up
   :type nproc: int
   :param nproc: number of processors to split `ntask` over
   :rtype: list of tuples
   :return: the start and stop `ntask` number for each processor `nproc`


.. py:data:: comm
   

   

