Processing w/ MPI
=================

The following example code processes (i.e., preprocess, window, generate
adjoint source) a number of synthetic waveforms in parallel using
Pyatoa and MPI (with `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_).

The script generates figures, adjoint sources and a text file with misfit
values for each source-receiver pair.

.. warning::

    This example requires `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_
    which is **not** a dependency of Pyatoa. To install mpi4py to the Conda
    environment created from the `installation <index.html#installation>`_ section, 
    you can run:

    .. code:: bash

        $ conda activate pyatoa
        $ conda install mpi4py

Example Data
------------

This example requires data stored in the Pyatoa repository. From
a working directory, you can grab this data using the following commands:

.. code:: bash

    cd ${PATH_TO_WORKDIR}
    cp -r ${PATH_TO_PYATOA}/pyatoa/tests/test_data/process_data_w_mpi.tar.gz .
    tar xf process_data_w_mpi.tar.gz

This will create a directory called ``data`` which contains synthetic waveforms
generated from `SeisFlows Example 2 
<https://seisflows.readthedocs.io/en/devel/specfem2d_example.html#example-2-checkerboard-inversion-w-pyaflowa-l-bfgs>`_.

In Example 2, which uses SPECFEM2D as a numerical solver, synthetic data are 
generated using a homogeneous halfspace model, while observed data are 
generated using a checkerboard perturbation model.

Script
------

To run this example on 4 cores, copy the script below to the filename
``process_data_w_mpi.py`` and run the following (note: you will need to
download the example data beforehand).

.. code:: bash

    mpirun -n 4 process_data_w_mpi.py

.. note::

    This script is also in the repo here:
    ``pyatoa/pyatoa/scripts/process_data_w_mpi.py``


.. code:: python

    #!/usr/bin/env python3
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from mpi4py import MPI
    from obspy import UTCDateTime

    from pysep.utils.io import read_sem
    from pyatoa import Config, Manager
    from pyatoa import logger as ptlogger
    from pyflex import logger as pflogger

    # Avoid overwhelming parallel log statements
    for logger in [ptlogger, pflogger]:
        logger.setLevel("CRITICAL")

    def get_rank_start_stop(ntask, nproc):
        """
        Determine the share of work that each processor gets based on the total
        number of tasks `ntask`, and the total number of processors `nproc`

        https://stackoverflow.com/questions/15658145/how-to-share-work-roughly-\
                              evenly-between-processes-in-mpi-despite-the-array-size

        :type ntask: int
        :param ntask: total number of tasks to divvy up
        :type nproc: int
        :param nproc: number of processors to split `ntask` over
        :rtype: list of tuples
        :return: the start and stop `ntask` number for each processor `nproc`
        """
        count = ntask // nproc
        remainder = ntask % nproc
        indices = []
        for rank in range(0, nproc):
            if rank < remainder:
                # The first 'remainder' ranks get 'count + 1' tasks each
                start = rank * (count + 1)
                stop = start + count
            else:
                # The remaining 'size - remainder' ranks get 'count' task each
                start = rank * count + remainder
                stop = start + (count - 1)
            indices.append((start, stop))

        return indices


    if __name__ == "__main__":
        # Initialize MPI
        comm = MPI.COMM_WORLD

        # Set up data structure and configuration parameters in rank 0
        if comm.rank == 0:
            # Define paths to data and output
            data_path = "./data/{ev}/{choice}/{sta}.semd"
            adjsrc_path = "./data/{ev}/adj"
            fig_path = "./figures"
            results_fid = "./misfit_results.txt"

            # Provides origin time for synthetic data which has none
            dummy_time = UTCDateTime("2000-01-01")

            # Create 30 unique event and station pairs
            _event_names = ["001", "002", "003"]
            _station_names = [f"AA.S{i:0>6}.BXY" for i in range(10)]
            evsta_pairs = []
            for event_name in _event_names:
                for sta_name in _station_names:
                    evsta_pairs.append((event_name, sta_name))

            # Generate paths for output results
            if not os.path.exists(fig_path):
                os.mkdir(fig_path)
            for ev in _event_names:
                if not os.path.exists(adjsrc_path.format(ev=ev)):
                    os.mkdir(adjsrc_path.format(ev=ev))

            # Determine how to divvy up the event-station pairs among processors
            indices = get_rank_start_stop(ntask=len(evsta_pairs), nproc=comm.size)

            # Generate Config object that controls processing
            config = Config(min_period=10., max_period=100., component_list=["Y"],
                            pyflex_preset="default", adj_src_type="cc_traveltime",
                            st_obs_type="syn", st_syn_type="syn"
                            )
        else:
            evsta_pairs = None
            data_path = None
            fig_path = None
            adjsrc_path = None
            results_fid = None
            dummy_time = None
            indices = None
            config = None

        # Broadcast generated data and config to each rank
        evsta_pairs = comm.bcast(evsta_pairs, root=0)
        data_path = comm.bcast(data_path, root=0)
        fig_path = comm.bcast(fig_path, root=0)
        adjsrc_path = comm.bcast(adjsrc_path, root=0)
        indices = comm.bcast(indices, root=0)
        dummy_time = comm.bcast(dummy_time, root=0)
        config = comm.bcast(config, root=0)

        if comm.rank == 0:
            print(f"{len(evsta_pairs)} total tasks to be accomplished with "
                  f"{comm.size} processors")

        # Each rank will process a different part of the event-station list
        start, stop = indices[comm.rank]

        # Initiate empty array to store misfit and measurement windows
        sendbuf = np.empty([stop - start + 1, 3], dtype=float)

        # Main processing for each rank: read data, process, write adjoint sources
        for i, evsta_pair in enumerate(evsta_pairs[start: stop + 1]):
            ev, sta = evsta_pair

            # Read in synthetic example data
            st_obs = read_sem(data_path.format(ev=ev, choice="obs", sta=sta),
                              origintime=dummy_time)
            st_syn = read_sem(data_path.format(ev=ev, choice="syn", sta=sta),
                              origintime=dummy_time)

            # Standard Pyatoa processing workflow
            mgmt = Manager(config=config, st_obs=st_obs, st_syn=st_syn)
            mgmt.standardize()
            mgmt.preprocess()
            mgmt.window()
            mgmt.measure()

            # Generate plot and adjoint source
            mgmt.write_adjsrcs(path=adjsrc_path.format(ev=ev), write_blanks=True)
            mgmt.plot(choice="wav", save=f"{fig_path}/{ev}_{sta}.png", show=False)
            plt.close()

            # Save misfit results to send buffer, which will be broadcast to Rank0
            sendbuf[i] = np.array([start + i, mgmt.stats.misfit, mgmt.stats.nwin],
                                  dtype=float)

        # Empty receiving buffer to collect results from all other ranks
        if comm.rank == 0:
            recvbuf = np.empty([len(evsta_pairs), 3], dtype=float)
        else:
            recvbuf = None

        # Gather misfit results from all ranks on Rank 0
        comm.Gather(sendbuf, recvbuf, root=0)

        # Use main rank to write out misfit information and number of windows
        if comm.rank == 0:
            with open(results_fid, "w") as f:
                for result in recvbuf:
                    idx, misfit, nwin = result
                    ev, sta = evsta_pairs[int(idx)]
                    f.write(f"{ev} {sta} {misfit:.2f} {int(nwin)}\n")


Results
-------

After successfully running the script, the results of the processing workflow
will be stored in multiple directories/files.

- ``figures/``: Contains waveform figures showing observed and synthetic traces, 
  misfit windows and adjoint sources for each source receiver pair.
- ``data/*/adj/*.adj``: Adjoint source files that are formatted as two-column 
  ASCII files (time v. amplitude), which are ready to be used by SPECFEM.
- ``misfit_results.txt``: A text file with information about misfit and number 
  of measurement windows used for each source-receiver pair.





