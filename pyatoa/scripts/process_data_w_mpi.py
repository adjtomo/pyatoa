#!/usr/bin/env python3
"""
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
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
from obspy import UTCDateTime

from pysep.utils.io import read_sem
from pyatoa import Config, Manager
from pyatoa import logger as ptlogger
from pyflex import logger as pflogger

# Avoid log statements which can be overly numerous
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

        dummy_time = UTCDateTime("2000-01-01")

        # Create unique event and station pairs
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

    # Partition data between the number of chosen processors
    start, stop = indices[comm.rank]

    # Misfit and Number of windows will be gathered by Rank 0, initiate empty
    sendbuf = np.empty([stop - start + 1, 3], dtype=float)

    # Main processing for each rank, read data, process, write adjoint sources
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

    # Collect all results and write to a single text file
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
