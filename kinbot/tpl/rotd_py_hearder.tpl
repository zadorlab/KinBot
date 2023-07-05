from rotd_py.fragment.nonlinear import Nonlinear
from rotd_py.system import Surface
from rotd_py.new_multi import Multi
from rotd_py.sample.multi_sample import MultiSample
from rotd_py.flux.fluxbase import FluxBase

import numpy as np
# # Parallel debugging stuff
# import os
# from mpi4py import MPI
# import pydevd_pycharm
#
# print(os.getpid())
# size = MPI.COMM_WORLD.Get_size()
# rank = MPI.COMM_WORLD.Get_rank()
# port_mapping = [58244, 58269]
# pydevd_pycharm.settrace('localhost', port=port_mapping[rank],
#                         stdoutToServer=True, stderrToServer=True)
# if (rank == 0):
#     print(size, rank)
#     print('Master')
# else:
#     print(size, rank)
#     print('Slave')


def generate_grid(start, interval, factor, num_point):
    """Return a grid needed for the simulation of length equal to num_point

    @param start:
    @param interval:
    @param factor:
    @param num_point:
    @return:
    """
    i = 1
    grid = [start]
    for i in range(1, num_point):
        start += interval
        grid.append(start)
        interval = interval * factor
    return np.array(grid)

