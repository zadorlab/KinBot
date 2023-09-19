from rotd_py.fragment.nonlinear import Nonlinear
from rotd_py.system import Surface
from rotd_py.new_multi import Multi
from rotd_py.sample.multi_sample import MultiSample
from rotd_py.flux.fluxbase import FluxBase
from ase.atoms import Atoms

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


# temperature, energy grid and angular momentum grid
temperature = generate_grid(10, 10, 1.05, 2)
energy = generate_grid(0, 10, 1.05, 2)
angular_mom = generate_grid(0, 1, 1.1, 2)

# fragment info
{Fragments_block}

# Setting the dividing surfaces
# This creates an output file, as much as the given dividing surface.
# The output file contains temperature and the flux based on the temperature.
# 
{Surfaces_block}

# how to sample the two fragments
# calc = 'amp.amp'
{calc_block}

r_inf = -79.47971696  # RS2/cc-pvtz
_{job_name} = MultiSample(fragments={frag_names}, inf_energy=r_inf,
                         energy_size=1, min_fragments_distance={min_dist})

# the flux info
{flux_block}

flux_base = FluxBase(temp_grid=temperature,
                     energy_grid=energy,
                     angular_grid=angular_mom,
                     flux_type='MICROCANONICAL',
                     flux_parameter=flux_parameter)

# start the final run
multi = Multi(sample=_{job_name}, dividing_surfaces=divid_surf,
              fluxbase=flux_base, calculator=calc)
multi.run()
# multi.total_flux['0'].flux_array[0].run(50)
# multi.total_flux['0'].save_file(0)
print(multi.total_flux['0'].flux_array[0].acct_smp())
print(multi.total_flux['0'].flux_array[0].fail_smp())
print(multi.total_flux['0'].flux_array[0].face_smp())
print(multi.total_flux['0'].flux_array[0].close_smp())

print(multi.total_flux['1'].flux_array[0].acct_smp())
print(multi.total_flux['1'].flux_array[0].fail_smp())
print(multi.total_flux['1'].flux_array[0].face_smp())
print(multi.total_flux['1'].flux_array[0].close_smp())
