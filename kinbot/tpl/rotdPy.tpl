from rotd_py.fragment.nonlinear import Nonlinear
from rotd_py.system import Surface
from rotd_py.new_multi import Multi
from rotd_py.sample.multi_sample import MultiSample
from rotd_py.flux.fluxbase import FluxBase
from ase.atoms import Atoms

import numpy as np

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
temperature = generate_grid(10, 10, 1.05, 10)
energy = generate_grid(0, 10, 1.05, 10)
angular_mom = generate_grid(0, 1, 1.1, 10)

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

{scan_ref}
{scan_trust}
{scan_sample}

inf_energy = {inf_energy}

_{job_name} = MultiSample(fragments={frag_names}, inf_energy=r_inf,
                         energy_size=1, min_fragments_distance={min_dist},
                         x_sample=x_sample, y_sample=y_sample,
                         x_trust=x_trust, y_trust=y_trust, scan_ref=scan_ref)

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
