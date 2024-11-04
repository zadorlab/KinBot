from rotd_py.fragment.nonlinear import Nonlinear
from rotd_py.fragment.linear import Linear
from rotd_py.fragment.monoatomic import Monoatomic
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
temperature = generate_grid(10, 10, 1.05, 70)
energy = generate_grid(0, 10, 1.05, 190)
angular_mom = generate_grid(0, 1, 1.1, 80)

# fragment info
# Coordinates in Angstrom
{f1}
{f2}

# Setting the dividing surfaces

# Pivot_points and distances are in Bohr
divid_surf = [
{Surfaces_block}
             ]

selected_faces = {selected_faces}
faces_weights = {faces_weights}


{calc_block}

{corrections_block}

inf_energy = {inf_energy} # Hartree

kb_sample = MultiSample(fragments={frag_names}, inf_energy=inf_energy,
                         energy_size=1, min_fragments_distance={min_dist},
                         corrections=corrections,
                         name='kb_{job_name}')

# Flux parameters:
#flux_rel_err: flux accuracy (1=99% certitude, 2=98%, ...)
#pot_smp_max: maximum number of sampling for each facet
#pot_smp_min: minimum number of sampling for each facet
#tot_smp_max: maximum number of total sampling per surface
#tot_smp_min: minimum number of total sampling per surface
#smp_len: Number of valid sample asked of each subprocess

flux_parameter = {{'pot_smp_max': 1000, 'pot_smp_min': 100,
                  'tot_smp_max': 15000, 'tot_smp_min': 100,
                  'flux_rel_err': 5, 'smp_len': 1}}

flux_base = FluxBase(temp_grid=temperature,
                     energy_grid=energy,
                     angular_grid=angular_mom,
                     flux_type='MICROCANONICAL',
                     flux_parameter=flux_parameter)

# start the final run
# Will read from the restart db
multi = Multi(sample=kb_sample,
              dividing_surfaces=divid_surf,
              selected_faces=selected_faces,
              fluxbase=flux_base,
              calculator=calc)
multi.run()
multi.print_results(
dynamical_correction=1.0,
faces_weights=faces_weights)
