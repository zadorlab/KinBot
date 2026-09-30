from rotd_py.fragment.nonlinear import Nonlinear
from rotd_py.fragment.linear import Linear
from rotd_py.fragment.monoatomic import Monoatomic
from rotd_py.system import Surface
from rotd_py.new_multi import Multi
from rotd_py.sample.multi_sample import MultiSample
from rotd_py.flux.fluxbase import FluxBase
from ase.atoms import Atoms

import json
import hashlib
import os
from pathlib import Path
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
temperature = generate_grid(*{temperature_grid})
energy = generate_grid(*{energy_grid})
angular_mom = generate_grid(*{angular_grid})

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

flux_parameter = {flux_parameters}

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
dynamical_correction={dynamical_correction},
faces_weights=faces_weights)

# This manifest is written only after sampling and result generation finish.
# KinBot reads and hash-checks the MESS-facing number-of-states file and the
# native surface flux output before allowing the workflow to proceed.
result_root = Path(kb_sample.name)
result_paths = sorted(
    list(result_root.glob('Ne_*.out')) +
    list(result_root.glob('output/surface_*.dat')))
if not any(path.name.startswith('Ne_') for path in result_paths):
    raise RuntimeError('rotdPy did not write a number-of-states output.')
if not any(path.parent.name == 'output' for path in result_paths):
    raise RuntimeError('rotdPy did not write a surface flux output.')
result_files = [str(path) for path in result_paths]
result_sha256 = {{
    str(path): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in result_paths
}}
Path('{result_file}').write_text(json.dumps({{
    'schema': 2,
    'status': 'complete',
    'reaction': '{job_name}',
    'surface_count': len(getattr(multi, 'total_flux', {{}})),
    'result_files': result_files,
    'result_sha256': result_sha256,
}}, indent=2) + '\n')
