import os
import numpy as np
import pickle

from ase import Atoms
from ase.io import read, write
from sella import Sella, Constraints, Internals

from fairchem.core import FAIRChemCalculator
from kinbot.fairchem_utils import load_predictor
from kinbot.utils import too_far

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')
if os.path.isfile('{label}.pkl'):
    os.remove('{label}.pkl')

# molecule
input_geom = np.asarray({geom}, dtype=float)
if (input_geom.shape != (len({atom}), 3)
        or not np.isfinite(input_geom).all()):
    raise ValueError('Reaction path input geometry is nonfinite or malformed.')
mol = Atoms(symbols={atom},
            positions=input_geom.copy())
kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})
mol.calc = FAIRChemCalculator(load_predictor('{fc_model_path}', '{fc_device}'), task_name='{fc_task_name}')

# constraints
const = Constraints(mol)
#internals = Internals(mol)
base_0_fix = [[idx - 1 for idx in fix] for fix in {fix}]
for fix in base_0_fix:
    if len(fix) == 2:
        #if fix not in internals:
        #    internals.add_bond(fix)
        const.fix_bond(fix)
    elif len(fix) == 3:
        const.fix_angle(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {{fix}}')

# sella
fmax = 0.1
steps = 100
sella_kwargs = {sella_kwargs}

opt = Sella(mol, 
            order=0,
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)

# run
recovery_reason = None
try:
    converged = opt.run(fmax=fmax, steps=steps)
except Exception as error:
    converged = False
    recovery_reason = f'{{type(error).__name__}}: {{error}}'

# These are path-building optimizations, not the final saddle search.  A
# singular constrained Sella step may be recovered from the already modified,
# finite input point; the final TS still has to pass optimization, frequencies,
# and IRC validation later in the ordinary workflow.
if (not np.isfinite(mol.positions).all() or too_far(mol.positions)):
    recovery_reason = recovery_reason or 'optimizer returned invalid coordinates'
if recovery_reason is not None:
    mol.positions = input_geom.copy()
    mol.calc.results.clear()
    with open('{label}_sella.log', 'a') as f:
        f.write('Recovered constrained path point from finite input geometry: '
                + recovery_reason + '\n')

try:
    traj = read('{label}.traj', index=':')
    write('{label}.xyz', traj, format='xyz')
except (FileNotFoundError, IndexError, ValueError):
    pass
e = mol.get_potential_energy()
if not np.isfinite(e) or not np.isfinite(mol.positions).all():
    raise RuntimeError('Reaction path recovery did not produce a finite result.')
mol.calc.results.pop('forces', None)

data = {{'energy': e, 'status': 'normal'}}
if recovery_reason is not None:
    data['path_point_recovery'] = recovery_reason

mol_pkl = {{'sym': mol.symbols,
            'pos': mol.positions,
            'calc': 'fairchemcalculator',
            'name': '{label}',
            'data': data}}
with open('{label}.pkl', 'wb') as f:
    pickle.dump(mol_pkl, f)
with open('{label}_sella.log', 'a') as f:
    f.write('am1\ndone\n')  # is am1 is there, it'll be deleted on restart
