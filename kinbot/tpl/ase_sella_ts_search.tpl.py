from kinbot.optimization_budget import run_optimizer
import os
import numpy as np
import shutil

from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.utils import too_far

db = connect('{working_dir}/kinbot.db')
input_geom = np.asarray({geom}, dtype=float)
if (input_geom.shape != (len({atom}), 3)
        or not np.isfinite(input_geom).all()):
    raise ValueError('Reaction path input geometry is nonfinite or malformed.')
mol = Atoms(symbols={atom}, 
            positions=input_geom.copy())

const = Constraints(mol)
base_0_fix = [[idx - 1 for idx in fix] for fix in {fix}]
for fix in base_0_fix:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 3:
        const.fix_angle(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {{fix}}')

# Coordinate modification supplies an approximate starting geometry, not the
# definition of a quantitative scan. Constrain to the requested coordinate.
constraint_targets = {constraint_targets}
for target in constraint_targets:
    indices = [int(i)-1 for i in target[:-1]]
    setter = {{2: const.fix_bond, 3: const.fix_angle, 4: const.fix_dihedral}}[len(indices)]
    setter(indices, target=float(target[-1]))

kwargs = {kwargs}
if '{Code}' == 'ORCA':
    from kinbot.ase_modules.calculators.orca import OrcaProfile
    kwargs['profile'] = OrcaProfile(command=kwargs['profile'])

mol.calc = {Code}(**kwargs)
if '{Code}' == 'Gaussian':
    mol.get_potential_energy()
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

sella_kwargs = {sella_kwargs}

opt = Sella(mol, 
            order=0,
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)

# Intermediate steps don't need to fully converge.  If a constrained
# optimizer becomes singular, retain the finite coordinate-modified input;
# the final saddle must still pass its independent TS checks.
recovery_reason = None
strict_scan = {strict_scan}
try:
    converged = run_optimizer(opt, fmax={fmax} if strict_scan else 0.1,
                              steps={steps} if strict_scan else 100)
except Exception as error:
    if strict_scan:
        raise
    converged = False
    recovery_reason = f'{{type(error).__name__}}: {{error}}'
if (not np.isfinite(mol.positions).all() or too_far(mol.positions)):
    recovery_reason = recovery_reason or 'optimizer returned invalid coordinates'
if strict_scan and (not converged or recovery_reason is not None):
    raise RuntimeError('Constrained dissociation point did not converge; no scan energy accepted.')
if recovery_reason is not None:
    mol.positions = input_geom.copy()
    mol.calc.results.clear()
    with open('{label}_sella.log', 'a') as f:
        f.write('Recovered constrained path point from finite input geometry: '
                + recovery_reason + '\n')
e = mol.get_potential_energy()
if not np.isfinite(e) or not np.isfinite(mol.positions).all():
    raise RuntimeError('Reaction path recovery did not produce a finite result.')

data = {{'energy': e, 'status': 'normal'}}
if recovery_reason is not None:
    data['path_point_recovery'] = recovery_reason
db.write(mol, name='{label}', data=data)

if os.path.isdir('{label}'):
    shutil.rmtree('{label}')

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
