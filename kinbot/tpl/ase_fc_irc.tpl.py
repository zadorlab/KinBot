import os
import numpy as np
import pickle
import shutil

from ase import Atoms
from ase.io import read, write

from sella import Sella, IRC

from fairchem.core import FAIRChemCalculator
from kinbot.fairchem_utils import load_predictor


def finite_geometry(atoms):
    values = np.asarray(atoms.positions, dtype=float)
    return (values.shape == (len(atoms), 3)
            and np.isfinite(values).all())


def finite_frames(path):
    try:
        frames = read(path, index=':')
    except Exception:
        return []
    return [frame for frame in frames if finite_geometry(frame)]


if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# molecule
mol = Atoms(symbols={atom},
            positions={geom})
initial_positions = mol.positions.copy()
kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})
mol.calc = FAIRChemCalculator(load_predictor('{fc_model_path}', '{fc_device}'), task_name='{fc_task_name}')

# irc
irc = IRC(mol, 
         trajectory='{label}.traj', 
         dx=0.1, 
         eta=1e-4, 
         gamma=0, 
         keep_going=True,
         logfile='{label}_sella.log')
if '{label}'.endswith('F'):
    direction = 'forward'
elif '{label}'.endswith('R'):
    direction = 'reverse'

# run
try:
    converged_irc = irc.run(fmax=0.01, steps=100, direction=direction)
except Exception as exc:
    converged_irc = False
    with open('{label}_sella.log', 'a') as f:
        f.write(f'IRC integration raised: {{exc}}\n')
# write xyz
traj = finite_frames('{label}.traj')
if traj:
    write('{label}.xyz', traj, format='xyz')
    if not finite_geometry(mol):
        mol.positions = traj[-1].positions.copy()
# An IRC exception can leave the unchanged transition state in ``mol``.  That
# is finite, but it is not an endpoint and must never be optimized as one.
endpoint_moved = (finite_geometry(mol)
                  and np.linalg.norm(mol.positions - initial_positions) > 1e-4)
endpoint_valid = endpoint_moved and (converged_irc or len(traj) > 1)

# save results
if endpoint_valid:
    try:
        e = mol.get_potential_energy()
        data = {{'energy': e, 'status': 'normal'}}
        mol.calc.results.pop('forces', None)
    except Exception as exc:
        data = {{'status': 'error', 'error': str(exc)}}
else:
    mol.positions = initial_positions
    data = {{'status': 'error',
            'error': 'IRC produced no distinct finite endpoint.'}}
mol_pkl = {{'sym': mol.symbols,
            'pos': mol.positions,
            'calc': 'fairchemcalculator',
            'name': '{label}',
            'data': data}}
with open('{label}.pkl', 'wb') as f:
    pickle.dump(mol_pkl, f)
with open('{label}_sella.log', 'a') as f:
    f.write('done\n')

# product optimization
#prod_kwargs = {prod_kwargs}

if data['status'] == 'normal':
    # sella
    sella_kwargs = {sella_kwargs}
    sella_kwargs.setdefault('internal', True)
    if sella_kwargs['internal'] and len(mol.symbols) < 5:
        sella_kwargs['internal'] = False
    opt = Sella(mol,
                order=0,
                trajectory='{label}_prod.traj',
                logfile='{label}_prod_sella.log',
                **sella_kwargs)

    try:
        converged_opt = opt.run(fmax={fmax}, steps={steps})
    except Exception as exc:
        converged_opt = False
        with open('{label}_prod_sella.log', 'a') as f:
            f.write(f'Endpoint optimization raised: {{exc}}\n')
    prod_traj = finite_frames('{label}_prod.traj')
    if prod_traj:
        write('{label}_prod.xyz', prod_traj, format='xyz')
        if not finite_geometry(mol):
            mol.positions = prod_traj[-1].positions.copy()
    if finite_geometry(mol):
        try:
            e = mol.get_potential_energy()
            mol.calc.results.pop('forces', None)
            # IRC endpoint identification only reads this optimized geometry.
            # Product frequencies are calculated later for the accepted
            # product at its requested level, so a numerical L1 Hessian here
            # is redundant and can dominate a short scheduler allocation.
            data = {{'energy': e, 'status': 'normal'}}
        except Exception as exc:
            data = {{'status': 'error', 'error': str(exc)}}
    else:
        data = {{'status': 'error',
                'error': 'Endpoint optimization produced no finite geometry.'}}

    mol_pkl = {{'sym': mol.symbols,
                'pos': mol.positions,
                'calc': 'fairchemcalculator',
                'name': '{label}_prod',
                'data': data}}
    with open('{label}_prod.pkl', 'wb') as f:
        pickle.dump(mol_pkl, f)

if os.path.isdir('{label}_prod_vib'):
    shutil.rmtree('{label}_prod_vib')

if data['status'] == 'normal' or os.path.isfile('{label}_prod_sella.log'):
    with open('{label}_prod_sella.log', 'a') as f:
        f.write('done\n')
