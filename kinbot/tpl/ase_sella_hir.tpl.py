import os
import pickle
import shutil

from ase import Atoms
from ase.io import read, write
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.{code} import {Code}

mol = Atoms(symbols={atom}, positions={geom})

const = Constraints(mol)
base_0_fix = [idx - 1 for idx in {fix}]
const.fix_dihedral(base_0_fix)

kwargs = {kwargs}
if '{Code}' == 'ORCA':
    from kinbot.ase_modules.calculators.orca import OrcaProfile
    kwargs['profile'] = OrcaProfile(command=kwargs['profile'])

basename = os.path.basename('{label}')
if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

sella_kwargs = {sella_kwargs}
sella_kwargs.setdefault('internal', True)
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False

converged = False
failure = None
try:
    mol.calc = {Code}(**kwargs)
    if '{Code}' == 'Gaussian':
        mol.get_potential_energy()
        kwargs['guess'] = 'Read'
        mol.calc = {Code}(**kwargs)
    opt = Sella(mol,
                order={order},
                constraints=const,
                trajectory='{label}.traj',
                logfile='{label}_sella.log',
                **sella_kwargs)
    converged = opt.run(fmax={fmax}, steps={steps})
except Exception as error:
    failure = error
    with open('{label}_sella.log', 'a') as stream:
        stream.write(f'Constrained optimization failed: '
                     f'{{type(error).__name__}}: {{error}}\n')

if os.path.isfile('{label}.traj'):
    try:
        trajectory = read('{label}.traj', index=':')
    except Exception as error:
        failure = failure or error
    else:
        if trajectory:
            write('{label}.xyz', trajectory, format='xyz')

if converged:
    try:
        energy = mol.get_potential_energy()
    except Exception as error:
        failure = failure or error
        data = {{'status': 'error'}}
    else:
        data = {{'energy': energy, 'status': 'normal'}}
else:
    data = {{'status': 'error'}}

# Workers publish one atomic handoff. The controller alone writes kinbot.db,
# avoiding concurrent SQLite writes from independent compute nodes.
payload = {{'sym': mol.symbols, 'pos': mol.positions,
            'calc': '{code}', 'name': '{label}', 'data': data}}
pkl_file = '{label}.pkl'
temporary = f'{{pkl_file}}.tmp.{{os.getpid()}}'
with open(temporary, 'wb') as stream:
    pickle.dump(payload, stream)
os.replace(temporary, pkl_file)

if os.path.isdir(f'{{basename}}'):
    shutil.rmtree(f'{{basename}}')

with open('{label}_sella.log', 'a') as stream:
    stream.write('done\n')
