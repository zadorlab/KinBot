import os
import pickle

from ase import Atoms
from ase.io import read, write
from sella import Sella, Constraints

from fairchem.core import FAIRChemCalculator
from kinbot.fairchem_utils import load_predictor

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

mol = Atoms(symbols={atom}, positions={geom})
kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})
mol.calc = FAIRChemCalculator(
    load_predictor('{fc_model_path}', '{fc_device}'),
    task_name='{fc_task_name}')

const = Constraints(mol)
base_0_fix = [idx - 1 for idx in {fix}]
const.fix_dihedral(base_0_fix)

sella_kwargs = {sella_kwargs}
sella_kwargs.setdefault('internal', True)
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
opt = Sella(mol,
            order={order},
            constraints=const,
            trajectory='{label}.traj',
            logfile='{label}_sella.log',
            **sella_kwargs)

converged = False
failure = None
try:
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
        mol.calc.results.pop('forces', None)
        data = {{'energy': energy, 'status': 'normal'}}
else:
    data = {{'status': 'error'}}

payload = {{'sym': mol.symbols, 'pos': mol.positions,
            'calc': 'fairchemcalculator', 'name': '{label}', 'data': data}}
pkl_file = '{label}.pkl'
temporary = f'{{pkl_file}}.tmp.{{os.getpid()}}'
with open(temporary, 'wb') as stream:
    pickle.dump(payload, stream)
os.replace(temporary, pkl_file)

with open('{label}_sella.log', 'a') as stream:
    stream.write('done\n')
