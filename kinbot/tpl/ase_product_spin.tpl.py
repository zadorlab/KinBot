"""Single-point energy of one product multiplicity at the IRC geometry."""
import os
import traceback
import numpy as np
from ase import Atoms
from ase.db import connect
from kinbot.bimol import create_calculator

label = '{label}'
mol = Atoms({atom}, positions={geom})
mol.info.update(charge={charge}, spin={mult})
data = {{'backend': '{qc}', 'method': {method!r}, 'basis': {basis!r},
         'charge': {charge}, 'multiplicity': {mult}}}
try:
    mol.calc = create_calculator(
        '{qc}', 'kinbot.ase_modules.calculators.{code}', '{Code}',
        {kwargs}, '{qc_command}', {ppn}, {model!r}, {task!r}, {device!r},
        compile_model={fc_compile})
    energy = mol.get_potential_energy()
    if not np.isfinite(energy):
        raise ValueError('Nonfinite single-point energy')
    data.update(energy=float(energy), status='normal')
except Exception:
    traceback.print_exc()
    data.update(status='error', error=traceback.format_exc())
connect('{working_dir}/kinbot.db').write(mol, name=label, data=data)
with open(label + '.done.tmp', 'w') as handle:
    handle.write('done\n')
os.replace(label + '.done.tmp', label + '.done')
