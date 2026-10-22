"""Queued L0 IRC validation for one unique bimolecular saddle candidate."""

import os
import pickle
import traceback

import numpy as np
from ase import Atoms

from kinbot.bimol import create_calculator, validate_candidate
from kinbot.stationary_pt import StationaryPoint


label = '{label}'
done_file = '{done_file}'
geom = np.array({geom}, dtype=float)
frag_a = StationaryPoint('frag_a', {charge_a}, {mult_a}, atom={atom_a},
                         geom=np.array({geom_a}, dtype=float))
frag_b = StationaryPoint('frag_b', {charge_b}, {mult_b}, atom={atom_b},
                         geom=np.array({geom_b}, dtype=float))
frag_a.characterize()
frag_b.characterize()
species = StationaryPoint('bimolecular', {charge}, {mult}, atom={atom},
                          geom=geom, fragA=frag_a, fragB=frag_b)

try:
    calculator = create_calculator(
        '{l0_qc}', '{calculator_module}', '{calculator_class}',
        {calculator_kwargs}, '{qc_command}', {ppn}, '{model}', '{task}',
        '{device}', compile_model={fc_compile})
    validation = validate_candidate(
        species, geom, {formed_bonds}, {broken_bonds}, calculator,
        '{workdir}', {search_kwargs})
    status = 'normal'
except Exception:
    validation = {{'valid': False, 'traceback': traceback.format_exc()}}
    status = 'error'

atoms = Atoms({atom}, positions=geom)
atoms.info.update({{'charge': {charge}, 'spin': {mult}}})
data = {{'status': status, 'validation': validation}}
pkl_file = label + '.pkl'
pkl_temporary = pkl_file + '.tmp'
with open(pkl_temporary, 'wb') as handle:
    pickle.dump({{'sym': atoms.symbols, 'pos': atoms.positions,
                 'calc': '{l0_qc}', 'name': label, 'data': data}}, handle)
    handle.flush()
    os.fsync(handle.fileno())
os.replace(pkl_temporary, pkl_file)
done_temporary = done_file + '.tmp'
with open(done_temporary, 'w') as handle:
    handle.write('done\n')
    handle.flush()
    os.fsync(handle.fileno())
os.replace(done_temporary, done_file)
