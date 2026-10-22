"""Queued L0 candidate search for one bimolecular reaction channel."""

import os
import pickle
import traceback
from importlib import import_module

import numpy as np
from ase import Atoms

from kinbot.bimol import create_calculator, find_ts_candidates, _json
from kinbot.stationary_pt import StationaryPoint


label = '{label}'
done_file = '{done_file}'
conformers_a = np.array({conformers_a}, dtype=float)
conformers_b = np.array({conformers_b}, dtype=float)
frag_a = StationaryPoint('frag_a', {charge_a}, {mult_a}, atom={atom_a},
                         geom=conformers_a[0])
frag_b = StationaryPoint('frag_b', {charge_b}, {mult_b}, atom={atom_b},
                         geom=conformers_b[0])
frag_a.characterize()
frag_b.characterize()
species = StationaryPoint('bimolecular', {charge}, {mult}, atom={atom},
                          geom=np.array({reactant_geom}, dtype=float),
                          fragA=frag_a, fragB=frag_b)

try:
    # Use the controller's separated reactants for atom equivalence. Individual
    # fragment conformers can overlap when stacked without their placement.
    species.characterize()
    calculator = create_calculator(
        '{l0_qc}', '{calculator_module}', '{calculator_class}',
        {calculator_kwargs}, '{qc_command}', {ppn}, '{model}', '{task}',
        '{device}', compile_model={fc_compile})
    encounter_generator = getattr(import_module('{family_module}'),
                                  '{encounter_function}')
    status, geometries, metadata = find_ts_candidates(
        species, conformers_a, conformers_b, {instance}, {formed_bonds},
        {broken_bonds}, encounter_generator, calculator,
        os.path.join('bimol', '{reaction_name}'), {search_kwargs},
        start_after={start_after})
except Exception:
    status = -1
    geometries = []
    metadata = {{'status': 'error', 'traceback': traceback.format_exc()}}

geometry = geometries[0] if geometries else species.geom
atoms = Atoms({atom}, positions=geometry)
atoms.info.update({{'charge': {charge}, 'spin': {mult}}})
data = {{'status': 'normal' if status == 0 else 'error',
        'search': metadata}}
_json(label + '_search.json', metadata)
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
