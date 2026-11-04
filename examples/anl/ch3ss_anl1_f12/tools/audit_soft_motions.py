"""Reproduce the w3/w7/w14/w15 motion audit from saved native L2 evidence.

Read only. Nothing here selects a population model or submits calculations.
Paths are explicit arguments so the script works on exports or the live run.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.db import connect
from kinbot.gaussian_modes import read_normal_modes
from kinbot.motion_character import characterize_modes
from kinbot.stationary_pt import StationaryPoint


CASES = [
    ('w3', '791671450790040000002', '_well_high'),
    ('w7', '791671560700020000002-sedd8f7729b69b19b', '_well_high'),
    ('w14', '791671560700020000002-s893815ff03c7331d', '_well_high'),
    ('w15', '791671560700020000002-sdc19d1516ced03fc', '_well_high'),
    ('rxn22', '791671560700020000002-sedd8f7729b69b19b', '_R_Addition_MultipleBond_2_1_3_high'),
    ('rxn29', '791671560700020000002-s893815ff03c7331d', '_12_shift_S_R_2_1_3_high'),
]


def audit(databases, native):
    result = {}
    for label, worker, suffix in CASES:
        job = worker + suffix
        dbpath = databases/worker/'kinbot.db'
        if not dbpath.is_file():
            raise FileNotFoundError(dbpath)
        rows = list(connect(dbpath, serial=True).select(name=job))
        row = rows[-1]
        if row.data.get('status') != 'normal':
            raise ValueError(f'{job}: latest saved result is not normally completed.')
        point = StationaryPoint(label, 0, 2, atom=row.symbols, geom=row.positions)
        point.characterize()
        observation = read_normal_modes(native/worker/(job+'.log'), row.symbols, row.positions)
        if not np.allclose(row.data['frequencies'], observation['frequencies_cm1'], rtol=0., atol=1.e-4):
            raise ValueError(f'{job}: native spectrum differs from the selected database row.')
        character = characterize_modes(point, observation['modes'])
        result[label] = dict(atoms=row.symbols, geometry=row.positions.tolist(),
            charge=0, multiplicity=2, job=job, database=str(dbpath), row_id=row.id,
            database_sha256=hashlib.sha256(dbpath.read_bytes()).hexdigest(),
            native=observation, character=character)
    return dict(schema=1, scope='Selected native L2 local motion diagnostics only',
                points=result, new_qc_calculations=0, production_models_changed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--databases', required=True, type=Path, help='CH3SS case folder with worker databases')
    parser.add_argument('--native', required=True, type=Path, help='CH3SS case folder with native Gaussian logs')
    args = parser.parse_args()
    print(json.dumps(audit(args.databases.resolve(), args.native.resolve()), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
