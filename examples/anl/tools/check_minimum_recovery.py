"""Submit an isolated L1/L2 check using the original, undisplaced worker geometry.

The production QC generator and minimum-recovery code perform all optimization,
mode descent and frequency validation. This helper changes no existing run.
It is a native validation convenience, not a required production restart step.
"""
import argparse
import ast
import json
import os
from pathlib import Path

import numpy as np

from kinbot.execution import atomic_json, file_sha256
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.run_format import ensure_current_run
from kinbot.stationary_pt import StationaryPoint


def starting_geometry(worker):
    """Read literal generated input without executing the old worker script."""
    tree = ast.parse(worker.read_text())
    for node in tree.body:
        if (isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'mol'
                for t in node.targets) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name) and node.value.func.id == 'Atoms'):
            values = {k.arg: ast.literal_eval(k.value) for k in node.value.keywords}
            symbols, positions = values['symbols'], np.asarray(values['positions'], dtype=float)
            if positions.shape != (len(symbols), 3) or not np.isfinite(positions).all():
                raise ValueError('Invalid starting geometry in worker.')
            return symbols, positions
    raise ValueError('No literal initial Atoms geometry found in the generated worker.')


def worker_state(worker):
    """Use the fragment's own electronic state, not its parent PES's state."""
    for node in ast.parse(worker.read_text()).body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == 'kwargs' for t in node.targets):
            options = ast.literal_eval(node.value)
            return int(options['charge']), int(options['mult'])
    raise ValueError('No literal charge/multiplicity found in the generated worker.')


def submit_check(input_path, worker, destination, level='l2'):
    input_path, worker = (Path(p).expanduser().resolve(strict=True) for p in (input_path, worker))
    destination = Path(destination).expanduser().resolve()
    if destination.exists():
        raise FileExistsError(f'{destination} already exists; do not submit duplicate checks.')
    symbols, positions = starting_geometry(worker)
    charge, multiplicity = worker_state(worker)
    parameters = Parameters(str(input_path)).par
    species = StationaryPoint('minimum_check', charge, multiplicity,
                              atom=symbols, geom=positions, natom=len(symbols))
    species.characterize()
    previous = Path.cwd()
    ensure_current_run(destination, create=True)
    try:
        os.chdir(destination)
        qc = getattr(QuantumChemistry(parameters), level)
        if qc is None or not (qc.qc == 'gauss' and qc.use_sella
                              and qc.frequency_mode == 'native_hessian'):
            raise ValueError('This check requires Gaussian/Sella/native_hessian.')
        if qc.par['minimum_recovery_attempts'] <= 0:
            raise ValueError('Minimum recovery is disabled in the input.')
        atomic_json(destination/'input.json', parameters)
        record = {'source_input': str(input_path), 'source_worker': str(worker),
                  'source_worker_sha256': file_sha256(worker),
                  'starting_geometry_angstrom': positions.tolist(), 'symbols': list(symbols),
                  'level': level, 'charge': charge, 'multiplicity': multiplicity,
                  'minimum_recovery_attempts': qc.par['minimum_recovery_attempts'],
                  'minimum_recovery_displacement': qc.par['minimum_recovery_displacement'],
                  'status': 'preparing'}
        atomic_json(destination/'minimum_check.json', record)
        qc.qc_opt(species, positions, high_level=int(level == 'l2'))
        record.update(status='submitted', job_ids=qc.job_ids)
        atomic_json(destination/'minimum_check.json', record)
        return record
    finally:
        os.chdir(previous)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('worker', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--level', choices=('l1', 'l2'), default='l2')
    args = parser.parse_args()
    print(json.dumps(submit_check(args.input, args.worker, args.destination, args.level), indent=2))
