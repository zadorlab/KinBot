"""Continue verified, step-exhausted ASE tasks without changing graph identity."""
from copy import deepcopy
import json
from pathlib import Path
import uuid

import numpy as np
from ase.io import read, write

from kinbot.execution import atomic_json, file_sha256
from kinbot.optimization_budget import step_retries


def _last_point(directory, task, original):
    trajectory = directory/'optimization.traj'
    rows = [s.split() for s in (directory/'optimization.log').read_text().splitlines()
            if s.startswith(('Sella ', 'BFGS:'))]
    if not rows:
        raise RuntimeError('Optimizer continuation lacks a step record.')
    row = rows[-1]
    index = int(row[-1]) if row[0] == 'Sella' and len(row) == 9 else -1
    first, last = read(trajectory, index=0), read(trajectory, index=index)
    if (first.get_chemical_symbols() != original.get_chemical_symbols()
            or last.get_chemical_symbols() != original.get_chemical_symbols()
            or not np.isfinite(last.positions).all()):
        raise RuntimeError('Optimizer trajectory atom identity/coordinates are invalid.')
    # A previous continuation starts at its explicitly bound seed.
    record = json.loads((directory/'task.json').read_text())
    expected = original.copy()
    load_seed(directory, record, expected)
    from kinbot.anl.optimization_progress import restore_attempt_start
    restore_attempt_start(directory,record,expected)
    if not np.allclose(first.positions, expected.positions, atol=1.e-6, rtol=0):
        raise RuntimeError('Optimizer trajectory does not start at its staged geometry.')
    energy, forces = last.get_potential_energy(), last.get_forces()
    if not np.isfinite(energy) or not np.isfinite(forces).all():
        raise RuntimeError('Optimizer continuation needs finite evaluated energy/forces.')
    if int(row[1]) < task.get('optimizer', {}).get('steps', 100):
        raise RuntimeError('Optimizer did not reach its stated step budget.')
    if abs(float(row[3]) - energy) > 2.e-6:
        raise RuntimeError('Optimizer final trajectory and log energies disagree.')
    return last


def load_seed(directory, record, atoms):
    seed = record.get('optimizer_continuation')
    if not seed:
        return 0
    path = directory/'continuation.xyz'
    if file_sha256(path) != seed['sha256']:
        raise RuntimeError('Optimizer continuation seed changed.')
    saved = read(path)
    if (saved.get_chemical_symbols() != atoms.get_chemical_symbols()
            or not np.isfinite(saved.positions).all()):
        raise RuntimeError('Invalid optimizer continuation seed.')
    atoms.positions = saved.positions
    return seed['chunks_used']


def stage_seed(graph, directory, task, state, atoms):
    seed = state.get('optimizer_continuations', {}).get(task['id'])
    if not seed:
        return None
    saved = Path(seed['source'])
    if file_sha256(saved) != seed['sha256']:
        raise RuntimeError('Archived optimizer continuation seed changed.')
    (directory/'continuation.xyz').write_bytes(saved.read_bytes())
    record = dict(seed)
    load_seed(directory, {'optimizer_continuation': record}, atoms.copy())
    return record


def resume(graph):
    marker = Path(graph)/'optimization_recovery.json'
    if marker.exists():
        record = json.loads(marker.read_text())
        if record['phase'] == 'preparing':
            from kinbot.run_recover import _apply_graph
            _apply_graph(record['plan'], Path(graph).resolve(), Path(record['archive']))
            record['phase'] = 'complete'
            atomic_json(marker, record)


def recover(graph, spec, state):
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id
    tasks = {t['id']: t for t in spec['tasks']}
    selected, seeds = [], {}
    for ident, entry in state['tasks'].items():
        task = tasks[ident]
        budget = 1 + task.get('optimizer', {}).get('step_retries', step_retries(graph))
        if entry['status'] not in ('failed', 'unavailable') or task['kind'] != 'ase_optimize':
            continue
        directory = graph/'tasks'/ident
        outcome = directory/'execution.json'
        if not outcome.is_file():
            continue
        execution = json.loads(outcome.read_text())
        if execution.get('error') != 'ASE geometry optimization did not converge.':
            continue
        marker = directory/'optimization_budget.json'
        prior = json.loads(marker.read_text()) if marker.exists() else {'chunks_used': 1}
        used = prior['chunks_used']
        if type(used) is not int or used < 1:
            raise RuntimeError('Invalid optimizer continuation budget record.')
        if used >= budget:
            continue
        job = _job_id(directory, ident, entry)
        if not job or dispatch._job_active(job):
            continue
        dispatch._verify_stage_files(graph, task, entry)
        dispatch._verify_execution(graph, task, entry, execution)
        if any(state['tasks'].get(n, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged') for n in descendants(spec, ident)):
            raise RuntimeError(f'{ident}: optimizer continuation has downstream work; review provenance.')
        original = read(directory/'geometry.xyz')
        last = _last_point(directory, task, original)
        seeds[ident] = (last, used, file_sha256(directory/'optimization.traj'))
        selected.append(ident)
    if not selected:
        return state
    archive = graph/'attempts'/('optimization_'+uuid.uuid4().hex)
    archive.mkdir(parents=True)
    updated = deepcopy(state)
    for ident, (last, used, trajectory_hash) in seeds.items():
        path = archive/(ident+'.xyz')
        write(path, last)
        updated.setdefault('optimizer_continuations', {})[ident] = dict(
            source=str(path), sha256=file_sha256(path), chunks_used=used,
            source_trajectory_sha256=trajectory_hash)
    plan = dict(graph=str(graph), spec=spec, state=updated, reset=selected,
                before={n: file_sha256(graph/n) for n in ('workflow.json', 'state.json')},
                reasons={n: 'Continue from verified last evaluated geometry after step exhaustion.' for n in selected},
                complete=[n for n,e in state['tasks'].items() if e['status']=='complete'],
                active={n:e['job_id'] for n,e in state['tasks'].items()
                        if e['status']=='submitted' and n not in selected})
    marker = graph/'optimization_recovery.json'
    if marker.exists():
        atomic_json(archive/'previous_recovery.json', json.loads(marker.read_text()))
    atomic_json(marker, dict(phase='preparing', archive=str(archive), plan=plan))
    resume(graph)
    return dispatch._load(graph)[2]
