"""Explicit repair of an orphaned automatic graph, retaining live ownership.

This is a controller-side migration, never an automatic native-error bypass.
Only tasks using their intended geometry retain their task/input hashes.
"""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import uuid

from kinbot.execution import atomic_json, file_sha256
from kinbot.anl import dispatch
from kinbot.anl.validation import configure_task_convergence, require_primary_geometries


def descendants(spec, ident):
    found = {ident}
    while True:
        added = {t['id'] for t in spec['tasks']
                 if set(t.get('depends_on', ())).union(
                     (t.get('geometry_from'), t.get('geometry_fallback'))) & found}
        if added <= found:
            return found - {ident}
        found |= added


def geometry_invalidations(spec, state):
    """Find tasks whose native work used a fallback or a replaced geometry.

    VPT2 and explicitly declared SP@L2 diagnostics use initial by design.
    They are not descendants of the failed L3 optimization.
    """
    invalid = set()
    for task in spec['tasks']:
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        source = task.get('geometry_from', 'initial')
        actual = entry.get('geometry_from')
        wrong = source != 'initial' and actual is not None and actual != source
        if (source != 'initial' and entry.get('geometry_sha256')
                and state['tasks'].get(source, {}).get('status') == 'complete'):
            expected = state['tasks'][source].get('final_geometry_sha256')
            wrong |= expected is not None and entry['geometry_sha256'] != expected
        failed_geometry = (task['kind'] == 'ase_optimize'
                           and entry.get('status') in ('failed', 'unavailable', 'skipped', 'blocked'))
        if wrong or failed_geometry:
            invalid.add(ident)
            invalid.update(descendants(spec, ident))
    return invalid


def _rollback(graph):
    journal = graph / '.salvage_transaction.json'
    if not journal.exists():
        return
    record = json.loads(journal.read_text())
    for path in (graph / 'tasks').iterdir():
        if path.is_dir() and path.name not in record['original_directories']:
            shutil.rmtree(path)
    for ident in record['changed']:
        directory = graph / 'tasks' / ident
        saved = Path(record['archive']) / ident
        if saved.exists():
            if directory.exists():
                shutil.rmtree(directory)
            saved.rename(directory)
    for name in ('workflow.json', 'state.json'):
        shutil.copy2(Path(record['archive']) / name, graph / name)
    journal.unlink()


def repair_graph(graph, archive, *, scf=1.e-8, cc=1.e-7, l3_max_retries=1):
    """Reconcile orphans, reparse successes and stage bounded targeted retries."""
    graph = Path(graph).resolve()
    _rollback(graph)
    marker = graph / 'salvage.json'
    if marker.exists():
        previous = json.loads(marker.read_text())
        if (previous['scf'], previous['cc'], previous['l3_max_retries']) != (scf, cc, l3_max_retries):
            raise ValueError('Salvage thresholds differ from the recorded migration.')
        return previous
    # advance queries each submitted allocation before looking at a possibly
    # intermediate execution record. It never submits with submit=False.
    dispatch.advance(graph, submit=False)
    graph, spec, state = dispatch._load(graph)
    by_id = {t['id']: t for t in spec['tasks']}
    invalid = geometry_invalidations(spec, state)
    # A cancelled allocation may never publish execution.json. advance keeps
    # its submitted state during the normal result-publication grace period;
    # that state alone does not mean the allocation is still alive. Invalid
    # geometry results will be archived regardless of outcome, so they need
    # no publication wait, but they DO need a fresh successful scheduler query.
    # Valid jobs retain the ordinary grace period and their ownership.
    retired_jobs = {}
    for ident in sorted(invalid):
        entry = state['tasks'].get(ident, {})
        if entry.get('status') == 'submitting':
            raise RuntimeError(f'{graph}/{ident}: ambiguous submission must be reconciled before archive.')
        if entry.get('status') != 'submitted':
            continue
        if dispatch._job_active(entry['job_id']):
            raise RuntimeError(f'{graph}/{ident}: incompatible live job must leave the queue before archive.')
        retired_jobs[ident] = entry['job_id']
    reset = set(invalid)
    retries = set()
    for ident, entry in state['tasks'].items():
        if entry['status'] not in ('failed', 'unavailable'):
            continue
        task = by_id[ident]
        if entry.get('failure_phase') == 'preflight':
            raise RuntimeError(f'{graph}/{ident}: resolve the preflight failure before salvage.')
        if task['kind'] == 'ase_optimize':
            reset.update(descendants(spec, ident))
        elif task.get('result_parser', {}).get('kind') == 'gaussian_vpt2':
            native = graph / 'tasks' / ident / task['result_parser']['file']
            text = native.read_text(errors='replace')
            if 'Inconsistency found in framework group definition:' in text:
                # Existing narrowly checked migration archives the complete
                # failed attempt; it also handles an unassigned old group.
                dispatch.recover_gaussian_vpt2_symmetry(graph, ident)
                graph, spec, state = dispatch._load(graph)
                by_id = {t['id']: t for t in spec['tasks']}
                reset.add(ident)
                retries.add(ident)
                continue
            if 'Molecule is non-linear but NTrRo=5' not in text:
                raise RuntimeError(f'{native}: unreviewed VPT2 failure.')
            from ase.io import read
            import numpy as np
            atoms = read(native.parent / 'geometry.xyz')
            singular = np.linalg.svd(atoms.positions-atoms.positions.mean(axis=0),
                                     compute_uv=False)
            if len(atoms) < 2 or singular[1] >= 1.e-7:
                raise RuntimeError(f'{native}: linearity failure at a nonlinear geometry.')
        elif ident not in invalid:
            execution = graph / 'tasks' / ident / 'execution.json'
            phase = (json.loads(execution.read_text()).get('failure_phase')
                     if execution.exists() else entry.get('failure_phase'))
            if phase not in ('execution', 'scheduler'):
                raise RuntimeError(f'{graph}/{ident}: unreviewed validation failure.')
        retries.add(ident)
        reset.add(ident)

    protected = {n for n, e in state['tasks'].items()
                 if e['status'] in ('complete', 'submitted', 'submitting')} - invalid
    if any(e['status'] == 'submitting' for e in state['tasks'].values()):
        raise RuntimeError('Ambiguous submission must be reconciled before salvage.')
    replacements = {}
    revised = deepcopy(spec)
    require_primary_geometries(revised, preserve_tasks=protected)
    for task in revised['tasks']:
        ident = task['id']
        if ident in protected:
            continue
        configure_task_convergence(task, scf_convergence=scf, cc_convergence=cc)
        replacements[ident] = task
        if state['tasks'].get(ident, {}).get('status') == 'staged':
            reset.add(ident)
    spec['intent'] = revised['intent']
    spec['intent']['l3_max_retries'] = l3_max_retries
    archive = Path(archive) / uuid.uuid4().hex
    archive.mkdir(parents=True, exist_ok=False)
    for name in ('state.json', 'workflow.json'):
        shutil.copy2(graph / name, archive / name)
    record = {'archive': str(archive), 'changed': sorted(reset),
              'original_directories': [p.name for p in (graph / 'tasks').iterdir() if p.is_dir()]}
    atomic_json(graph / '.salvage_transaction.json', record)
    attempts = {n: state['tasks'].get(n, {}).get('attempt', 1) for n in reset}
    try:
        for ident in sorted(reset):
            directory = graph / 'tasks' / ident
            if directory.exists():
                directory.rename(archive / ident)
            state['tasks'].pop(ident, None)
        spec['tasks'] = [replacements.get(t['id'], t) for t in spec['tasks']]
        dispatch.validate_spec(spec, allow_obsolete=True)
        atomic_json(graph / 'workflow.json', spec)
        state['spec_sha256'] = file_sha256(graph / 'workflow.json')
        dispatch._clear_blocked_tasks(state)
        atomic_json(graph / 'state.json', state)
        state = dispatch.advance(graph, submit=False)
        for ident in reset:
            if ident in state['tasks']:
                state['tasks'][ident].update(
                    attempt=attempts[ident] + int(ident in retries),
                    previous_attempt=str(archive / ident))
        atomic_json(graph / 'state.json', state)
    except BaseException:
        _rollback(graph)
        raise
    result = {'schema': 1, 'scf': scf, 'cc': cc, 'l3_max_retries': l3_max_retries, 'archive': str(archive),
              'retries': sorted(retries), 'geometry_invalidations': sorted(invalid),
              'retired_jobs': retired_jobs,
              'adopted_jobs': {n: e['job_id'] for n, e in state['tasks'].items()
                               if e['status'] == 'submitted'},
              'complete': [n for n, e in state['tasks'].items() if e['status'] == 'complete']}
    (graph / '.salvage_transaction.json').unlink()
    atomic_json(marker, result)
    return result
