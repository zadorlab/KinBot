"""Explicit, journalled retirement of resource-only VRC duplicate chains.

Used inside run_refresh's stopped-controller/worker-lock transaction. No ANL
graph or unrelated allocation is selected. Original graph paths are retained.
"""
import json
from pathlib import Path

from kinbot.execution import atomic_json, file_sha256
from kinbot.vrc_graph import identity, _rank, bind


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from _strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _strings(item)


def verify_complete(graph):
    from kinbot.anl import dispatch
    graph, spec, state = dispatch._load(graph)
    for task in spec['tasks']:
        entry = state['tasks'].get(task['id'], {})
        if entry.get('status') != 'complete':
            raise RuntimeError(f'{graph}: original calculation is no longer complete.')
        dispatch._verify_stage_files(graph, task, entry)
        execution = json.loads((graph/'tasks'/task['id']/'execution.json').read_text())
        if execution.get('status') != 'executed':
            raise RuntimeError(f'{graph}: original calculation has no successful execution.')
        actual = dispatch._verify_execution(graph, task, entry, execution)
        if actual != entry.get('final_geometry_sha256'):
            raise RuntimeError(f'{graph}: original geometry provenance changed.')
    return spec


def plan(root, jobs):
    from kinbot.anl import dispatch
    root = Path(root).resolve()
    retired = set()
    marker = root/'vrc_restart.json'
    if marker.exists():
        retired = {str(Path(x['graph']).resolve()) for x in
                   json.loads(marker.read_text()).get('retired_scans', [])}
    records, groups = {}, {}
    for folder in ('molpro_scan', 'molpro_pivots'):
        for p in sorted(root.glob(f'*/[0-9]*/vrctst/{folder}/*/workflow.json')):
            if str(p.parent.resolve()) in retired:
                continue
            if p.parent.is_symlink() or not p.parent.resolve().is_relative_to(root):
                raise RuntimeError(f'VRC graph is outside the run: {p.parent}')
            item = dispatch._load(p.parent)
            path, spec, state = item
            records[str(path)] = item
            key = (str(path.parent), path.name.rsplit('_', 1)[0], identity(spec))
            groups.setdefault(key, []).append(item)
    keep, remove, reasons = {}, set(), {}
    for items in groups.values():
        if len(items) < 2:
            continue
        original = min(items, key=_rank)
        entries = original[2]['tasks'].values()
        if not all(e.get('status') == 'complete' for e in entries):
            raise RuntimeError(f'{original[0]}: duplicate group has no completed '
                               'original to verify; review before retiring attempts.')
        path, spec, _ = original
        keep[str(path)] = file_sha256(path/'workflow.json')
        for other, _, _ in items:
            if other != path:
                remove.add(str(other))
                reasons[str(other)] = f'Resource-only duplicate of {path}'
    # Continue through checkpoint dependencies and the fixed-fragment pivot
    # provenance records. Close-coordinate comparison is NEVER a reuse rule.
    artifacts = {}
    evidence = []
    for p in sorted(root.glob('*/[0-9]*/vrctst/pivot_cubes/**/*.orbitals.json')):
        data = json.loads(p.read_text())
        evidence.append((p, data))
    def refers(text, graph):
        return graph+'/' in text or text == graph
    while True:
        before = set(remove)
        for name, (_, spec, _) in records.items():
            if name in remove:
                continue
            source = next((g for g in remove if any(refers(s, g) for s in _strings(spec))), None)
            if source:
                keep.pop(name, None)
                remove.add(name)
                reasons[name] = f'Depends on retired duplicate {source}'
        for path, data in evidence:
            source = str(Path(data.get('geometry_source', {}).get('dispatch', '/')).resolve())
            graph = str(Path(data.get('run_dir', '/')).resolve())
            if source in remove:
                if graph not in records:
                    raise RuntimeError(f'{path}: duplicate-derived pivot graph is unavailable.')
                if graph not in keep:
                    remove.add(graph)
                    reasons[graph] = f'Pivot geometry came from retired duplicate {source}'
                artifacts[str(path)] = file_sha256(path)
                cube = path.with_suffix('').with_suffix('.cube')
                if cube.exists():
                    artifacts[str(cube)] = file_sha256(cube)
        if remove == before:
            break
    if remove & set(keep):
        raise RuntimeError('Duplicate closure would retire an original calculation.')
    for parent in {records[g][0].parent for g in remove}:
        for path in parent.glob('.graph_binding.*.json'):
            record = json.loads(path.read_text())
            if str((parent/record['graph']).resolve()) in remove:
                artifacts[str(path)] = file_sha256(path)
    # Corrections or sampling already built from an obsolete chain require a
    # broader scientific migration, not this targeted scheduling repair.
    workers = {records[g][0].parent.parent.parent for g in remove}
    for worker in workers:
        if list((worker/'vrctst').glob('corr_*.json')):
            raise RuntimeError(f'{worker}: correction potentials already exist; review their provenance first.')
        for p in (worker/'vrctst/molpro').glob('*/workflow.json'):
            if any(refers(s, g) for s in _strings(json.loads(p.read_text())) for g in remove):
                raise RuntimeError(f'{p}: correction work depends on a duplicate; explicit review required.')
        for p in (worker/'vrctst').glob('*.log'):
            if p.stat().st_size < 8192 and any(g in p.read_text() for g in remove):
                artifacts[str(p)] = file_sha256(p)
    selected = []
    live = {str(Path(j['directory']).resolve()): j for j in jobs}
    for name in sorted(remove):
        graph, spec, state = records[name]
        for task in spec['tasks']:
            entry = state['tasks'].get(task['id'], {})
            if entry.get('status') == 'submitting':
                raise RuntimeError(f'{graph}: uncertain submission must be reconciled before retirement.')
            job = live.get(str(graph/'tasks'/task['id']))
            if job:
                if str(entry.get('job_id')) != str(job['job_id']):
                    raise RuntimeError(f'{graph}: duplicate job ownership changed.')
                selected.append(job)
    for graph in keep:
        verify_complete(graph)
    return {'keep': keep, 'remove': {g: {'workflow_sha256': file_sha256(Path(g)/'workflow.json'),
                                       'reason': reasons[g]} for g in sorted(remove)},
            'artifacts': artifacts, 'selected_jobs': selected}


def apply(root, record, archive, journal, marker):
    """Validate, cancel exact duplicate IDs, then move results into archive."""
    from kinbot.vrc_restart import cancel_scans
    root, archive = Path(root).resolve(), Path(archive).resolve()
    committed = record.get('phase') == 'archived'
    for name, expected in record['keep'].items():
        graph = Path(name)
        if file_sha256(graph/'workflow.json') != expected:
            raise RuntimeError(f'{graph}: original manifest changed during repair.')
        verify_complete(graph)
    # Check every move before cancelling anything; support interrupted moves.
    moves = []
    for name, details in record['remove'].items():
        path = Path(name)
        if not path.is_relative_to(root) or path.is_symlink():
            raise RuntimeError(f'Unsafe duplicate graph: {path}')
        destination = archive/path.relative_to(root)
        actual = destination if committed else path if path.exists() else destination
        if file_sha256(actual/'workflow.json') != details['workflow_sha256']:
            raise RuntimeError(f'{path}: duplicate manifest changed during repair.')
        moves.append((path, destination))
    for name, expected in record['artifacts'].items():
        path = Path(name)
        if not path.is_relative_to(root) or path.is_symlink():
            raise RuntimeError(f'Unsafe duplicate artifact: {path}')
        destination = archive/path.relative_to(root)
        actual = destination if committed else path if path.exists() else destination
        if file_sha256(actual) != expected:
            raise RuntimeError(f'{path}: duplicate artifact changed during repair.')
        moves.append((path, destination))
    cancel_scans(record, journal, marker)
    # No surviving allocation may still write into a directory being moved.
    from kinbot.run_refresh import _owned_jobs
    for job in _owned_jobs(root):
        work = Path(job['directory']).resolve()
        if any(work.is_relative_to(Path(g)) for g in record['remove']):
            raise RuntimeError(f'Duplicate allocation {job["job_id"]} is still present; repeat refresh.')
    for path, destination in moves:
        if path.exists() and not committed:
            if destination.exists():
                raise RuntimeError(f'Archive destination already exists: {destination}')
            destination.parent.mkdir(parents=True, exist_ok=True)
            path.rename(destination)
    previous = root/'vrc_deduplication.json'
    if not committed and previous.exists() and not (archive/'previous_deduplication.json').exists():
        atomic_json(archive/'previous_deduplication.json', json.loads(previous.read_text()))
    record['phase'] = 'archived'
    record['archive'] = str(archive)
    atomic_json(marker, journal)
    for name in record['keep']:
        graph = Path(name)
        spec = json.loads((graph/'workflow.json').read_text())
        bind(graph.parent, spec, graph)
    atomic_json(previous, record)
