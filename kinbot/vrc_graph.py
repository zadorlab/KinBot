"""Stable VRC graph identity: resource requests do not change the science.

Callers hold the worker/controller lock. Legacy graph paths are retained so
orbital checkpoint paths and completed downstream scan points remain valid.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from kinbot.execution import atomic_json


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'Ambiguous VRC specification key after JSON serialization: {key!r}')
        result[key] = value
    return result


def scientific_spec(spec):
    value = deepcopy(spec)
    value.pop('limits', None)
    from kinbot.manual_fallback import logical_task
    value['tasks'] = [logical_task(t) for t in value['tasks']]
    for task in value['tasks']:
        task.pop('resources', None)
    intent = value.get('intent', {})
    for key in ('vrc_scan_fingerprint', 'cap_walltime_to_explicit_partition'):
        intent.pop(key, None)
    if not intent:
        value.pop('intent', None)
    # Use the same representation as workflow.json. In particular, validated
    # MR models use integer multiplicity keys, but JSON objects use strings.
    # Do not round coordinates, normalize physical values or omit state data.
    return json.loads(json.dumps(value, allow_nan=False),
                      object_pairs_hook=_unique_object)


def identity(spec):
    return hashlib.sha256(json.dumps(scientific_spec(spec), sort_keys=True,
        separators=(',', ':')).encode()).hexdigest()


def _rank(item):
    """Prefer accepted work, then the earliest recorded Slurm allocation.

    Resource-only duplicates can have slightly different optimized results.
    Choosing the original accepted attempt preserves the original WFU chain.
    A persisted binding prevents later completions from changing this choice.
    """
    path, spec, state = item
    entries = list(state['tasks'].values())
    complete = bool(entries) and all(e.get('status') == 'complete' for e in entries)
    jobs = [int(e['job_id']) for e in entries if str(e.get('job_id', '')).isdecimal()]
    return (not complete, min(jobs) if jobs else float('inf'),
            (path/'workflow.json').stat().st_mtime_ns, path.name)


def binding(parent, spec):
    return Path(parent)/('.graph_binding.'+identity(spec)+'.json')


def bind(parent, spec, graph):
    parent, graph = Path(parent).resolve(), Path(graph).resolve()
    if graph.parent != parent:
        raise ValueError('VRC graph binding must remain in its graph directory.')
    atomic_json(binding(parent, spec), {'schema': 1, 'identity': identity(spec),
                                      'graph': graph.name})


def select_graph(parent, prefix, spec):
    """Find/create one graph, including resource-dependent legacy names.

    Only resources/graph concurrency and execution-policy metadata are ignored.
    Geometry, method, basis, state/root/active orbitals, checkpoint identity,
    constraints, convergence and every native input remain exact comparisons.
    Existing graph manifests and staged inputs are never rewritten here.
    """
    from kinbot.vrc_parallel import guard
    with guard(Path(parent) / '.graph_selection'):
        return _select_graph(parent, prefix, spec)


def _select_graph(parent, prefix, spec):
    from kinbot.anl.dispatch import _load, prepare
    parent = Path(parent).resolve()
    key = identity(spec)
    marker = binding(parent, spec)
    if marker.exists():
        record = json.loads(marker.read_text())
        name = record['graph']
        if (record.get('identity') != key or Path(name).name != name
                or not name.startswith(prefix+'_')):
            raise RuntimeError(f'Invalid VRC graph binding: {marker}')
        path, saved, _ = _load(parent/name)
        if scientific_spec(saved) != scientific_spec(spec):
            raise RuntimeError(f'Bound VRC graph changed scientific settings: {path}')
        return path
    matches = []
    for path in sorted(parent.glob(prefix+'_*/workflow.json')):
        saved = json.loads(path.read_text())
        if scientific_spec(saved) == scientific_spec(spec):
            matches.append(_load(path.parent))  # Require intact manifest hash.
    if matches:
        graph, _, _ = min(matches, key=_rank)
    else:
        graph = parent/(prefix+'_'+key[:16])
        if graph.exists():
            raise RuntimeError(f'VRC graph identity collision: {graph}')
        prepare(spec, graph)
    bind(parent, spec, graph)
    return graph


def select_correction_graph(parent, spec):
    """Reuse a legacy correction batch even if some outputs were copied back.

    The remaining input names/decks must all match. Never create a new batch
    for a subset of the same allocated work just because other tasks finished.
    The old absolute directory and Slurm task records remain authoritative.
    """
    from kinbot.anl.dispatch import _load, prepare
    from kinbot.vrc_parallel import guard
    parent = Path(parent).resolve()
    wanted = {t['id']: t for t in scientific_spec(spec)['tasks']}
    context = {k: v for k, v in spec['molecule'].items() if k != 'positions'}
    geometry_tokens = ('{{CARTESIAN}}', '{{XYZ}}', '{{MRCC_XYZ}}')
    uses_initial_geometry = any(token in t['input_template']
                                for t in wanted.values() for token in geometry_tokens)
    with guard(parent / '.correction_selection'):
        matches = []
        for path in sorted(parent.glob('dispatch*/workflow.json')):
            graph, saved, state = _load(path.parent)
            present = {t['id']: t for t in scientific_spec(saved)['tasks']}
            same_context = context == {k: v for k, v in saved['molecule'].items() if k != 'positions'}
            same_geometry = (not uses_initial_geometry or
                             saved['molecule']['positions'] == spec['molecule']['positions'])
            if (same_context and same_geometry
                    and all(present.get(key) == value for key, value in wanted.items())):
                matches.append((graph, saved, state))
        if matches:
            return min(matches, key=_rank)[0]
        graph = parent / 'dispatch'
        if graph.exists():
            digest = identity(spec)[:16]
            graph = parent / f'dispatch_{digest}'
        if graph.exists():
            raise RuntimeError(f'VRC correction identity collision: {graph}')
        prepare(spec, graph)
        return graph
