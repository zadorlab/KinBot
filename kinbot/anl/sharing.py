"""Explicit directional sharing of whole ANL graphs across PES boundaries.

The receiving owner adopts an immutable alias only before local native work
starts. The canonical graph keeps its scheduler owner and attempt history.
There is no transfer of completed flags, no component splicing, and no implicit
sharing based on a molecule name alone.
"""
from copy import deepcopy
from contextlib import contextmanager
import fcntl
import json
from pathlib import Path

from kinbot.execution import atomic_json, file_sha256
from kinbot.anl import pes_ts_reuse as reuse

CONFIG = '.kinbot_anl_sharing.json'
RECEIPT = 'cross_pes_anl_reuse.json'


@contextmanager
def locked(graph):
    """Serialize adoption with edits of its canonical scientific contract."""
    from kinbot import live_control
    control = live_control.find(graph)
    if control is None:
        yield
        return
    with (control.parent/'.kinbot_anl_sharing.lock').open('a+') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def configure(root, source, target, *, apply=False):
    from kinbot import live_control
    root = Path(root).resolve(strict=True)
    source, target = (Path(p).resolve(strict=True) for p in (source, target))
    if source == target or any(not p.is_relative_to(root) for p in (source, target)):
        raise ValueError('Choose distinct source and target PES directories under the control root.')
    if any(not (p/'pes_exploration.json').is_file() for p in (source, target)):
        raise ValueError('Both sharing endpoints must be PES owners.')
    with live_control.locked(root):
        path = root/CONFIG
        old = json.loads(path.read_text()) if path.exists() else dict(schema=1, revision=0, routes=[])
        route = dict(source=str(source), target=str(target))
        routes = old['routes'] if route in old['routes'] else [*old['routes'], route]
        edges = {}
        for row in routes:
            edges.setdefault(row['source'], set()).add(row['target'])
        def visit(node, stack):
            if node in stack:
                raise ValueError('Cross-PES sharing routes must be acyclic.')
            for other in edges.get(node, ()):
                visit(other, stack | {node})
        for node in edges:
            visit(node, set())
        result = dict(schema=1, revision=old['revision']+int(routes != old['routes']), routes=routes)
        if apply and result != old:
            atomic_json(root/'.kinbot_controls/sharing_history'/f'{result["revision"]:06d}.json', result)
            atomic_json(path, result)
        return dict(result, status='applied' if apply else 'preview',
                    adoption='owning ANL driver; only before local native execution')


def configuration(case):
    for root in (case, *case.parents):
        path = root/CONFIG
        if path.exists():
            return json.loads(path.read_text())
    return {'routes': []}


def _science(spec):
    from kinbot.manual_fallback import logical_task
    from kinbot.anl.protocol import omitted, sources
    tasks = []
    for task in spec['tasks']:
        if task['id'] in omitted(spec):
            continue
        task = deepcopy(logical_task(task))
        for key in ('resources', 'setup', 'runtime_setup'):
            task.pop(key, None)
        tasks.append(task)
    intent = deepcopy(spec.get('intent', {}))
    for key in ('live_protocol',):
        intent.pop(key, None)
    geometry, harmonic = sources(spec)
    return dict(tasks=tasks, intent=intent, geometry_source=geometry, harmonic_sources=harmonic)


def _bound(item, spec):
    """A catalog identity must describe the graph that actually ran."""
    return spec['molecule'] == item['record']['molecule']


def discover(graph):
    """Called under the receiving drive lock, before ordinary dispatch."""
    with locked(graph):
        return _discover(graph)


def _discover(graph):
    from kinbot.anl import dispatch
    graph = Path(graph).resolve()
    if (graph/'.cross_pes_sharing_disabled.json').exists():
        return None
    base = graph.parents[2]
    case = reuse.scope(base)
    if case is None:
        return None
    if not (base/RECEIPT).exists() and not configuration(case)['routes']:
        return None
    graph, spec, state = dispatch._load(graph)
    if spec.get('name') != 'anl-automatic-post-l2':
        return None
    saved = load(base)
    local = next((v for v in reuse.catalog(case).values() if Path(v['graph']).resolve() == graph), None)
    if local is None or not _bound(local, spec):
        return None
    if local['name'] in saved:
        return saved[local['name']]
    # Staged files are harmless; completed/failed/running native attempts are
    # never silently retired to favor another PES.
    if any(e.get('status') not in ('staged', 'waiting', 'blocked', 'skipped')
           or e.get('job_id') or e.get('previous_attempt') for e in state['tasks'].values()):
        return None
    sources = [Path(r['source']) for r in configuration(case)['routes'] if r['target'] == str(case)]
    limits = dict(reuse.ts_reuse.DEFAULTS, rmsd_angstrom=1.e-4)
    for source in sources:
        for item in reuse.catalog(source).values():
            candidate = Path(item['graph'])
            if not (candidate/'workflow.json').exists():
                continue
            if item['name'] in load(source/item['base']) or item['name'] in reuse.load(source/item['base']):
                continue  # Keep aliases flat, with one canonical calculation.
            _, native_spec, native_state = dispatch._load(candidate)
            if not _bound(item, native_spec) or _science(spec) != _science(native_spec):
                continue
            proof = reuse.compare(local, item, limits)
            if proof is None:
                continue
            if not any(e['status'] in ('complete', 'submitted') for e in native_state['tasks'].values()):
                continue
            # MESS must also be able to transfer the accepted atom mapping.
            import numpy as np
            a, b = (v['record']['molecule'] for v in (local, item))
            x, y = np.array(a['positions']), np.array(b['positions'])
            if reuse.ts_reuse.aligned_permutation(a['symbols'], x-x.mean(0),
                    b['symbols'], y-y.mean(0), 1.e-4) is None:
                continue
            record = dict(schema=1, receiver=local, source=item, limits=limits,
                          measurements=proof, source_science=_science(native_spec))
            p = base/RECEIPT
            data = json.loads(p.read_text()) if p.exists() else dict(schema=1, aliases={})
            data['aliases'][local['name']] = record
            atomic_json(p, data)
            return load(base)[local['name']]
    return None


def load(base):
    """Recheck scientific identity on every use, not mutable scheduler metadata."""
    from kinbot.anl import dispatch
    base = Path(base).resolve()
    path = base/RECEIPT
    if not path.exists():
        return {}
    answer = {}
    for name, row in json.loads(path.read_text())['aliases'].items():
        target = row['source']; graph = Path(target['graph'])
        current = []
        for saved in (row['receiver'], target):
            case = reuse.scope(Path(saved['graph']).parents[2])
            found = next((v for v in reuse.catalog(case).values()
                          if v['graph'] == saved['graph']), None)
            if found is None or reuse.signature(found['record']) != reuse.signature(saved['record']):
                raise ValueError('Cross-PES accepted source/receiver provenance changed.')
            current.append(found)
        _, spec, _ = dispatch._load(graph)
        _, receiver_spec, _ = dispatch._load(row['receiver']['graph'])
        if (not _bound(current[0], receiver_spec) or not _bound(current[1], spec)
                or _science(spec) != row['source_science'] or reuse.compare(
                *current, row['limits']) is None or _science(receiver_spec) != row['source_science']):
            raise ValueError('Cross-PES source science changed; review its dependent receivers.')
        answer[name] = dict(target, graph=graph, cross_pes_receipt=str(path))
    return answer


def consumers(graph):
    """Find registered receivers before allowing a canonical scientific edit."""
    from kinbot import live_control
    graph = Path(graph).resolve()
    control = live_control.find(graph)
    if control is None:
        return []
    result = []
    for route in configuration(control.parent)['routes']:
        target = Path(route['target'])
        for item in reuse.catalog(target).values():
            path = target/item['base']/RECEIPT
            if not path.exists():
                continue
            for row in json.loads(path.read_text())['aliases'].values():
                if Path(row['source']['graph']) == graph:
                    result.append(row['receiver']['graph'])
    return sorted(set(result))


def detach(graph, *, reason):
    """Receiving owner stops borrowing before applying its own protocol edit.

    The canonical graph and its allocations are untouched. The receiver's
    unexecuted graph resumes independently; the old alias remains archived.
    """
    graph = Path(graph).resolve()
    path = graph.parents[2]/RECEIPT
    if not path.exists():
        return
    data = json.loads(path.read_text())
    rows = {k:v for k,v in data['aliases'].items() if v['receiver']['graph'] == str(graph)}
    if not rows:
        return
    atomic_json(graph/'.cross_pes_sharing_disabled.json', dict(reason=reason, aliases=rows))
    data['aliases'] = {k:v for k,v in data['aliases'].items() if k not in rows}
    atomic_json(path, data)
