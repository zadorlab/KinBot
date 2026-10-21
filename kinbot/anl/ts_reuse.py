"""Conservative, opt-in sharing of indistinguishable accepted L2 saddles.

An alias consumes one entire canonical graph. Components from independently
optimized geometries are never spliced together. Existing graphs are retired
only by the explicit stopped-controller migration in network_revisit.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
from itertools import islice

import numpy as np

from kinbot.execution import atomic_json


DEFAULTS = dict(rmsd_angstrom=5.e-4, energy_hartree=1.e-6,
                zpe_hartree=1.e-5, frequency_cm_inverse=1.)


def options(value):
    if not isinstance(value, dict) or set(value) - set(DEFAULTS):
        raise ValueError('anl_ts_reuse must be an object of TS comparison tolerances.')
    result = {**DEFAULTS, **value}
    for k, v in result.items():
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) or v <= 0:
            raise ValueError(f'anl_ts_reuse.{k} must be finite and positive.')
    return result


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def directory(base, name):
    return Path(base)/'states'/('ts_'+hashlib.sha256(name.encode()).hexdigest()[:16])/'composite'


def aligned_permutation(symbols_x, px, symbols_y, py, tolerance):
    """Element/distance graph mapping, followed by a proper-rotation RMSD.

    Distance matching alone would merge enantiomers. The final determinant +1
    alignment prevents that. Cap enumeration conservatively: missing a match
    never grants reuse. Isotope-bearing states are rejected by the caller.
    """
    import networkx as nx
    if sorted(symbols_x) != sorted(symbols_y):
        return None
    cutoff = 2 * np.sqrt(len(px)) * tolerance + 1.e-12
    def graph(symbols, positions):
        g = nx.Graph()
        g.add_nodes_from((i, {'element': s}) for i,s in enumerate(symbols))
        g.add_edges_from((i,j,{'distance':float(np.linalg.norm(positions[i]-positions[j]))})
                         for i in range(len(positions)) for j in range(i))
        return g
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph(symbols_x,px), graph(symbols_y,py),
        node_match=lambda a,b:a['element']==b['element'],
        edge_match=lambda a,b:abs(a['distance']-b['distance'])<=cutoff)
    for mapping in islice(matcher.isomorphisms_iter(), 10000):
        permutation = [mapping[i] for i in range(len(px))]
        target = py[permutation]
        u, _, vt = np.linalg.svd(px.T@target)
        correction = np.eye(3)
        correction[-1,-1] = np.linalg.det(u@vt)
        rmsd = float(np.sqrt(np.mean(np.sum((px@u@correction@vt-target)**2,axis=1))))
        if rmsd <= tolerance:
            return rmsd, permutation
    return None


def equivalent(plan, left, right, limits, *, stationary_point='transition_state'):
    """Return measured evidence, or None. Never use endpoints alone as identity."""
    reactions = {r['transition_state']: r for r in plan['network']['reactions']
                 if r.get('transition_state')}
    a, b = reactions[left], reactions[right]
    endpoints = lambda r: sorted(((r['reactant'],), tuple(sorted(r['products']))))
    if endpoints(a) != endpoints(b):
        return None
    if a.get('path_class') != b.get('path_class'):
        return None
    x, y = (plan['network']['transition_states'][n] for n in (left, right))
    if any(x.get(k) != y.get(k) for k in ('l2_theory', 'highest_tier', 'fallback_tiers')):
        return None
    # Legacy records omit isotope/site metadata. Do not infer it from formula.
    if any(re.search(r'\[\d', z['smiles']) for z in (x,y)):
        return None
    mx, my = x['molecule'], y['molecule']
    if any(mx.get(k) != my.get(k) for k in ('charge', 'multiplicity')):
        return None
    px, py = (np.asarray(m['positions'], dtype=float) for m in (mx, my))
    if px.shape != py.shape or px.shape != (len(mx['symbols']), 3):
        return None
    if not np.isfinite(px).all() or not np.isfinite(py).all():
        return None
    px, py = px-px.mean(axis=0), py-py.mean(axis=0)
    alignment = aligned_permutation(mx['symbols'],px,my['symbols'],py,limits['rmsd_angstrom'])
    if alignment is None:
        return None
    rmsd, permutation = alignment
    fx, fy = (np.asarray(z.get('frequencies_cm_inverse', []), float) for z in (x, y))
    if not len(fx) or fx.shape != fy.shape or not np.isfinite(fx).all() or not np.isfinite(fy).all():
        return None
    order = 1 if stationary_point == 'transition_state' else 0
    if np.count_nonzero(fx < 0) != order or np.count_nonzero(fy < 0) != order:
        return None
    measurements = dict(rmsd_angstrom=rmsd,
        energy_hartree=abs(x['energy_hartree']-y['energy_hartree']),
        zpe_hartree=abs(x['harmonic_zpe_hartree']-y['harmonic_zpe_hartree']),
        frequency_cm_inverse=float(np.max(np.abs(np.sort(fx)-np.sort(fy)))))
    if any(not np.isfinite(v) or v > limits[k] for k, v in measurements.items()):
        return None
    measurements['atom_permutation'] = permutation
    return measurements


def discover(plan, base, limits, previous=None, pinned=()):
    """Prefer the graph with most completed tasks, then most retained allocations."""
    records = plan['network']['transition_states']
    result = deepcopy(previous or dict(schema=1, aliases={}, evidence={}, records={}))
    for name, signature in result['records'].items():
        if name not in records or digest(records[name]) != signature:
            raise ValueError(f'{name}: accepted TS changed since reuse was established.')
    def rank(name):
        p = directory(base, name)/'state.json'
        entries = json.loads(p.read_text()).get('tasks', {}).values() if p.exists() else []
        return (-sum(e.get('status') == 'complete' for e in entries),
                -sum(e.get('status') == 'submitted' for e in entries), -int(p.exists()), name)
    canonical = sorted(set(result['aliases'].values()) | set(pinned), key=rank)
    for name in sorted(records, key=rank):
        if name in result['aliases'] or name in canonical:
            continue
        match = next(((other, evidence) for other in canonical
                      if (evidence := equivalent(plan, name, other, limits)) is not None), None)
        if match:
            other, evidence = match
            result['aliases'][name] = other
            result['evidence'][name] = dict(canonical=other, measurements=evidence,
                                           limits=limits, atom_mapping='element_permutation_proper_rotation')
        else:
            canonical.append(name)
    result['records'] = {n: digest(records[n]) for n in
                         set(result['aliases']) | set(result['aliases'].values())}
    return result


def load(plan, base):
    p = Path(base)/'ts_reuse.json'
    if not p.exists():
        return {}
    result = json.loads(p.read_text())
    records = plan['network']['transition_states']
    for n, signature in result['records'].items():
        if n not in records or digest(records[n]) != signature:
            raise ValueError(f'{n}: TS reuse provenance differs from accepted network.')
    aliases = result['aliases']
    if any(k == v or v in aliases for k, v in aliases.items()):
        raise ValueError('TS reuse must be a flat map of aliases to canonical saddles.')
    return aliases


def install_new(plan, base, value):
    """Share not-yet-created graphs; existing allocations require explicit retirement."""
    path=Path(base)/'ts_reuse.json'
    previous=json.loads(path.read_text()) if path.exists() else None
    result=discover(plan,base,options(value),previous,protected(base))
    for alias in set(result['aliases'])-set((previous or {}).get('aliases',{})):
        if directory(base,alias).exists():
            raise ValueError('TS graph already exists; use --deduplicate-ts to retire it with scheduler verification.')
    atomic_json(path,result)
    return result


def protected(base):
    """PES ownership is immutable even when a later local duplicate appears."""
    from kinbot.anl.pes_ts_reuse import scope, NAME
    case=scope(base)
    if case is None or not (case/NAME).exists():
        return ()
    record=json.loads((case/NAME).read_text())
    prefix=str(Path(base).resolve().relative_to(case))+'/'
    return tuple(key[len(prefix):] for key in record['records']
                 if key.startswith(prefix) and not key[len(prefix):].startswith(('stable:', 'complex:')))
