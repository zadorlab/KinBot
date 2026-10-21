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


def equivalent(plan, left, right, limits):
    """Return measured evidence, or None. Never use endpoints alone as identity."""
    reactions = {r['transition_state']: r for r in plan['network']['reactions']
                 if r.get('transition_state')}
    a, b = reactions[left], reactions[right]
    if (a['reactant'], sorted(a['products'])) != (b['reactant'], sorted(b['products'])):
        return None
    if a.get('path_class') != b.get('path_class'):
        return None
    x, y = (plan['network']['transition_states'][n] for n in (left, right))
    if any(x.get(k) != y.get(k) for k in ('l2_theory', 'highest_tier', 'fallback_tiers', 'smiles')):
        return None
    # Legacy records omit isotope/site metadata. Do not infer it from formula.
    if re.search(r'\[\d', x['smiles']):
        return None
    mx, my = x['molecule'], y['molecule']
    if any(mx.get(k) != my.get(k) for k in ('symbols', 'charge', 'multiplicity')):
        return None
    px, py = (np.asarray(m['positions'], dtype=float) for m in (mx, my))
    if px.shape != py.shape or px.shape != (len(mx['symbols']), 3):
        return None
    if not np.isfinite(px).all() or not np.isfinite(py).all():
        return None
    px, py = px-px.mean(axis=0), py-py.mean(axis=0)
    u, _, vt = np.linalg.svd(px.T@py)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u@vt)
    rmsd = float(np.sqrt(np.mean(np.sum((px@u@correction@vt-py)**2, axis=1))))
    fx, fy = (np.asarray(z.get('frequencies_cm_inverse', []), float) for z in (x, y))
    if not len(fx) or fx.shape != fy.shape or not np.isfinite(fx).all() or not np.isfinite(fy).all():
        return None
    if np.count_nonzero(fx < 0) != 1 or np.count_nonzero(fy < 0) != 1:
        return None
    measurements = dict(rmsd_angstrom=rmsd,
        energy_hartree=abs(x['energy_hartree']-y['energy_hartree']),
        zpe_hartree=abs(x['harmonic_zpe_hartree']-y['harmonic_zpe_hartree']),
        frequency_cm_inverse=float(np.max(np.abs(np.sort(fx)-np.sort(fy)))))
    if any(not np.isfinite(v) or v > limits[k] for k, v in measurements.items()):
        return None
    return measurements


def discover(plan, base, limits, previous=None):
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
                -sum(e.get('status') == 'submitted' for e in entries), name)
    canonical = sorted(set(result['aliases'].values()), key=rank)
    for name in sorted(records, key=rank):
        if name in result['aliases'] or name in canonical:
            continue
        match = next(((other, evidence) for other in canonical
                      if (evidence := equivalent(plan, name, other, limits)) is not None), None)
        if match:
            other, evidence = match
            result['aliases'][name] = other
            result['evidence'][name] = dict(canonical=other, measurements=evidence,
                                           limits=limits, atom_mapping='same_order_proper_rotation')
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
    result=discover(plan,base,options(value),previous)
    for alias in set(result['aliases'])-set((previous or {}).get('aliases',{})):
        if directory(base,alias).exists():
            raise ValueError('TS graph already exists; use --deduplicate-ts to retire it with scheduler verification.')
    atomic_json(path,result)
    return result
