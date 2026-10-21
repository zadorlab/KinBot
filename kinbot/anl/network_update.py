"""Explicit append-only expansion of a frozen automatic network."""
from copy import deepcopy
import json
from pathlib import Path

from kinbot.execution import atomic_json
from kinbot.anl.ts_reuse import digest


def extend(previous, current):
    """Permit additional work, never change the scientific input of saved work."""
    for key in set(previous) | set(current):
        if key not in ('network', 'cbh') and previous.get(key) != current.get(key):
            raise ValueError(f'Network extension changed {key}; scientific migration is required.')
    current = deepcopy(current)
    a, b = previous['network'], current['network']
    if a['initial'] != b['initial']:
        raise ValueError('Network extension changed the initial well.')

    def retain(old, new, label, *, roles=False):
        for key, value in old.items():
            x, y = deepcopy(value), deepcopy(new.get(key))
            if roles and y is not None:
                if not set(x.pop('roles', [])) <= set(y.pop('roles', [])):
                    raise ValueError(f'{label}/{key}: an existing network role disappeared.')
            if x != y:
                raise ValueError(f'{label}/{key}: saved scientific record changed or disappeared; '
                                 'existing calculations were retained for review.')
    retain(a['route_to_smiles'], b['route_to_smiles'], 'species identity')
    retain({r['routing_name']: r for r in a['stable_species']},
           {r['routing_name']: r for r in b['stable_species']}, 'stable species', roles=True)
    retain(a['transition_states'], b['transition_states'], 'transition state')
    retain({r['name']: r for r in a['reactions']},
           {r['name']: r for r in b['reactions']}, 'reaction')
    retain(a.get('product_complexes', {}), b.get('product_complexes', {}), 'product complex')
    for key in ('targets', 'terminal_atct_species'):
        retain(previous['cbh'][key], current['cbh'][key], 'CBH '+key)
    identity = lambda r: (r['smiles'], r['charge'], r['multiplicity'])
    for key in ('calculated_states', 'thermochemistry_states'):
        old = {identity(r): r for r in previous['cbh'][key]}
        for record in current['cbh'][key]:
            prior = old.pop(identity(record), None)
            if prior is None:
                continue
            before, after = deepcopy(prior), deepcopy(record)
            revision = before.pop('graph_revision', None)
            after.pop('graph_revision', None)
            # Role expansion is bookkeeping, never a changed geometry input.
            for value in (before, after):
                source = value.get('network_record')
                if isinstance(source, dict):
                    source.pop('roles', None)
            if before == after:
                if revision is not None:
                    record['graph_revision'] = revision
                continue
            # A newly discovered network state can have existed as a CBH
            # reference. Its geometry source changes; use a new graph rather
            # than silently attaching the old reference's single points.
            if (before.pop('source', None) == 'automatic_reference_l2'
                    and after.pop('source', None) == 'network_l2'):
                before.pop('network_record', None)
                source = after.pop('network_record', None)
                if before == after:
                    if key == 'calculated_states':
                        record['graph_revision'] = digest(source)[:16]
                    continue
            raise ValueError(f'CBH {key}/{identity(record)}: saved scientific input changed.')
        if old:
            raise ValueError(f'CBH {key}: saved states disappeared.')
    return current


def archive_assembly(base, archive):
    """Retain derived evidence without moving any native QC graph or cache."""
    base, archive = Path(base), Path(archive)
    archive.mkdir(parents=True, exist_ok=True)
    for name in ('handoff.json', 'execution.json', 'final_report.json',
                 'cbh_thermochemistry_models.json', 'profiles', 'formations'):
        source = base/name
        if source.exists():
            dest = archive/name
            if dest.exists():
                raise ValueError(f'{source}: changed during network extension.')
            source.rename(dest)


def publish(previous, current, base):
    """Consume only the exact revision authorized after branch quiescence."""
    base = Path(base)
    path = base/'network_revisit.json'
    if not path.exists():
        return False
    authorization = json.loads(path.read_text())
    if authorization.get('previous_sha256') != digest(previous):
        raise ValueError('Network revisit authorization does not match the frozen plan.')
    current = extend(previous, current)
    archive = Path(authorization['archive'])
    archive.mkdir(parents=True, exist_ok=True)
    saved = archive/'network_request.json'
    if saved.exists() and json.loads(saved.read_text()) != previous:
        raise ValueError('Network revision archive already contains different inputs.')
    atomic_json(saved, previous)
    # No complete kinetics may be reused after extending its network.
    archive_assembly(base, archive)
    atomic_json(base/'network_request.json', current)
    authorization.update(status='extended', current_sha256=digest(current))
    atomic_json(path, authorization)
    return current
