"""Evidence for distinguishing parallel saddles from repeated discoveries.

Endpoint and stereochemical classes are checked by the caller. A duplicate
also requires a proper atom-mapped geometry fit and its reactive displacement.
Missing evidence preserves the observation and is reported by the selector.
"""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import base64
import hashlib
import json
import re

import numpy as np

from kinbot import constants
from kinbot.anl.ts_reuse import aligned_permutations, MappingLimitReached
from kinbot.spatial import proper_alignment, mode_overlap


LIMITS = dict(rmsd_angstrom=.01, maximum_atom_error_angstrom=.03,
              minimum_mode_overlap=.99, energy_kcal_mol=.1,
              frequency_cm_inverse=5.)
MARKER = re.compile(r'^! kinbot_saddle_identity (\S+)$', re.M)


def evidence(species):
    """Freeze the selected geometry and a matching native reactive mode."""
    record = dict(symbols=list(map(str, species.atom)),
                  observation=getattr(species, 'name', getattr(species, 'source_job', None)),
                  positions=np.asarray(species.geom, float).tolist(),
                  charge=int(species.charge), multiplicity=int(species.mult),
                  isotopes=list(map(int, getattr(species, 'isotopes', [0]*len(species.atom)))),
                  frequencies=list(map(float, getattr(species, 'freq', ()))),
                  energy_hartree=getattr(species, 'energy', None),
                  source_job=getattr(species, 'source_job', None),
                  source_row_id=getattr(species, 'source_row_id', None))
    from kinbot.transition_symmetry import _documented_reactive_mode
    documented = _documented_reactive_mode(species)
    if documented is not None:
        mode, source = documented
        record['mode'] = np.asarray(mode).tolist()
        record['mode_source'] = source
    else:
        from kinbot.gaussian_modes import read_normal_modes
        for job in dict.fromkeys((getattr(species, 'hessian_source_job', None),
                                  record['source_job'])):
            if not job or not Path(str(job)+'.log').is_file():
                continue
            try:
                native = read_normal_modes(str(job)+'.log', species.atom, species.geom)
            except (OSError, ValueError):
                continue
            negative = [i for i, f in enumerate(native['frequencies_cm1']) if f < -1.]
            if len(negative) == 1:
                record['mode'] = native['modes'][negative[0]]
                record['mode_source'] = {k:v for k,v in native.items() if k != 'modes'}
                break
    indices = getattr(species, 'conformer_index', ())
    if len(indices) > 1:
        inventory = {r.index:r for r in getattr(species, 'conformer_inventory', ())}
        inventory.update(getattr(species, 'conformer_records', {}))
        record['members'] = []
        for offset, index in enumerate(indices):
            if index < 0:
                continue
            member = inventory.get(index)
            if member is None:
                record['members'].append({})
                continue
            reference = member.hessian_reference or {}
            point = SimpleNamespace(atom=species.atom, charge=species.charge, mult=species.mult,
                geom=species.conformer_geom[offset], freq=species.conformer_freq[offset],
                isotopes=record['isotopes'], energy=member.electronic_energy_hartree,
                source_job=member.source_job, source_row_id=reference.get('source_row_id'),
                hess=member.hessian, optical_hessian_reference=reference)
            record['members'].append(evidence(point))
    return record


def row_evidence(row, directory, charge, multiplicity):
    """Use the exact accepted database row selected for PES energy assembly."""
    point = SimpleNamespace(name=row.name, atom=row.symbols, geom=row.positions,
        charge=row.data.get('charge', charge), mult=row.data.get('mult', multiplicity),
        isotopes=row.data.get('isotopes', [0]*len(row.symbols)),
        freq=row.data.get('frequencies', ()),
        energy=row.data['energy']*constants.EVtoHARTREE,
        source_job=str(Path(directory)/row.name), source_row_id=row.id)
    from kinbot.conformer_records import hessian_record
    hessian = hessian_record(None, point.source_job, row.positions, row.symbols, row=row)
    if hessian:
        point.hess = hessian['hessian']
        point.optical_hessian_reference = hessian['hessian_reference']
    return evidence(point)


def compare(left, right):
    """Return equivalent/distinct/unresolved with all measured tolerances.

All admissible atom permutations are tried. The rotation fitted to positions
also rotates the displacement. The sign of an eigenvector is arbitrary.
Reflections are explicit caller operations for a demonstrated mirror orbit.
"""
    result = dict(status='unresolved', limits=dict(LIMITS))
    if not left or not right:
        return dict(result, reason='missing_saddle_evidence')
    if any(left.get(k) != right.get(k) for k in ('charge', 'multiplicity')):
        return dict(result, status='distinct', reason='different_electronic_state')
    try:
        x, y = (np.asarray(r['positions'], float) for r in (left, right))
        labels = [[(s, i) for s, i in zip(r['symbols'], r['isotopes'])] for r in (left, right)]
        if (x.shape != (len(labels[0]), 3) or y.shape != (len(labels[1]), 3)
                or not np.isfinite(x).all() or not np.isfinite(y).all()):
            return dict(result, reason='invalid_geometry')
        mappings = aligned_permutations(labels[0], x, labels[1], y, LIMITS['rmsd_angstrom'], report_limit=True)
        found = False
        for rmsd, order in mappings:
            fit = proper_alignment(y[order], x)
            if max(fit.errors) > LIMITS['maximum_atom_error_angstrom']:
                continue
            found = True
            if any(any(r['isotopes']) for r in (left, right)):
                return dict(result, reason='isotope_specific_mode_masses_unavailable')
            fx, fy = (np.asarray(r['frequencies'], float) for r in (left, right))
            if (not len(fx) or fx.shape != fy.shape or not np.isfinite(fx).all()
                    or not np.isfinite(fy).all() or any(np.count_nonzero(f < -1.) != 1 for f in (fx, fy))):
                return dict(result, reason='missing_first_order_spectrum')
            if not left.get('mode') or not right.get('mode'):
                return dict(result, reason='missing_reactive_mode')
            overlap = mode_overlap(left['mode'], np.asarray(right['mode'])[order] @ fit.rotation,
                                  masses=[constants.exact_mass[s] for s in left['symbols']])
            delta = abs(left['energy_hartree']-right['energy_hartree'])*constants.AUtoKCAL
            frequency = float(max(abs(np.sort(fx)-np.sort(fy))))
            measurements = dict(rmsd_angstrom=rmsd, atom_mapping=order,
                maximum_atom_error_angstrom=float(max(fit.errors)), mode_overlap=overlap,
                energy_kcal_mol=delta, frequency_cm_inverse=frequency)
            result['measurements'] = measurements
            if (np.isfinite(delta) and delta <= LIMITS['energy_kcal_mol']
                    and frequency <= LIMITS['frequency_cm_inverse']
                    and overlap >= LIMITS['minimum_mode_overlap']):
                if 'members' in left or 'members' in right:
                    import networkx as nx
                    a, b = left.get('members', []), right.get('members', [])
                    if not a or len(a) != len(b):
                        return dict(result, reason='different_or_unverified_conformer_ensembles')
                    graph = nx.Graph()
                    graph.add_nodes_from(range(len(a)))
                    graph.add_nodes_from(range(len(a), len(a)+len(b)))
                    graph.add_edges_from((i, len(a)+j) for i, x in enumerate(a) for j, y in enumerate(b)
                                         if compare(x,y)['status'] == 'equivalent')
                    matching = nx.algorithms.bipartite.maximum_matching(graph, top_nodes=set(range(len(a))))
                    if len(matching) != 2*len(a):
                        return dict(result, reason='different_or_unverified_conformer_ensembles')
                    result['member_mapping'] = [matching[i]-len(a) for i in range(len(a))]
                return dict(result, status='equivalent', reason='mapped_geometry_and_reactive_mode')
        return dict(result, status='distinct', reason='different_geometry_or_mode' if found
                    else 'geometry_outside_tolerance_or_mapping_limit')
    except MappingLimitReached:
        return dict(result, reason='atom_mapping_limit')
    except (KeyError, TypeError, ValueError, IndexError):
        return dict(result, reason='incomplete_saddle_evidence')


def annotate(block, record):
    """Carry selected saddle evidence through cached and mirror MESS assembly."""
    block = re.sub(MARKER.pattern+'\n?', '', block, flags=re.M)
    record = dict(record)
    record.setdefault('origin', hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest())
    payload = base64.urlsafe_b64encode(json.dumps(record, sort_keys=True).encode()).decode()
    header = re.search(r'^\s*Barrier[^\n]*\n', block, re.M)
    if header is None:
        raise ValueError('Saddle evidence requires one Barrier header.')
    return block[:header.end()]+'! kinbot_saddle_identity '+payload+'\n'+block[header.end():]


def block_evidence(block):
    match = MARKER.search(block)
    return json.loads(base64.urlsafe_b64decode(match[1])) if match else None


def reflect(block):
    record = block_evidence(block)
    if record is None:
        return block
    return annotate(block, _reflected_record(record))


def _reflected_record(record):
    record = deepcopy(record)
    for key in ('positions', 'mode'):
        if record.get(key):
            record[key] = (np.asarray(record[key])*[-1., 1., 1.]).tolist()
    if 'members' in record:
        record['members'] = [_reflected_record(member) for member in record['members']]
    return record


def same_blocks(left, right, *, mirror=False):
    """Compare already classified MESS routes, optionally in a mirror orbit."""
    a, b = block_evidence(left), block_evidence(right)
    if a and b and a.get('origin') and a['origin'] == b.get('origin'):
        # Explicitly reflected/copied versions of one native observation.
        # Independent calculations have separate source-bound origins.
        return True
    if mirror and b is not None:
        b = _reflected_record(b)
    return compare(a, b)['status'] == 'equivalent'
