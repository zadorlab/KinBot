"""Compare candidate VRC interaction potentials using a common energy zero.

This is a model-development aid, not a chemical active-space oracle. Inputs
contain matched physical state labels, explicit reference cuts and evidence.
A radial correction cannot repair the angular residual computed here.
"""
import math
import json
from pathlib import Path

from kinbot.constants import AUtoKCAL
from kinbot.execution import atomic_json, file_sha256


def _validate_observation(point, row):
    distance = point['distance_angstrom']
    if not math.isfinite(distance) or distance <= 0:
        raise ValueError('Positive finite point separations are required')
    geometry = point['geometry_angstrom']
    if (not geometry or any(len(atom) != 4 or not isinstance(atom[0], str)
            or not all(math.isfinite(v) for v in atom[1:]) for atom in geometry)):
        raise ValueError('Finite labelled Cartesian geometries are required')
    energies = row['energies_hartree']
    if not energies or any(not math.isfinite(v) for v in energies.values()):
        raise ValueError('Finite nonempty state energies are required')
    if row['geometry_angstrom'] != geometry:
        raise ValueError('Potential comparison geometries differ')


def assess_path(points, observations, *, path, asymptote, return_pairs,
                return_tolerance_kcal_mol, spin_pairs=()):
    """Measure radial changes and hysteresis without inventing state identities.

    Paths, return pairs and target/reference spin labels are supplied as data.
    A secant is diagnostic, not proof of smoothness between sampled points.
    Missing points never cause a secant across an uncomputed interval.
    """
    if (not math.isfinite(return_tolerance_kcal_mol)
            or return_tolerance_kcal_mol <= 0):
        raise ValueError('A positive finite return tolerance is required')
    if len(path) < 2 or len(set(path)) != len(path):
        raise ValueError('A radial path needs distinct point identifiers')
    if any(len(pair) != 2 or pair[0] == pair[1] for pair in return_pairs):
        raise ValueError('Return pairs require distinct forward/return identifiers')
    required = list(dict.fromkeys([asymptote, *path,
                                   *(p for pair in return_pairs for p in pair)]))
    if any(p not in points for p in required):
        raise ValueError('Path and return points must be declared')
    labels = [pair['label'] for pair in spin_pairs]
    if len(labels) != len(set(labels)):
        raise ValueError('Spin-gap labels must be unique')
    for pair in return_pairs:
        if points[pair[0]]['geometry_angstrom'] != points[pair[1]]['geometry_angstrom']:
            raise ValueError('Return checks require exactly identical nuclear coordinates')
    available = {p: observations[p] for p in required if p in observations}
    states = None
    for p, row in available.items():
        _validate_observation(points[p], row)
        keys = set(row['energies_hartree'])
        if states is not None and keys != states:
            raise ValueError('Path state manifolds differ; supply a physical state mapping')
        states = keys
        if any(pair['target'] not in states or pair['reference'] not in states
               or pair['target'] == pair['reference'] for pair in spin_pairs):
            raise ValueError('Spin gaps require distinct declared target/reference states')
    result = {'missing': [p for p in required if p not in available],
              'relative_energies_kcal_mol': {}, 'spin_gaps': {},
              'radial_secants': [], 'return_checks': [],
              'production_approved': False, 'state_tracking_verified': False,
              'interpretation': 'Discrete potential and return diagnostics; no interpolation or automatic chemical selection.'}
    if asymptote in available:
        far = available[asymptote]['energies_hartree']
        zero = min(far.values())
        for p, row in available.items():
            energies = row['energies_hartree']
            result['relative_energies_kcal_mol'][p] = {
                s: (e-zero)*AUtoKCAL for s, e in energies.items()}
            result['spin_gaps'][p] = {pair['label']: {
                'gap_kcal_mol': (energies[pair['target']]-energies[pair['reference']])*AUtoKCAL,
                'change_from_far_kcal_mol': ((energies[pair['target']]-energies[pair['reference']])
                    -(far[pair['target']]-far[pair['reference']]))*AUtoKCAL}
                for pair in spin_pairs}
    for first, second in zip(path, path[1:]):
        dr = points[second]['distance_angstrom']-points[first]['distance_angstrom']
        if not math.isfinite(dr) or abs(dr) < 1.e-8:
            raise ValueError('Radial secants require distinct finite separations')
        if first not in available or second not in available:
            continue
        a, b = (available[p]['energies_hartree'] for p in (first, second))
        result['radial_secants'].append({'from': first, 'to': second,
            'kcal_mol_per_angstrom': {s: (b[s]-e)*AUtoKCAL/dr for s, e in a.items()}})
    for first, returned in return_pairs:
        if first not in available or returned not in available:
            continue
        a, b = (available[p]['energies_hartree'] for p in (first, returned))
        delta = {s: (b[s]-e)*AUtoKCAL for s, e in a.items()}
        maximum = max(map(abs, delta.values()))
        result['return_checks'].append({'forward': first, 'returned': returned,
            'returned_minus_forward_kcal_mol': delta,
            'maximum_absolute_difference_kcal_mol': maximum,
            'within_requested_tolerance': maximum <= return_tolerance_kcal_mol})
    result['status'] = ('incomplete' if result['missing'] else
        'return_test_failed' if any(not r['within_requested_tolerance'] for r in result['return_checks']) else
        'return_test_passed' if result['return_checks'] else 'path_observed')
    return result


def compare_potentials(points, candidate, reference, *, asymptote,
                       angular_reference, tolerance_kcal_mol):
    """Return raw errors and the residual after a radial electronic correction.

``points`` maps IDs to geometry and separation. Both models map point IDs to
``energies_hartree`` dictionaries keyed by physical state labels. Each model
uses ONE zero: the lowest included asymptotic energy. Missing observations are
reported, never replaced by a smaller comparator or fitted-away root offsets.
``angular_reference`` maps off-axis points to a reference cut at the SAME R.
Sorted root numbers are usable only as explicitly provisional spectrum labels;
the caller must not confuse this comparison with state tracking.
"""
    if not math.isfinite(tolerance_kcal_mol) or tolerance_kcal_mol <= 0:
        raise ValueError('A positive finite comparison tolerance is required')
    if asymptote not in points:
        raise ValueError('Asymptote must be a declared point')
    for name, anchor in angular_reference.items():
        if name not in points or anchor not in points or name == anchor:
            raise ValueError('Angular comparisons need two distinct declared points')
        if abs(points[name]['distance_angstrom']-points[anchor]['distance_angstrom']) > 1.e-8:
            raise ValueError('An angular correction comparison must use the same separation')
    missing = [{'model': name, 'point': point} for name, model in
               [('candidate', candidate), ('reference', reference)]
               for point in points if point not in model]
    result = {'missing': missing, 'raw_comparisons': [], 'angular_residuals': [],
              'status': 'incomplete', 'production_approved': False,
              'zero_policy': 'one lowest included asymptotic energy per model, shared by all states',
              'state_tracking_verified': False, 'reference_adequacy_verified': False}
    for model in (candidate, reference):
        for point, row in model.items():
            if point not in points:
                continue
            _validate_observation(points[point], row)
    if asymptote not in candidate or asymptote not in reference:
        return result
    if set(candidate[asymptote]['energies_hartree']) != set(reference[asymptote]['energies_hartree']):
        raise ValueError('Different state manifolds need an explicit physical state mapping')
    zeros = [min(m[asymptote]['energies_hartree'].values()) for m in (candidate, reference)]
    raw = {}
    for point in points.keys() & candidate.keys() & reference.keys():
        lo = candidate[point]['energies_hartree']; hi = reference[point]['energies_hartree']
        if set(lo) != set(hi) or set(lo) != set(candidate[asymptote]['energies_hartree']):
            raise ValueError('Point state labels differ from the asymptotic manifold')
        values = {state: ((lo[state]-zeros[0])-(hi[state]-zeros[1]))*AUtoKCAL for state in lo}
        raw[point] = values
        result['raw_comparisons'].append({'point': point,
            'candidate_minus_reference_kcal_mol': values,
            'maximum_absolute_difference_kcal_mol': max(map(abs, values.values()))})
    for point, anchor in angular_reference.items():
        if point not in raw or anchor not in raw:
            continue
        values = {s: raw[point][s]-raw[anchor][s] for s in raw[point]}
        maximum = max(map(abs, values.values()))
        result['angular_residuals'].append({'point': point, 'reference_cut_point': anchor,
            'corrected_candidate_minus_reference_kcal_mol': values,
            'maximum_absolute_residual_kcal_mol': maximum,
            'within_requested_tolerance': maximum <= tolerance_kcal_mol})
    result['raw_comparisons'].sort(key=lambda x: x['point'])
    if not missing and angular_reference:
        result['status'] = ('angular_test_passed' if all(
            x['within_requested_tolerance'] for x in result['angular_residuals']) else 'angular_test_failed')
    result['next_required_evidence'] = [
        'Matched physical state and orbital-manifold identities along the contributing paths',
        'Justified reference: expansion or independent method, basis and intruder/shift sensitivity',
        'Forward/reverse and independent-start checks at radial and angular neighbors']
    result['subsequent_production_checks'] = [
        'KinBot relaxed/frozen scan correction and checks of its angular dependence',
        'Corrected potential and flux sensitivity over the requested kinetic conditions',
        'Achieved RotdPy Monte Carlo uncertainty and validated MESS handoff']
    result['workflow_scope'] = ('This report aids active-space selection. Production relaxation and '
        'kinetics validation are performed by the normal KinBot workflow after selection; '
        'they do not require a duplicate reconnaissance campaign.')
    return result


def assess(source, output):
    """Read backend-independent, hash-bound observations and explicit comparisons."""
    source = Path(source).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    if output.exists():
        raise FileExistsError('Preserving existing assessment: choose a fresh output')
    data = json.loads(source.read_text())
    if data.get('schema') != 1 or not (data.get('comparisons') or data.get('path_tests')):
        raise ValueError('VRC assessment requires schema 1 and explicit comparisons or path tests')
    models = data['models']
    for observations in models.values():
        for row in observations.values():
            if not row.get('evidence'):
                raise ValueError('Every observation needs hash-bound native evidence')
            for artifact in row['evidence']:
                path = (source.parent / artifact['path']).resolve()
                if not path.is_file() or file_sha256(path) != artifact['sha256']:
                    raise ValueError(f'Observation evidence changed or is missing: {path}')
    results = []
    identifiers = set()
    for comparison in data.get('comparisons', []):
        if comparison['id'] in identifiers:
            raise ValueError('Comparison identifiers must be unique')
        identifiers.add(comparison['id'])
        if comparison['candidate'] == comparison['reference']:
            raise ValueError('A candidate cannot validate itself')
        points = {p: data['points'][p] for p in comparison['points']}
        result = compare_potentials(points, models.get(comparison['candidate'], {}),
            models.get(comparison['reference'], {}), asymptote=comparison['asymptote'],
            angular_reference=comparison['angular_reference'],
            tolerance_kcal_mol=comparison['tolerance_kcal_mol'])
        results.append({'id': comparison['id'], 'candidate': comparison['candidate'],
            'reference': comparison['reference'], 'question': comparison['question'],
            'state_assignment': comparison['state_assignment'], **result})
    paths = []
    for test in data.get('path_tests', []):
        if test['id'] in identifiers:
            raise ValueError('Assessment identifiers must be unique')
        identifiers.add(test['id'])
        result = assess_path(data['points'], models.get(test['model'], {}),
            path=test['path'], asymptote=test['asymptote'], return_pairs=test['return_pairs'],
            return_tolerance_kcal_mol=test['return_tolerance_kcal_mol'],
            spin_pairs=test.get('spin_pairs', []))
        paths.append({'id': test['id'], 'model': test['model'],
            'state_assignment': test['state_assignment'], **result})
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output, {'schema': 1, 'source': str(source), 'source_sha256': file_sha256(source),
        'comparisons': results, 'path_tests': paths, 'selection': None,
        'interpretation': 'Measured potential tests aid selection; they do not establish reference adequacy, state identity or a completed kinetics run.'})
    return output
