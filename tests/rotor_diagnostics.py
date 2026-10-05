"""Read-only diagnostics for branch-tracked torsion/inversion observations.

These functions report coordinates and evidence, not a symmetry number or a
partition function. Fixed atom labels are retained across every scan point.
"""
import numpy as np
from itertools import combinations
from ase import Atoms
from ase.units import mol, kcal
from kinbot.molecular_symmetry import proper_rmsd


def signed_pyramidalization(geometry, center, neighbors):
    """Signed center-to-neighbor-plane distance in the input length unit."""
    xyz = np.asarray(geometry, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.all(np.isfinite(xyz)):
        raise ValueError('Pyramidalization requires finite Cartesian coordinates.')
    a, b, c = xyz[list(neighbors)]
    normal = np.cross(b-a, c-a)
    length = np.linalg.norm(normal)
    if length < 1.e-12:
        raise ValueError('Pyramidalization requires three non-collinear neighbors.')
    return float(np.dot(xyz[center]-a, normal)/length)


def hydroxymethyl_coordinates(geometry):
    """Coordinates for atom order C,O,H(C),H(C),H(O).

    The transverse H-H direction reverses under carbon-H exchange, so its
    torsional coordinate advances by 180 degrees while signed wag reverses.
    The ordinary H-C-O-H dihedral generally does not advance by 180 degrees.
    """
    xyz = np.asarray(geometry, dtype=float)
    if xyz.shape != (5, 3) or not np.all(np.isfinite(xyz)):
        raise ValueError('CH2OH coordinates require five finite Cartesian positions.')
    axis = xyz[1]-xyz[0]
    if np.linalg.norm(axis) < 1.e-12:
        raise ValueError('The C-O axis must have nonzero length.')
    axis /= np.linalg.norm(axis)
    hh = xyz[2]-xyz[3]
    oh = xyz[4]-xyz[1]
    hh -= np.dot(hh, axis)*axis
    oh -= np.dot(oh, axis)*axis
    if min(np.linalg.norm(hh), np.linalg.norm(oh)) < 1.e-12:
        raise ValueError('Torsion is undefined for an axial transverse vector.')
    torsion = np.degrees(np.arctan2(np.dot(axis, np.cross(hh, oh)), np.dot(hh, oh))) % 360.
    return {'hcoh_dihedral_degrees': float(Atoms('COHHH', positions=xyz).get_dihedral(2, 0, 1, 4)),
            'hh_frame_torsion_degrees': float(torsion),
            'signed_pyramidalization_angstrom': signed_pyramidalization(xyz, 0, [1, 2, 3])}


def continuation_diagnostics(points, energy_tolerance_kcal_mol=.02):
    """Assess one 30-degree, 0..360-degree continuation with eV energies.

    Missing/failed points leave periodicity unknown. Closure refers to fixed
    atom labels; a branch jump is evidence, not an optimization failure.
    """
    if not np.isfinite(energy_tolerance_kcal_mol) or energy_tolerance_kcal_mol < 0:
        raise ValueError('Energy tolerance must be finite and nonnegative.')
    by_index = {point['index']: point for point in points}
    if len(by_index) != len(points):
        raise ValueError('Pass one continuation at a time; point indices must be unique.')
    successful = {index: point for index, point in by_index.items()
                  if point.get('status') == 'normal'}
    invalid = []
    for index, point in successful.items():
        try:
            if not all(np.isfinite(float(point[key])) for key in
                       ('energy_ev', 'angle_target_degrees', 'angle_actual_degrees')):
                raise ValueError('Nonfinite energy or angle')
            hydroxymethyl_coordinates(point['geometry'])
        except (KeyError, TypeError, ValueError) as error:
            invalid.append({'index': index, 'reason': str(error)})
    complete = set(successful) == set(range(13)) and not invalid
    result = {'complete': complete, 'sampled_half_turn_periodicity': None,
              'maximum_half_turn_difference_kcal_mol': None,
              'closure_energy_difference_kcal_mol': None,
              'closure_fixed_label_rmsd_angstrom': None,
              'constraint_max_error_degrees': None,
              'branch_sign_changes': None, 'invalid_points': invalid, 'rate_model_ready': False}
    if not complete:
        return result
    ordered = [successful[index] for index in range(13)]
    angles = [point['angle_target_degrees'] for point in ordered]
    steps = np.diff(angles)
    if not (np.allclose(steps, 30.) or np.allclose(steps, -30.)):
        raise ValueError('This diagnostic requires a signed 30-degree continuation grid.')
    errors = [abs((point['angle_actual_degrees'] - point['angle_target_degrees'] + 180.)
                  % 360. - 180.) for point in ordered]
    result['constraint_max_error_degrees'] = max(errors)
    if max(errors) > .2:
        return result
    energies = np.array([point['energy_ev'] for point in ordered]) / (kcal/mol)
    maximum = float(np.max(np.abs(energies[:6]-energies[6:12])))
    signs = [np.sign(hydroxymethyl_coordinates(point['geometry'])[
        'signed_pyramidalization_angstrom']) for point in ordered]
    result.update(maximum_half_turn_difference_kcal_mol=maximum,
        sampled_half_turn_periodicity=maximum <= energy_tolerance_kcal_mol,
        closure_energy_difference_kcal_mol=float(energies[-1]-energies[0]),
        closure_fixed_label_rmsd_angstrom=proper_rmsd(ordered[0]['geometry'], ordered[-1]['geometry']),
        branch_sign_changes=int(sum(a*b < 0 for a, b in zip(signs[:-1], signs[1:]))))
    return result


def cross_continuation_diagnostics(points):
    """Compare saved continuations at matching angles, excluding closure repeats.

    RMSD uses proper rotations and fixed atom labels to expose direction/seed
    dependence; it neither counts branches nor establishes equilibration.
    """
    groups = {}
    for point in points:
        groups.setdefault((point['branch_seed'], point['direction']), []).append(point)
    checked = {key: continuation_diagnostics(group) for key, group in groups.items()}
    invalid = [list(key) for key, value in checked.items()
               if not value['complete'] or value['sampled_half_turn_periodicity'] is None]
    result = {'comparisons': [], 'seed_starts': [], 'invalid_continuations': invalid,
              'lowest_observed_envelope_half_turn_difference_kcal_mol': None,
              'independent_branch_coverage_established': False, 'rate_model_ready': False}
    if invalid:
        return result
    ordered = {key: sorted(groups[key], key=lambda p: p['index'])
               for key in sorted(groups, key=lambda k: (k[0], -k[1]))}
    for key, group in ordered.items():
        first = group[0]
        start = first.get('start_geometry')
        result['seed_starts'].append({'branch_seed': key[0], 'direction': key[1],
            'start_wag_angstrom': (hydroxymethyl_coordinates(start)['signed_pyramidalization_angstrom']
                                   if start is not None else None),
            'accepted_wag_angstrom': hydroxymethyl_coordinates(first['geometry'])['signed_pyramidalization_angstrom']})
    for left_key, right_key in combinations(ordered, 2):
        pairs = []
        for first in ordered[left_key][:12]:
            distances = [abs((p['angle_target_degrees'] - first['angle_target_degrees'] + 180.)
                             % 360. - 180.) for p in ordered[right_key][:12]]
            index = int(np.argmin(distances))
            if distances[index] > .2:
                raise ValueError('Cross-continuation grids do not share the same torsional angles.')
            second = ordered[right_key][index]
            pairs.append({'left_index': first['index'], 'right_index': second['index'],
                'angle_degrees': first['angle_target_degrees'] % 360.,
                'energy_left_minus_right_kcal_mol': (first['energy_ev']-second['energy_ev'])/(kcal/mol),
                'proper_fixed_label_rmsd_angstrom': proper_rmsd(first['geometry'], second['geometry']),
                'left_wag_angstrom': hydroxymethyl_coordinates(first['geometry'])['signed_pyramidalization_angstrom'],
                'right_wag_angstrom': hydroxymethyl_coordinates(second['geometry'])['signed_pyramidalization_angstrom']})
        result['comparisons'].append({'left': list(left_key), 'right': list(right_key), 'points': pairs})
    # This is only the lowest *observed* energy at each angle, not a claim of
    # global constrained minima or a validated one-coordinate partition sum.
    if ordered:
        reference = next(iter(ordered.values()))[:12]
        envelope = []
        for first in reference:
            matched = [p['energy_ev'] for group in ordered.values() for p in group[:12]
                       if abs((p['angle_target_degrees']-first['angle_target_degrees']+180.)%360.-180.) < .2]
            envelope.append(min(matched)/(kcal/mol))
        result['lowest_observed_envelope_half_turn_difference_kcal_mol'] = float(
            max(abs(envelope[i]-envelope[i+6]) for i in range(6)))
    return result
