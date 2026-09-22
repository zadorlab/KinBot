"""Local harmonic midpoint estimates and an explicit optical approximation.

The midpoint is a Cartesian construction, not a saddle or an inversion path.
Only the original Hessian is used. No QC, relaxation, or overlap integral is
computed. The counting rule is a heuristic, not proof of coverage.
"""
from numbers import Real
import copy
import networkx as nx
import numpy as np
from ase import Atoms, units

from kinbot import constants, frequencies, geometry
from kinbot.calculation import array_fingerprint
from kinbot.molecular_symmetry import _physical_role_graphs, _preserves_resonance


METHOD = 'torsion-aligned-harmonic-midpoint-v1'
OPTICAL_MIDPOINT_CUTOFF_KCAL_MOL = 4.


def apply_midpoint_heuristic(counting, diagnostic):
    """Resolve only geometric uncertainty under the agreed 4 kcal/mol rule.

    Fixed populations, measured coverage, conflicting calculations and explicit
    assumptions are authoritative. Missing or incomplete estimates keep the
    existing warned weight-one fallback. No physical mirror count is inferred.
    """
    result = dict(counting, harmonic_midpoint=diagnostic)
    if (counting.get('status') != 'unresolved'
            or counting.get('fallback') != 'unresolved_symmetry'
            or counting.get('reason') not in (
                'Rigid mirror comparison is numerically undetermined.',
                'The represented motions do not establish mirror coverage or retained handedness.')
            or counting.get('population_scope', {}).get('mirror_allowed') is not True
            or diagnostic.get('status') != 'complete'
            or diagnostic.get('mapping_search_complete') is not True):
        return result
    energy = diagnostic.get('stable_midpoint_energy_kcal_mol')
    if (not isinstance(energy, Real) or isinstance(energy, (bool, np.bool_))
            or not np.isfinite(energy) or energy < 0):
        return result
    cutoff = OPTICAL_MIDPOINT_CUTOFF_KCAL_MOL
    weight = 1. if energy <= cutoff else 2.
    result.pop('fallback')
    result.update(status='assumed', remaining_multiplier=weight,
        heuristic='harmonic_midpoint',
        reason=(f'Harmonic midpoint approximation: {float(energy)} kcal/mol '
                f'{"<=" if weight == 1 else ">"} {cutoff:g} kcal/mol; optical factor {weight:g}. '
                'This is a counting heuristic, not an inversion barrier or measured mirror coverage.'),
        harmonic_midpoint=dict(diagnostic, cutoff_kcal_mol=cutoff,
            optical_weight_changed=weight != counting.get('remaining_multiplier'),
            used_for_optical_counting=True))
    return result


def harmonic_midpoint_energy(species, aligned_mirror, hessian, *, hessian_unit):
    """Evaluate the original quadratic model at half the supplied displacement.

    Coordinates are angstrom. Explicit units distinguish native Q-Chem's
    mass-weighted Hessian from Cartesian Gaussian/ASE Hessians. Rigid motions
    are removed consistently with the existing linearity assessment. Negative
    vibrational modes are reported separately, never subtracted from the
    stable-mode energy. Internal rotor tangents are not projected a second time.
    """
    original = np.asarray(species.geom, float)
    other = np.asarray(aligned_mirror, float)
    n = len(species.atom)
    hessian = np.asarray(hessian, float)
    if (original.shape != (n, 3) or other.shape != original.shape
            or hessian.shape != (3*n, 3*n)
            or not all(np.all(np.isfinite(a)) for a in (original, other, hessian))):
        raise ValueError('Midpoint energy requires finite matching geometry and full Hessian arrays.')
    if any(getattr(species, 'isotopes', [])):
        raise ValueError('Isotope-specific vibrational masses are not supplied by this diagnostic.')
    masses = np.repeat([constants.exact_mass[str(atom)] for atom in species.atom], 3)
    if hessian_unit == 'eV / angstrom^2':
        hessian = hessian * units.Bohr**2 / units.Hartree
    elif hessian_unit not in ('hartree / bohr^2', 'hartree / (bohr^2 * amu)'):
        raise ValueError('Unknown Hessian units for the harmonic midpoint diagnostic.')
    weighted = hessian_unit == 'hartree / (bohr^2 * amu)'
    hmw = hessian if weighted else hessian / np.sqrt(np.outer(masses, masses))
    hmw = (hmw + hmw.T) / 2.
    centered = original - geometry.get_center_of_mass(original, species.atom)
    verdict = frequencies.assess_linearity(species, hmw, centered)
    rigid_geometry = frequencies.linearize(centered, species.atom) if verdict['linear'] else centered
    translations, rotations = frequencies.rigid_body_vectors(rigid_geometry, species.atom)
    largest = max(np.linalg.norm(r) for r in rotations)
    external = list(translations) + [r / np.linalg.norm(r) for r in rotations
                                     if np.linalg.norm(r) > 1.e-8 * largest]
    _, singular, vectors = np.linalg.svd(external, full_matrices=True)
    rank = int(np.sum(singular > 1.e-10 * singular[0]))
    vibrational = vectors[rank:]
    values, modes = np.linalg.eigh(vibrational @ hmw @ vibrational.T)
    # d is the displacement to the midpoint, not to the aligned mirror.
    d = (other-original).ravel() * np.sqrt(masses) / (2. * units.Bohr)
    q = modes.T @ vibrational @ d
    terms = .5 * values * q**2 * constants.AUtoKCAL
    negative = values < 0.
    norm = float(q @ q)
    return dict(
        stable_midpoint_energy_kcal_mol=float(terms[~negative].sum()),
        signed_vibrational_energy_kcal_mol=float(terms.sum()),
        negative_mode_energy_kcal_mol=float(terms[negative].sum()),
        raw_cartesian_energy_kcal_mol=float(.5 * d @ hmw @ d * constants.AUtoKCAL),
        negative_mode_displacement_fraction=float(q[negative] @ q[negative] / norm) if norm else 0.,
        negative_mode_displacement_fraction_definition='fraction of squared mass-weighted vibrational displacement',
        negative_mode_count=int(negative.sum()), external_rank=rank,
        mode_contributions=[dict(frequency_cm1=frequencies.convert_to_wavenumbers(value),
                                 midpoint_energy_kcal_mol=float(term))
                            for value, term in zip(values, terms)],
        midpoint_geometry_angstrom=((original+other)/2.).tolist(),
        hessian_input_unit=hessian_unit,
        vibrational_masses='KinBot exact elemental masses, amu',
        internal_rotor_tangents_projected=False,
        transition_surface_preserved='not established')


def mirror_midpoint_diagnostic(species, hessian, *, hessian_unit, rotors=(), max_mappings=4096,
                              target_geometry=None):
    """Set represented mirror torsions to reference angles, then fit properly.

    Select by whole-structure RMSD, not by the resulting energy. Test chemical
    mappings preserving reacting roles, resonance and represented rotor domains.
    A capped search reports only the best tested construction, never a global
    optimum. No torsional angles are varied by a numerical optimizer.
    """
    from kinbot.optical import _Parts
    if not isinstance(max_mappings, int) or isinstance(max_mappings, bool) or max_mappings < 1:
        raise ValueError('max_mappings must be a positive integer.')
    parts = _Parts(species, rotors)
    graph = parts.graph.copy()
    for a, b, edge in graph.edges(data=True):
        edge['represented_axis'] = parts.axes.get(frozenset((a, b)))
    original = np.asarray(species.geom, float)
    target_geometry = original if target_geometry is None else np.asarray(target_geometry, float)
    if target_geometry.shape != original.shape or not np.all(np.isfinite(target_geometry)):
        raise ValueError('The reflected comparison structure has invalid coordinates.')
    reference = Atoms(species.atom, positions=original)
    instructions = []
    for rotor in rotors:
        dihed = list(rotor['dihedral'])
        if len(dihed) != 4 or len(set(dihed)) != 4 or set(dihed[1:3]) != set(rotor['axis']):
            raise ValueError('A rotor must supply four distinct atoms and its actual central axis.')
        cut = graph.copy()
        cut.remove_edge(*dihed[1:3])
        group = sorted(nx.node_connected_component(cut, dihed[2]))
        if dihed[0] in group or dihed[3] not in group:
            raise ValueError('The rotor outer references do not lie on opposite sides of its axis.')
        instructions.append((dihed, group, reference.get_dihedral(*dihed)))
    best, tested, seen, complete = None, 0, set(), True
    for target in _physical_role_graphs(graph, species):
        matcher = nx.algorithms.isomorphism.GraphMatcher(graph, target,
            node_match=lambda a, b: a['label'] == b['label'],
            edge_match=lambda a, b: (a['label'], a['represented_axis']) ==
                                     (b['label'], b['represented_axis']))
        for mapping in matcher.isomorphisms_iter():
            order = tuple(mapping[i] for i in range(len(original)))
            if order in seen:
                continue
            if tested == max_mappings:
                complete = False
                break
            seen.add(order)
            tested += 1
            if not _preserves_resonance(species, list(order)):
                continue
            mirror = Atoms(species.atom, positions=(target_geometry * [-1., 1., 1.])[list(order)])
            for dihed, group, angle in instructions:
                mirror.set_dihedral(*dihed, angle, indices=group)
            if any(abs((mirror.get_dihedral(*dihed)-angle+180.) % 360.-180.) > 1.e-6
                   for dihed, _, angle in instructions):
                raise ValueError('The constructed mirror does not retain all target dihedral angles.')
            left = mirror.positions - mirror.positions.mean(axis=0)
            right = original - original.mean(axis=0)
            u, _, vt = np.linalg.svd(left.T @ right)
            sign = np.eye(3)
            sign[-1, -1] = np.linalg.det(u @ vt)
            aligned = left @ u @ sign @ vt + original.mean(axis=0)
            rmsd = float(np.sqrt(np.mean(np.sum((aligned-original)**2, axis=1))))
            if best is None or rmsd < best[0]:
                best = rmsd, aligned, order
        if not complete:
            break
    if best is None:
        raise ValueError('No compatible chemical mapping was found for the mirror construction.')
    result = harmonic_midpoint_energy(species, best[1], hessian, hessian_unit=hessian_unit)
    midpoint = np.asarray(result['midpoint_geometry_angstrom'])
    bonds = []
    for a, b in graph.edges:
        initial = float(np.linalg.norm(original[a]-original[b]))
        distance = float(np.linalg.norm(midpoint[a]-midpoint[b]))
        if initial <= 0:
            raise ValueError('Coincident bonded atoms cannot define a midpoint distortion.')
        bonds.append(dict(atom_indices=[a, b], reference_angstrom=initial,
                          midpoint_angstrom=distance, fractional_change=distance/initial-1.))
    result.update(method=METHOD, status='complete' if complete else 'limited',
        purpose='diagnostic only; not an inversion barrier or an overlap integral',
        cutoff_kcal_mol=None, optical_weight_changed=False, additional_qc_calculations=0,
        mapping_search_complete=complete, mappings_tested=tested,
        atom_mapping=list(best[2]), atom_index_base=0,
        global_rmsd_angstrom=best[0],
        maximum_displacement_angstrom=float(np.max(np.linalg.norm(best[1]-original, axis=1))),
        aligned_mirror_geometry_angstrom=best[1].tolist(),
        represented_dihedrals=[entry[0] for entry in instructions],
        target_angles_degrees=[entry[2] for entry in instructions],
        midpoint_bond_distortions=bonds,
        maximum_relative_bond_length_change=max((abs(b['fractional_change']) for b in bonds), default=0.))
    return result


def selected_midpoint_diagnostic(species, rotors, *, target_geometry=None):
    """Use an already documented selected Hessian; never read jobs or launch QC."""
    unavailable = dict(method=METHOD, status='unavailable', cutoff_kcal_mol=None,
                       optical_weight_changed=False, additional_qc_calculations=0)
    reference = (getattr(species, 'rotor_projection', None) or {}).get('reference')
    hessian = getattr(species, 'hess', [])
    if not reference or not len(hessian):
        return dict(unavailable, reason='No selected Hessian with documented geometry and units is available.')
    if (reference.get('geometry_sha256') != array_fingerprint(species.geom)
            or reference.get('hessian_sha256') != array_fingerprint(hessian)
            or reference.get('atoms') != list(map(str, species.atom))
            or any(getattr(species, key, None) is None
                   or reference.get(key) != getattr(species, key) for key in ('source_job', 'source_row_id'))):
        return dict(unavailable, reason='The documented Hessian does not match the selected calculation.')
    expected_unit = ('hartree / (bohr^2 * amu)' if reference.get('hessian_massweighted') is True
                     else 'hartree / bohr^2' if reference.get('hessian_massweighted') is False else None)
    if expected_unit is None or reference.get('hessian_unit') != expected_unit:
        return dict(unavailable, reason='The Hessian units or mass-weighting convention are unknown.')
    try:
        result = mirror_midpoint_diagnostic(species, hessian, hessian_unit=expected_unit,
                                           rotors=rotors, target_geometry=target_geometry)
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError, np.linalg.LinAlgError) as error:
        return dict(unavailable, reason=str(error))
    result['calculation_reference'] = {key: reference.get(key) for key in (
        'source_job', 'source_row_id', 'geometry_sha256', 'hessian_sha256', 'hessian_source_job')}
    return result


def conformer_midpoint_diagnostic(species, record, *, target_geometry=None):
    """Use the conformer's own geometry and Hessian, including for explicit pairs."""
    view = copy.copy(species)
    view.geom = np.asarray(record.geometry)
    view.hess = record.hessian or ()
    view.source_job = record.source_job
    view.source_row_id = (record.hessian_reference or {}).get('source_row_id')
    view.rotor_projection = dict(reference=record.hessian_reference)
    result = selected_midpoint_diagnostic(view, (), target_geometry=target_geometry)
    if record.hessian is None:
        result['reason'] = 'This conformer record has no own Hessian; frequencies alone are insufficient.'
    return result


def conformer_pair_midpoint(species, left, right):
    """Resolve uncertain explicit mirror matching with both local quadratic models."""
    forward = conformer_midpoint_diagnostic(species, left, target_geometry=right.geometry)
    reverse = conformer_midpoint_diagnostic(species, right, target_geometry=left.geometry)
    result = dict(status='undetermined', heuristic='harmonic_midpoint_pair',
                  cutoff_kcal_mol=OPTICAL_MIDPOINT_CUTOFF_KCAL_MOL,
                  conformer_ids=[left.member_id, right.member_id],
                  forward=forward, reverse=reverse)
    if not all(d.get('status') == 'complete' and d.get('mapping_search_complete') is True
               for d in (forward, reverse)):
        result['reason'] = 'A complete own-Hessian estimate is unavailable for one or both conformers.'
        return result
    energies = [d['stable_midpoint_energy_kcal_mol'] for d in (forward, reverse)]
    if all(e <= OPTICAL_MIDPOINT_CUTOFF_KCAL_MOL for e in energies):
        result.update(status='match', reason='Both reflected-pair midpoint estimates are below the cutoff.')
    elif all(e > OPTICAL_MIDPOINT_CUTOFF_KCAL_MOL for e in energies):
        result.update(status='distinct', reason='Both reflected-pair midpoint estimates exceed the cutoff.')
    else:
        result['reason'] = 'The two conformer Hessians disagree across the midpoint cutoff.'
    return result
