"""Use actual abstraction rotor groups without inventing a chemical bond."""
import copy

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from kinbot.calculation import array_fingerprint, geometry_reference
from kinbot.counting_contract import optical_counting
from kinbot.optical_harmonic import selected_midpoint_diagnostic, mirror_midpoint_diagnostic
from kinbot.reaction_path import prepare_stereopath
from kinbot.stationary_pt import StationaryPoint


def thf_abstraction():
    # Rounded puckered-ring geometry with OH approaching H11 on C4. The unit
    # Hessian tests coordinate handling; no calculated modes or energies enter.
    geometry = np.array([
        [1.47, 2.44, -.54], [2.12, 1.94, -1.07], [-.39, -.54, -1.20],
        [-.19, -.49, .32], [1.24, .01, .41], [1.93, -.52, -.70],
        [.99, -.99, -1.68], [-1.18, -1.22, -1.51], [-.62, .46, -1.57],
        [-.26, -1.49, .75], [-.89, .16, .83], [1.23, 1.15, .33],
        [1.78, -.23, 1.33], [1.07, -2.08, -1.73], [1.26, -.57, -2.65]])
    species = StationaryPoint('nonbonded_rotor', 0, 2,
        atom=['O', 'H', 'C', 'C', 'C', 'O', 'C'] + ['H'] * 8,
        geom=geometry, wellorts=1)
    bond = np.zeros((15, 15), dtype=int)
    for a, b in [(0, 1), (0, 11), (2, 3), (3, 4), (4, 5), (5, 6),
                 (6, 2), (2, 7), (2, 8), (3, 9), (3, 10), (4, 11),
                 (4, 12), (6, 13), (6, 14)]:
        bond[a, b] = bond[b, a] = 1
    species.bond, species.bond01 = bond.copy(), bond.copy()
    species.bonds = [bond.copy()]
    species.characterize(bond_mx=bond)
    endpoints = []
    for missing in [(0, 11), (4, 11)]:
        matrix = bond.copy()
        matrix[missing] = matrix[missing[::-1]] = 0
        point = StationaryPoint('endpoint', 0, 2, atom=species.atom, geom=geometry)
        point.bond, point.bond01 = matrix.copy(), matrix.copy()
        point.bonds = [matrix]
        point.characterize(bond_mx=matrix)
        endpoints.append(point)
    species.reac_bond = endpoints[1].bond01 - endpoints[0].bond01
    prepare_stereopath(species, *endpoints)
    species.dihed = [[3, 4, 0, 1]]
    species.freq = [-500.] + [100.] * (3 * species.natom - 7)
    species.reduced_freqs = species.freq[:-1]
    species.hess = np.eye(3 * species.natom)
    species.source_job, species.source_row_id = 'nonbonded_rotor', 1
    species.energy, species.zpe = 0., 0.
    species.optical_population = 'specified'
    reference = dict(geometry_reference(species), dihedrals=species.dihed,
        hessian_sha256=array_fingerprint(species.hess),
        hessian_massweighted=False, hessian_unit='hartree / bohr^2')
    species.rotor_projection = dict(reference=reference, internal_rank=1,
        rotors=[dict(rotor_index=0, projected=True)])
    rotor = dict(index=0, axis=[4, 0], dihedral=species.dihed[0],
        partition=[list(range(2, 15)), [0, 1]], sigma_int=1, usable=True,
        represented_domain_degrees=[0., 360.], points=[])
    for i in range(2):
        atoms = Atoms(species.atom, positions=geometry)
        atoms.rotate_dihedral(*species.dihed[0], 180. * i, indices=[0, 1])
        rotor['points'].append(dict(index=i, status='successful',
            angle_offset_degrees=180. * i, geometry_angstrom=atoms.positions.tolist(),
            electronic_energy_hartree=0.))
    return species, dict(status='recorded', scan_reference=reference, rotors=[rotor])


def test_nonbonded_rotor_reaches_midpoint_without_changing_the_model():
    species, hir = thf_abstraction()
    assert species.bond[0, 11] == species.bond[4, 11] == 1
    assert species.reac_bond[0, 11] == 1
    assert species.reac_bond[4, 11] == -1
    assert np.count_nonzero(species.reac_bond) == 4
    original = copy.deepcopy((species.geom, species.hess, species.dihed, hir))
    axis = hir['rotors'][0]['axis']
    assert species.bond[tuple(axis)] == 0
    result = optical_counting(species, hir)
    assert result['heuristic'] == 'harmonic_midpoint'
    midpoint = result['harmonic_midpoint']
    assert midpoint['status'] == 'complete'
    assert midpoint['additional_qc_calculations'] == 0
    assert result['remaining_multiplier'] == 2.
    aligned = Atoms(species.atom, positions=midpoint['aligned_mirror_geometry_angstrom'])
    for dihedral, angle in zip(midpoint['represented_dihedrals'], midpoint['target_angles_degrees']):
        assert (aligned.get_dihedral(*dihedral) - angle + 180) % 360 - 180 == pytest.approx(0., abs=1.e-6)
    np.testing.assert_array_equal(species.geom, original[0])
    np.testing.assert_array_equal(species.hess, original[1])
    assert species.dihed == original[2]
    # Optical evidence may annotate a scan, but never change its observations.
    assert hir['rotors'][0]['points'] == original[3]['rotors'][0]['points']


def test_nonbonded_rotor_requires_recorded_complete_groups():
    species, hir = thf_abstraction()
    for groups in (None, [[0, 1], [1, 2]], [[0], list(range(1, species.natom - 1))]):
        rotors = copy.deepcopy(hir['rotors'])
        rotors[0]['partition'] = groups
        result = selected_midpoint_diagnostic(species, rotors)
        assert result['status'] == 'unavailable'
        assert 'recorded HIR groups' in result['reason']


def test_reversing_nonbonded_axis_and_recorded_groups_preserves_midpoint():
    species, hir = thf_abstraction()
    rotors = hir['rotors']
    first = selected_midpoint_diagnostic(species, rotors)
    reverse = copy.deepcopy(rotors)
    reverse[0]['axis'].reverse()
    reverse[0]['dihedral'].reverse()
    reverse[0]['partition'].reverse()
    second = selected_midpoint_diagnostic(species, reverse)
    assert second['status'] == 'complete'
    assert second['global_rmsd_angstrom'] == pytest.approx(first['global_rmsd_angstrom'])
    assert second['stable_midpoint_energy_kcal_mol'] == pytest.approx(first['stable_midpoint_energy_kcal_mol'])


def test_measured_mirror_coverage_precedes_nonbonded_midpoint(monkeypatch):
    species, hir = thf_abstraction()
    points = hir['rotors'][0]['points']
    points[0].update(status='successful', geometry_angstrom=species.geom.tolist(),
                     electronic_energy_hartree=species.energy)
    points[1].update(status='successful', geometry_angstrom=(species.geom * [-1., 1., 1.]).tolist(),
                     electronic_energy_hartree=species.energy)
    def forbidden(*args, **kwargs):
        raise AssertionError('Measured coverage must not be replaced by a midpoint guess')
    monkeypatch.setattr('kinbot.optical_harmonic.selected_midpoint_diagnostic', forbidden)
    result = optical_counting(species, hir)
    assert result['remaining_multiplier'] == 1.
    assert 'measured_coverage' in result


def test_later_dihedral_must_not_undo_an_earlier_target(monkeypatch):
    from test_optical_evaluator import loose_ts
    from kinbot.thermochemistry import hir_evidence
    species = loose_ts()
    rotors = hir_evidence(species)['rotors']
    first, last = [tuple(r['dihedral']) for r in rotors]
    setter = Atoms.set_dihedral
    def coupled_setter(atoms, a, b, c, d, angle, **kwargs):
        setter(atoms, a, b, c, d, angle, **kwargs)
        if (a, b, c, d) == last:
            # Model an overlapping rotation which also changes the first angle.
            setter(atoms, *first, atoms.get_dihedral(*first) + 5., indices=rotors[0]['partition'][1])
    monkeypatch.setattr(Atoms, 'set_dihedral', coupled_setter)
    result = selected_midpoint_diagnostic(species, rotors)
    assert result['status'] == 'unavailable'
    assert 'does not retain all target dihedral angles' in result['reason']


def test_ordinary_rotor_in_disconnected_assembly_keeps_original_refusal():
    atoms = molecule('CH3OH') + Atoms('H2', positions=[[10., 0., 0.], [10.74, 0., 0.]])
    species = StationaryPoint('methanol_and_hydrogen', 0, 1,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    species.characterize()
    dihedral = species.dihed[0]
    rotor = dict(axis=dihedral[1:3], dihedral=dihedral,
                 usable=True, sigma_int=1, represented_domain_degrees=[0., 360.])
    with pytest.raises(ValueError, match='Disconnected assemblies'):
        mirror_midpoint_diagnostic(species, np.eye(3 * species.natom),
            hessian_unit='hartree / bohr^2', rotors=[rotor])


def test_chemical_ring_bond_does_not_use_nonbonded_rotor_exception():
    from kinbot.optical import evaluate_optical
    atoms = molecule('C6H6')
    species = StationaryPoint('benzene', 0, 1,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    species.characterize()
    rotor = dict(axis=[0, 1], dihedral=[5, 0, 1, 2],
        usable=True, sigma_int=1, represented_domain_degrees=[0., 360.],
        partition=[[0, 5, 6, 11], [1, 2, 3, 4, 7, 8, 9, 10]])
    with pytest.raises(ValueError, match='cannot be separated'):
        mirror_midpoint_diagnostic(species, np.eye(3 * species.natom),
            hessian_unit='hartree / bohr^2', rotors=[rotor])
    with pytest.raises(ValueError, match='cannot be separated'):
        evaluate_optical(species, rotors=[rotor])
