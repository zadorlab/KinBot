"""Axial identities and downstream counting, without QC or calculated rates."""
import copy

import networkx as nx
import numpy as np
import pytest
from ase import Atoms

from tests import test_irc_stereo_identity as irc_checks
from tests.test_irc_stereo_identity import point
from kinbot.conformer_records import ConformerRecord
from kinbot.conformer_counting import evaluate_members
from kinbot.optical import evaluate_optical
from kinbot.reaction_path import prepare_ts_context, path_geometry_allowed
from kinbot.species_routing import routing_name, same_species
from kinbot.stereo_identity import canonical_identity, configured_geometry_allowed
from kinbot.stereochemistry import (stereotopic_relation, configuration_erased_graph,
                                    virtually_labelled)


def turn(p, center=2, end=3):
    """Reverse one axis while keeping the downstream fragment internally fixed."""
    q = copy.copy(p)
    q.geom = p.geom.copy()
    graph = nx.from_numpy_array(p.bond01)
    graph.remove_edge(center, end)
    group = list(nx.node_connected_component(graph, end))
    axis = p.geom[end] - p.geom[center]
    axis /= np.linalg.norm(axis)
    vectors = p.geom[group] - p.geom[center]
    q.geom[group] = 2*np.outer(vectors @ axis, axis)-vectors+p.geom[center]
    q.__dict__.pop('optical_reference', None)
    return q


@pytest.mark.parametrize('smiles', [
    'CC=C=CC', 'FC=C=CF', 'CCC=C=CC', 'CC=C=C[C@H](F)Cl',
    'CC=C=CC=C=CC', 'CC=C=C=C=CC', 'CC=C=C=CC',
    'Fc1ccccc1-c1ccccc1Cl', 'C1=C=CCCCCCC1', '[CH2]C(=C=CF)C',
    'CC=C=C(C=C=CC)C=C=CC'])
def test_axis_identity_survives_atom_order_rotation_and_reflection(smiles):
    p = point(smiles)
    if smiles.startswith('[CH2]'):
        p.mult = 2
    identity = canonical_identity(p)
    assert identity['status'] == 'assigned'
    mirror = canonical_identity(p, p.geom * [-1., 1., 1.])
    assert identity['mirror_id'] == mirror['id']
    assert identity['mirror_family_id'] == mirror['mirror_family_id']
    rng = np.random.default_rng(463)
    original = [b.copy() for b in p.bonds]
    for _ in range(8):
        order = rng.permutation(p.natom)
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        rotation[:, 0] *= np.linalg.det(rotation)
        q = copy.copy(p)
        q.atom = np.asarray(p.atom)[order]
        q.geom = p.geom[order] @ rotation + [2., -3., 4.]
        q.bonds = [b[np.ix_(order, order)] for b in reversed(p.bonds)]
        assert canonical_identity(q)['id'] == identity['id']
    for before, after in zip(original, p.bonds):
        np.testing.assert_array_equal(before, after)


def test_relative_axes_isotopes_and_reaction_sites_remain_distinct():
    # Two equivalent axes: RR/SS and meso, rather than one shared species.
    p = point('CC=C=CC=C=CC')
    variants = [p, turn(p), turn(p, 5, 6), turn(turn(p), 5, 6)]
    identities = [canonical_identity(v) for v in variants]
    assert len({i['id'] for i in identities}) == 3
    assert len({i['mirror_family_id'] for i in identities}) == 2
    # Side axes can make the central axis stereogenic only after their own
    # handedness is assigned. Retain the resulting diastereomers.
    p = turn(point('CC=C=C(C=C=CC)C=C=CC'), 5, 6)
    q = turn(p)
    assert not same_species(p, q)
    assert not configured_geometry_allowed(p, q.geom, 'racemic')
    assert canonical_identity(p)['mirror_id'] == canonical_identity(
        p, p.geom * [-1., 1., 1.])['id']
    # Changing the axis alone next to a fixed tetrahedral center gives a
    # diastereomer, not another conformer of that center's racemate.
    p = point('CC=C=C[C@H](F)Cl')
    q = turn(p)
    assert not same_species(p, q)
    assert not configured_geometry_allowed(p, q.geom, 'racemic')
    assert configuration_erased_graph(canonical_identity(p)) == configuration_erased_graph(canonical_identity(q))
    p = point('CCC=C=CC')
    hydrogens = [i for i, atom in enumerate(p.atom) if atom == 'H' and p.bond[1, i]]
    assert stereotopic_relation(p, *hydrogens) == 'diastereotopic'
    # One CH2 terminus is not an axis until isotope substitution separates
    # its two H atoms. Virtual labels stay confined to the identity view.
    p = point('C=C=CC')
    hydrogen = next(i for i, atom in enumerate(p.atom) if atom == 'H' and p.bond[0, i])
    assert not canonical_identity(p)['is_chiral_configuration']
    assert canonical_identity(virtually_labelled(p, {hydrogen: 0}))['is_chiral_configuration']
    q = copy.copy(p)
    q.isotopes = [0] * p.natom
    q.isotopes[hydrogen] = 2
    assert canonical_identity(q)['is_chiral_configuration']
    assert not hasattr(p, 'isotopes')


def test_axial_irc_ts_and_mc_use_the_requested_enantiomers():
    p = point('CC=C=CC')
    from unittest.mock import Mock
    from kinbot.stereo_identity import log_input_stereochemistry
    logger = Mock()
    log_input_stereochemistry(p, {'smiles': 'CC=C=CC'}, logger)
    assert logger.info.call_count == 1
    logger.warning.assert_not_called()
    mirror = copy.copy(p)
    mirror.geom = p.geom * [-1., 1., 1.]
    assert routing_name(p) != routing_name(mirror)
    assert irc_checks.TestIRCStereoIdentity().interpret(p, mirror) != 0
    assert irc_checks.TestIRCStereoIdentity().interpret(p, mirror, 'racemic') == 0
    assert evaluate_optical(p, population='specified')['remaining_multiplier'] == 1.
    assert evaluate_optical(p, population='racemic')['remaining_multiplier'] == 2.
    freq = [1000.] * (3*p.natom-6)  # example properties, not calculated spectra
    records = [ConformerRecord(str(i), i, str(i), 'valid', geometry=g,
                zero_energy_hartree=-100., frequencies_cm1=freq)
               for i, g in enumerate((p.geom, mirror.geom))]
    counted, groups = evaluate_members(p, records, population='specified')
    assert groups == [[0]]
    assert counted[1].exclusion_reason == 'different configured stereoisomer'
    counted, groups = evaluate_members(p, records, population='racemic')
    assert len(groups) == 1
    assert [r.remaining_optical_weight for r in counted] == [1., 1.]
    # A structural TS example preserves a nonreacting axis. This is not an
    # optimized saddle or proof that this abstraction has a barrier.
    q = copy.copy(p)
    q.bond = p.bond.copy()
    h = next(i for i, atom in enumerate(p.atom) if atom == 'H' and p.bond[0, i])
    q.bond[0, h] = q.bond[h, 0] = 0
    q.bonds = [q.bond.copy()]
    ts = copy.copy(p)
    ts.wellorts = 1
    prepare_ts_context(ts, p, q)
    assert not path_geometry_allowed(ts, ts.geom * [-1., 1., 1.], 'specified')
    assert path_geometry_allowed(ts, ts.geom * [-1., 1., 1.], 'racemic')


def test_biaryl_hir_respects_specified_and_racemic_models():
    p = point('Fc1ccccc1-c1ccccc1Cl')
    axis = next((i, j) for i, j in nx.bridges(nx.from_numpy_array(p.bond01))
                if p.cycle[i] and p.cycle[j])
    outer = [next(i for i in np.flatnonzero(p.bond[a]) if i != b) for a, b in (axis, axis[::-1])]
    dihedral = [outer[0], *axis, outer[1]]
    mirror = p.geom * [-1., 1., 1.]
    delta = (Atoms(p.atom, positions=mirror).get_dihedral(*dihedral)
             - Atoms(p.atom, positions=p.geom).get_dihedral(*dihedral)) % 360.
    rotor = dict(index=0, axis=axis, dihedral=dihedral, usable=True, sigma_int=1,
        represented_domain_degrees=[0., 360.], points=[dict(index=i, status='successful',
        geometry_angstrom=g, electronic_energy_hartree=0., angle_offset_degrees=angle)
        for i, (g, angle) in enumerate(((p.geom, 0.), (mirror, delta)))])
    result = evaluate_optical(p, rotors=[rotor], population='specified')
    assert result['invalid_rotor_index'] == 0
    assert evaluate_optical(p, rotors=[rotor], population='racemic')['remaining_multiplier'] == 1.
    # A real scan also crosses a planar geometry and usually misses the exact
    # mirror minimum. Construct angular samples without new QC calculations.
    graph = nx.from_numpy_array(p.bond01)
    graph.remove_edge(*axis)
    moving = list(nx.node_connected_component(graph, axis[1]))
    atoms = Atoms(p.atom, positions=p.geom)
    start = atoms.get_dihedral(*dihedral)
    samples = []
    for index, angle in enumerate([start+30*i for i in range(12)] + [0.]):
        atoms.set_positions(p.geom)
        atoms.set_dihedral(*dihedral, angle, indices=moving)
        samples.append(dict(index=index, status='successful',
            geometry_angstrom=atoms.get_positions(), electronic_energy_hartree=0.,
            angle_offset_degrees=(angle-start) % 360.))
    rotor['points'] = samples
    assert evaluate_optical(p, rotors=[rotor], population='racemic')['remaining_multiplier'] == 1.
    from kinbot.hindered_rotors import HIR
    hir = HIR.__new__(HIR)
    hir.species = p
    p.dihed = [dihedral]
    p.optical_population = 'racemic'
    planar = samples[-1]['geometry_angstrom']
    assert hir._configuration_error(planar, 0) is None
    assert not configured_geometry_allowed(p, planar, 'racemic')
    # The same scan-only exception works for a nonreacting axis on a TS.
    product = copy.copy(p)
    product.bond = p.bond.copy()
    hydrogen = next(i for i, atom in enumerate(p.atom) if atom == 'H')
    carbon = int(np.flatnonzero(p.bond[hydrogen])[0])
    product.bond[carbon, hydrogen] = product.bond[hydrogen, carbon] = 0
    product.bonds = [product.bond.copy()]
    ts = copy.copy(p)
    ts.wellorts = 1
    prepare_ts_context(ts, p, product)
    from kinbot.stereo_identity import rotor_geometry_allowed
    assert rotor_geometry_allowed(ts, planar, axis, 'racemic')
    assert not rotor_geometry_allowed(ts, planar, axis, 'specified')


def test_planar_radical_resonance_drawing_does_not_create_an_axis():
    p = point('C=C[C]=CC')
    p.mult = 2
    assert len(p.bonds) > 1
    # Planar structural control, not a computed radical minimum.
    p.geom[:, 2] = 0.
    identity = canonical_identity(p)
    assert identity['status'] == 'assigned'
    assert not any(':' in graph for graph in identity['canonical_graphs'])


@pytest.mark.parametrize('smiles', [
    'CC=C=Cc1ccccc1-c1ccccc1F', 'C[C@H](F)c1ccccc1-c1ccccc1Cl'])
def test_biaryl_rotation_cannot_supply_another_fixed_centers_mirror(smiles):
    p = point(smiles)
    axis = next((i, j) for i, j in nx.bridges(nx.from_numpy_array(p.bond01))
                if p.cycle[i] and p.cycle[j])
    outer = [next(i for i in np.flatnonzero(p.bond[a]) if i != b) for a, b in (axis, axis[::-1])]
    dihedral = [outer[0], *axis, outer[1]]
    graph = nx.from_numpy_array(p.bond01)
    graph.remove_edge(*axis)
    moving = list(nx.node_connected_component(graph, axis[1]))
    atoms = Atoms(p.atom, positions=p.geom)
    start = atoms.get_dihedral(*dihedral)
    atoms.set_dihedral(*dihedral, 0., indices=moving)
    planar = atoms.get_positions()
    from kinbot.stereo_identity import rotor_geometry_allowed
    assert rotor_geometry_allowed(p, planar, axis, 'racemic')
    assert not rotor_geometry_allowed(p, planar, (0, 1), 'racemic')
    atoms.set_dihedral(*dihedral, -start, indices=moving)
    assert not rotor_geometry_allowed(p, atoms.get_positions(), axis, 'racemic')
    rotor = dict(index=0, axis=axis, dihedral=dihedral, usable=True, sigma_int=1,
        represented_domain_degrees=[0., 360.], points=[dict(index=i, status='successful',
        geometry_angstrom=g, electronic_energy_hartree=0., angle_offset_degrees=angle)
        for i, (g, angle) in enumerate(((p.geom, 0.), (planar, -start % 360.)))])
    assert evaluate_optical(p, rotors=[rotor], population='racemic')['remaining_multiplier'] == 2.
