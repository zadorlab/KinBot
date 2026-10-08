"""Saved beta-ocimene product regressions; no QC is submitted.

Coordinates are from Judit's UMA-s-1p2p1 run of 2026-10-07. Synthetic
energies/Hessians below test result association, not thermochemistry.
"""
import copy
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot.calculation import load_calculation_record
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import (relaxed_homolytic_product, routing_name,
                                    configured_result_identity)
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.stereo_routing import guard_well_job, StereoRoutingError, relaxed_radical_configuration
from kinbot.stereochemistry import virtually_labelled
from kinbot.symmetry import calculate_symmetry
from tests.test_pah_stereo_scope import molecule, observation


FIXTURES = json.loads((Path(__file__).parent/'reference/ocimene_products.json').read_text())


def point(name):
    data = FIXTURES[name]
    p = StationaryPoint(name, 0, data['mult'], atom=np.array(data['atoms']),
                        geom=np.array(data['geometry']))
    p.characterize()
    return p


def embedded(smiles, mult=1):
    data = observation(molecule(smiles))
    p = StationaryPoint(smiles, 0, mult, atom=data.atom, geom=data.geom)
    p.characterize()
    return p


def test_terminal_allene_virtual_labels_do_not_reject_the_physical_product():
    p = point('allene_minimum')
    expected = canonical_identity(p)
    assert expected['status'] == 'assigned'
    # Both terminal CH2 hydrogens used to cause an unsupported-axis refusal
    # inside configured_external_labels, despite the ordinary identity passing.
    for h in (10, 11):
        assert canonical_identity(p, tagged_atom=h)['status'] == 'assigned'
        labelled = virtually_labelled(p, {h: 1})
        assert canonical_identity(labelled)['status'] == 'assigned'
        deuterated = copy.copy(p)
        deuterated.isotopes = [0] * p.natom
        deuterated.isotopes[h] = 2
        assert canonical_identity(deuterated)['status'] == 'unsupported'
    calculate_symmetry(p)
    assert p.sigma_ext == 1
    assert canonical_identity(p) == expected


def test_bent_vinyl_resonance_does_not_license_genuine_axes():
    parent = point('parent')
    removed = [3] + [i for i, atom in enumerate(parent.atom)
                     if atom == 'H' and parent.bond01[3, i]]
    keep = [i for i in range(parent.natom) if i not in removed]
    p = StationaryPoint('CH3 loss', 0, 2, atom=parent.atom[keep], geom=parent.geom[keep])
    p.characterize()
    matrices = [b.copy() for b in p.bonds]
    assert p.chemid == 1213614424683883102242
    assert len(p.bonds) == 2
    expected = canonical_identity(p)
    assert expected['status'] == 'assigned'
    assert canonical_identity(p, tagged_atom=0)['status'] == 'assigned'
    moved = copy.copy(p)
    moved.geom = p.geom @ np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]) + [2., 3., 4.]
    moved.bonds = list(reversed(moved.bonds))
    assert canonical_identity(moved) == expected
    for original, after in zip(matrices, p.bonds):
        np.testing.assert_array_equal(original, after)
    calculate_symmetry(p)
    for smiles, mult in [('CC=C=CC', 1), ('[CH2]C(=C=CF)C', 2),
                          ('C1=C=CCCCC1', 1)]:
        axial = embedded(smiles, mult)
        assert canonical_identity(axial)['status'] == 'unsupported', smiles
        axial.geom *= [-1., 1., 1.]
        assert canonical_identity(axial)['status'] == 'unsupported', smiles
    # A bent true axis with non-coplanar terminal groups must not qualify
    # merely because another bond arrangement contains a vinyl radical.
    axial = embedded('[CH2]C(=C=CF)C', 2)
    for twist in (10., 90.):
        atoms = Atoms(axial.atom, positions=axial.geom)
        indices = [3, 4, 8]  # terminal carbon, fluorine and hydrogen
        tail = atoms[indices]
        tail.rotate(twist, v=axial.geom[3]-axial.geom[2], center=axial.geom[2])
        atoms.positions[indices] = tail.positions
        atoms.set_angle(1, 2, 3, 140., indices=indices)
        bent = copy.copy(axial)
        bent.geom = atoms.positions
        assert canonical_identity(bent)['status'] == 'unsupported'
        bent.geom *= [-1., 1., 1.]
        assert canonical_identity(bent)['status'] == 'unsupported'


@pytest.mark.parametrize('target_exists', [False, True, 'reordered'])
def test_initial_homolytic_radical_relaxes_without_weakening_well_checks(tmp_path, monkeypatch, target_exists):
    monkeypatch.chdir(tmp_path)
    initial, final, parent = (point(name) for name in ('h_loss_input', 'h_loss_minimum', 'parent'))
    keep = [i for i in range(parent.natom) if i != 17]  # H18 loss from C6
    parent_bonds = [b[np.ix_(keep, keep)] for b in parent.bonds]
    assert canonical_identity(initial)['id'] != canonical_identity(final)['id']
    assert relaxed_radical_configuration(initial, final, parent_bonds)
    # A pre-existing double bond must still protect the same E/Z change.
    assert not relaxed_radical_configuration(initial, final, initial.bonds)
    import networkx as nx
    graph = nx.from_numpy_array(final.bond01)
    graph.remove_edge(2, 4)  # inherited C3=C5
    rotating = list(nx.node_connected_component(graph, 5))
    opposite = copy.copy(final)
    atoms = Atoms(final.atom, positions=final.geom)
    atoms.set_dihedral(1, 2, 4, 5, atoms.get_dihedral(1, 2, 4, 5) + 180., indices=rotating)
    opposite.geom = atoms.positions
    assert not relaxed_radical_configuration(initial, opposite, parent_bonds)
    chiral = embedded('C[C@H](F)C[CH]C=C', 2)
    mirrored = copy.copy(chiral)
    mirrored.geom = chiral.geom * [-1., 1., 1.]
    assert not relaxed_radical_configuration(chiral, mirrored, [chiral.bond01])
    # A product graph can keep exactly the same adjacency matrix after an
    # atom permutation which exchanges the PARENT single/double bond roles.
    symmetric = embedded('[CH2]/C=C/C=C', 2)
    order = np.array([4, 3, 2, 1, 0, 10, 11, 9, 8, 7, 6, 5])
    reordered = StationaryPoint('reordered', 0, 2,
        atom=np.asarray(symmetric.atom)[order], geom=symmetric.geom[order])
    reordered.characterize()
    np.testing.assert_array_equal(reordered.bond01, symmetric.bond01)
    inherited = next(b for b in symmetric.bonds if b[1, 2] == 2 and b[3, 4] == 2)
    reaction = SimpleNamespace(products=[reordered], product_parent_bonds=[
        (np.asarray(symmetric.atom), symmetric.geom, [inherited])])
    generator = ReactionGenerator(SimpleNamespace(reac_type=['hom_sci']), {}, None, '')
    with patch('kinbot.reaction_generator.relaxed_homolytic_product',
               side_effect=AssertionError('Parent indices changed')):
        assert generator.relaxed_initial_product(reaction, 0, 0) is reordered
    qc = QuantumChemistry.__new__(QuantumChemistry)
    qc.db = connect('kinbot.db')
    qc.par, qc.job_ids = {'optical_population': 'specified'}, {}
    qc.check_qc = Mock(return_value=0)
    qc.read_qc_hess = Mock(side_effect=AssertionError('No additional Hessian read'))
    job, target = routing_name(initial)+'_well', routing_name(final)+'_well'
    guard_well_job(qc, initial, initial.geom, job)
    data = {'status': 'normal', 'charge': 0, 'multiplicity': 2,
            'energy': -100., 'zpe': .01, 'frequencies': [100.]*(3*initial.natom-6),
            'hess': np.eye(3*initial.natom)}
    row = qc.db.write(Atoms(initial.atom, positions=final.geom), name=job, data=data)
    existing = None
    if target_exists:
        saved = copy.copy(final)
        if target_exists == 'reordered':
            order = np.arange(final.natom)
            order[[0, 1]] = order[[1, 0]]
            saved = StationaryPoint('previous product', 0, 2,
                atom=final.atom[order], geom=final.geom[order])
            saved.characterize()
        guard_well_job(qc, saved, saved.geom, target)
        existing = qc.db.write(Atoms(saved.atom, positions=saved.geom), name=target,
                               data=dict(data, energy=-101.))
    # Established wells must refuse this result, before AND after adoption.
    with pytest.raises(StereoRoutingError, match='different stereoisomer'):
        guard_well_job(qc, initial, initial.geom, job)
    # A reused input in a different coordinate/order frame cannot establish
    # which parent bond indices apply, even if its chemical identity matches.
    moved = copy.copy(initial)
    moved.geom = initial.geom + [1., 0., 0.]
    with pytest.raises(StereoRoutingError, match='different stereoisomer'):
        relaxed_homolytic_product(qc, moved, parent_bonds)
    adopted = relaxed_homolytic_product(qc, initial, parent_bonds)
    assert routing_name(adopted) == routing_name(final)
    if target_exists != 'reordered':
        assert adopted.source_job == job and adopted.source_row_id == row
        np.testing.assert_array_equal(adopted.geom, final.geom)
        np.testing.assert_array_equal(adopted.hess, data['hess'])
    target_row = list(qc.db.select(name=target))[-1]
    assert configured_result_identity(qc.db, routing_name(adopted), target_row) == canonical_identity(final)['id']
    guard_well_job(qc, adopted, adopted.geom, target)
    load_calculation_record(adopted, qc, target)
    if target_exists:
        assert qc.db.get(name=target).id == existing
        assert qc.db.get(name=target).data['energy'] == -101.
        np.testing.assert_array_equal(adopted.geom, saved.geom)
    else:
        published = list(qc.db.select(name=target))[-1]
        assert published.data['copied_from_row_id'] == row
        assert published.data['energy'] == data['energy']
        np.testing.assert_array_equal(published.data['hess'], data['hess'])
    with pytest.raises(StereoRoutingError, match='different stereoisomer'):
        guard_well_job(qc, initial, initial.geom, job)
    assert qc.db.get(row).data['energy'] == data['energy']
    qc.read_qc_hess.assert_not_called()
