"""Saved beta-ocimene product regressions; no QC is submitted.

Coordinates are from Judit's UMA-s-1p2p1 run of 2026-10-07.
"""
import copy
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from tests.test_pah_stereo_scope import molecule, observation

from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.stereochemistry import virtually_labelled
from kinbot.symmetry import calculate_symmetry


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
