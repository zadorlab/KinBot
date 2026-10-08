"""Saved beta-ocimene product regressions; no QC is submitted.

Coordinates are from Judit's UMA-s-1p2p1 run of 2026-10-07.
"""
import copy
import json
from pathlib import Path

import numpy as np

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
