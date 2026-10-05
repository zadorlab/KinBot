import copy
import numpy as np
from kinbot.reaction_path import prepare_stereopath, prepare_ts_context, path_geometry_allowed
from kinbot.stereo_identity import optical_scope
from kinbot.conformer_records import ConformerRecord
from kinbot.conformer_counting import evaluate_members
from test_reaction_paths import peroxy, transfer


def test_unsplit_ts_rejects_opposite_spectator_configuration_before_mc_counting():
    p = peroxy()
    h = next(i for i in range(p.natom) if p.atom[i] == 'H' and p.bond[0, i])
    ts, q = transfer(p, h, donor=0)
    assert ts.stereopath_id == 'ordinary'
    assert path_geometry_allowed(ts, ts.geom)
    assert not path_geometry_allowed(ts, ts.geom * [-1., 1., 1.])
    assert path_geometry_allowed(ts, ts.geom * [-1., 1., 1.], 'racemic')
    records = [ConformerRecord(str(i), i, str(i), 'valid', geometry=g.tolist(),
                               zero_energy_hartree=-100., frequencies_cm1=tuple(ts.freq))
               for i, g in enumerate((ts.geom, ts.geom * [-1., 1., 1.]))]
    counted, groups = evaluate_members(ts, records)
    assert groups == [[0]]
    assert counted[1].exclusion_reason == 'different stereochemical pathway'


def test_reacting_centre_is_not_frozen_as_a_spectator():
    p = peroxy()
    h = next(i for i in range(p.natom) if p.atom[i] == 'H' and p.bond[2, i])
    q = copy.deepcopy(p)
    q.bond[2,h] = q.bond[h,2] = 0
    q.bonds = [q.bond.copy()]
    ts = copy.copy(p)
    ts.wellorts = 1
    ts.geom = p.geom.copy()
    neighbours = np.flatnonzero(p.bond[2])
    # A planar reacting centre has no defined R/S assignment to freeze.
    ts.geom[[2, *neighbours], 2] = 0.
    prepare_ts_context(ts, p, q)
    assert 2 in ts.configuration_ignored_atoms
    assert path_geometry_allowed(ts, ts.geom * [-1., 1., 1.])


def test_ho2_elimination_keeps_both_path_labels_and_additional_bond_edits():
    p = peroxy()
    ids = []
    for h in [9, 10]:
        q = copy.deepcopy(p)
        q.bond[1,h] = q.bond[h,1] = q.bond[2,4] = q.bond[4,2] = 0
        q.bond[5,h] = q.bond[h,5] = 1
        q.bond[1,2] = q.bond[2,1] = 2
        q.bonds = [q.bond.copy()]
        q.calc_chemid(); q.find_atom_eqv()
        ts = copy.copy(p); ts.wellorts = 1
        prepare_stereopath(ts, p, q)
        ids.append(ts.stereopath_id)
        assert len(ts.stereopath_metadata['broken_bonds']) == 2
        assert len(ts.stereopath_metadata['formed_bonds']) == 2  # includes increasing bond order
    assert len(set(ids)) == 2
