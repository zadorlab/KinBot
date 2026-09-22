"""MC optical estimates use each retained conformer's own quadratic model."""
import copy
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot import constants
from kinbot.calculation import array_fingerprint
from kinbot.conformer_records import hessian_record, inventory, retain, update_member
from kinbot.conformer_counting import evaluate_members, writer_members, CountingError
from kinbot.optical_harmonic import conformer_midpoint_diagnostic, conformer_pair_midpoint
from kinbot.mess import MESS, finalize_mc_mess
from test_optical_evaluator import ethanol_near_planar, record


def own_hessian(p, r, scale=1.):
    hess = np.eye(3 * p.natom) * scale
    source = f'conf/{p.name}_{r.index:04d}'
    ref = dict(source_job=source, source_row_id=10+r.index, atoms=list(map(str, p.atom)),
               geometry_sha256=array_fingerprint(r.geometry), hessian_sha256=array_fingerprint(hess),
               hessian_unit='hartree / bohr^2', hessian_massweighted=False)
    return replace(r, source_job=source, hessian=tuple(map(tuple, hess)), hessian_reference=ref)


def ambiguous_pair():
    p = ethanol_near_planar(60.)
    other = Atoms(p.atom, positions=p.geom * [-1, 1, 1])
    other.rotate_dihedral(0, 1, 2, 8, 12.5, indices=[2, 8])
    return p, [own_hessian(p, record(p, 0)), own_hessian(p, record(p, 1, other.positions))]


@pytest.mark.parametrize('energy,weight', [(2., 1.), (8., 2.)])
def test_each_unresolved_self_mirror_uses_its_own_hessian_before_pruning(energy, weight):
    p = ethanol_near_planar()
    r = own_hessian(p, record(p, 0))
    base = conformer_midpoint_diagnostic(p, r)['stable_midpoint_energy_kcal_mol']
    r = own_hessian(p, r, energy / base)
    p.hess = np.eye(3 * p.natom) * 1.e9  # unrelated representative must not be used
    records, groups = evaluate_members(p, [r])
    assert groups == [[0]]
    assert records[0].remaining_optical_weight == weight
    evidence = records[0].optical_evidence
    assert evidence['heuristic'] == 'harmonic_midpoint'
    assert evidence['harmonic_midpoint']['stable_midpoint_energy_kcal_mol'] == pytest.approx(energy)
    p.energy, p.zpe, p.reduced_freqs = -100., 0., list(r.frequencies_cm1)
    writer = MESS(dict(multi_conf_tst=1, pes=0, freq_uq_ref=1000., freq_uq_max_exp=1.), p)
    text = writer._member_rrho(p, records[0], 1., 0.)
    assert 'Harmonic midpoint approximation' in text
    assert f'optical factor {weight:g}' in text
    assert next(float(line.split()[1]) for line in text.splitlines()
                if line.strip().startswith('SymmetryFactor')) == records[0].sigma_ext / weight
    assert 'Harmonic midpoint approximation' in finalize_mc_mess(text)


@pytest.mark.parametrize('energy,paired', [(2., True), (8., False)])
def test_ambiguous_explicit_pair_uses_both_hessians(energy, paired):
    p, pair = ambiguous_pair()
    d = conformer_pair_midpoint(p, *pair)
    pair = [own_hessian(p, r, energy/d[k]['stable_midpoint_energy_kcal_mol'])
            for r, k in zip(pair, ('forward', 'reverse'))]
    records, groups = evaluate_members(p, pair)
    assert groups == ([[0, 1]] if paired else [[0], [1]])
    assert [r.remaining_optical_weight for r in records] == ([1., 1.] if paired else [2., 2.])
    assert all(not r.optical_evidence.get('fallback') for r in records)
    if paired:
        comment = MESS._counting_comment(records[0].optical_evidence, records[0].member_id)
        assert 'explicit mirror' in comment and 'cutoff 4' in comment


def test_pair_disagreement_is_reported_and_cannot_double_both_conformers():
    p, pair = ambiguous_pair()
    d = conformer_pair_midpoint(p, *pair)
    pair = [own_hessian(p, r, energy/d[k]['stable_midpoint_energy_kcal_mol'])
            for r, k, energy in zip(pair, ('forward', 'reverse'), (2., 8.))]
    records, groups = evaluate_members(p, pair)
    assert groups == [[0], [1]]
    assert all(r.remaining_optical_weight == 1. for r in records)
    assert 'disagree' in records[0].optical_evidence['explicit_mirror_comparisons'][pair[1].member_id]['reason']


def test_explicit_mirrors_override_a_single_conformer_heuristic_of_two():
    p = ethanol_near_planar()
    first = own_hessian(p, record(p, 0))
    energy = conformer_midpoint_diagnostic(p, first)['stable_midpoint_energy_kcal_mol']
    first = own_hessian(p, first, 8./energy)
    second = own_hessian(p, record(p, 1, p.geom * [-1, 1, 1]), 8./energy)
    records, groups = evaluate_members(p, [first, second])
    assert groups == [[0, 1]]
    assert [r.remaining_optical_weight for r in records] == [1., 1.]
    assert all(r.optical_evidence.get('heuristic') is None for r in records)
    assert all(not r.optical_evidence['harmonic_midpoint']['used_for_optical_counting']
               for r in records)


def test_harmonic_pairing_warns_about_different_energies_without_multiplying_again(caplog):
    p, pair = ambiguous_pair()
    pair = [own_hessian(p, r, .01) for r in pair]
    pair[1] = replace(pair[1], zero_energy_hartree=pair[0].zero_energy_hartree+2/constants.AUtoKCAL)
    records, groups = evaluate_members(p, pair)
    assert groups == [[0, 1]]
    assert all(r.remaining_optical_weight == 1. for r in records)
    assert records[1].zero_energy_hartree == pair[1].zero_energy_hartree
    assert 'disagree' in caplog.text
    assert all(r.optical_evidence['warnings'] for r in records)


@pytest.mark.parametrize('weighted', [False, True])
def test_backend_hessian_units_and_geometry_are_saved_with_the_record(tmp_path, weighted):
    p = ethanol_near_planar()
    db = connect(str(tmp_path/'kinbot.db'))
    row_id = db.write(Atoms(p.atom, positions=p.geom), name='conf/test', data={'status': 'normal'})
    qc = SimpleNamespace(db=db, read_qc_hess=Mock(return_value=np.eye(3*p.natom)),
                         hessian_is_massweighted=lambda: weighted)
    data = hessian_record(qc, 'conf/test', p.geom, p.atom)
    assert data['hessian_reference']['source_row_id'] == row_id
    assert data['hessian_reference']['hessian_massweighted'] is weighted
    assert ('amu' in data['hessian_reference']['hessian_unit']) is weighted
    assert hessian_record(qc, 'conf/test', p.geom + .01, p.atom) == {}
    with patch('kinbot.species_routing.resolve_job', return_value='conf/test'):
        assert hessian_record(qc, 'legacy_alias', p.geom, p.atom)['hessian_reference']['source_job'] == 'legacy_alias'


def test_l1_inventory_l2_replacement_and_writer_preserve_own_hessian():
    p = ethanol_near_planar()
    first = own_hessian(p, record(p, 0))
    source = {'source_job': first.source_job, 'hessian': first.hessian,
              'hessian_reference': first.hessian_reference}
    records = inventory(p, [first.geometry], [-100.], [first.frequencies_cm1], [0], {0: source})
    retain(p, records, [0])
    p.conformer_index = [0]
    p.conformer_geom = [first.geometry]
    p.conformer_energy = [-100.01]
    p.conformer_zeroenergy = [-100.]
    p.conformer_freq = [first.frequencies_cm1]
    rendered = writer_members(p)[0]
    assert rendered.hessian == first.hessian
    assert rendered.hessian_reference == first.hessian_reference
    assert 'hessian' not in rendered.as_dict() and rendered.as_dict()['hessian_available']
    newref = dict(first.hessian_reference, source_job='L2', source_row_id=99)
    update_member(p, 0, source_job='L2', hessian_data=dict(hessian=first.hessian, hessian_reference=newref))
    assert writer_members(p)[0].hessian_reference == newref
    # A missing L2 matrix must never leave the old L1 matrix attached.
    update_member(p, 0, source_job='L2_without_hessian')
    assert writer_members(p)[0].hessian is None


@pytest.mark.parametrize('change', ['geometry', 'job', 'hessian', 'unit'])
def test_stale_conformer_hessian_cannot_be_used(change):
    p = ethanol_near_planar()
    r = own_hessian(p, record(p, 0))
    if change == 'geometry': r = replace(r, geometry=tuple(map(tuple, p.geom + .01)))
    elif change == 'job': r = replace(r, source_job='different')
    elif change == 'hessian': r = replace(r, hessian=tuple(map(tuple, np.eye(3*p.natom)*2)))
    else: r = replace(r, hessian_reference=dict(r.hessian_reference, hessian_unit='unknown'))
    assert conformer_midpoint_diagnostic(p, r)['status'] == 'unavailable'
