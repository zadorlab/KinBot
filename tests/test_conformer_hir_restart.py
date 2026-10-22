"""A restarted conformer search must retain a completed lower HIR refinement."""
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from kinbot import constants
from kinbot.calculation import load_calculation_record
from kinbot.conformers import Conformers
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity


@pytest.fixture(params=[0, 1], ids=['well', 'TS'])
def completed(tmp_path, monkeypatch, request):
    monkeypatch.chdir(tmp_path)
    Path('conf').mkdir()
    Path('input.json').write_text('{"barrier_threshold":50, "conformer_search":1, "rotor_scan":0}')
    par = Parameters('input.json', show_warnings=False).par
    par['qc'] = 'fc'
    atoms = molecule('CH3OH')
    sp = StationaryPoint('methanol', 0, 1,
                         atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    sp.characterize()
    sp.optical_reference = canonical_identity(sp)
    sp.wellorts = request.param
    sp.name = 'methanol_ts' if sp.wellorts else str(sp.chemid)
    sp.freq = ([-1000.] if sp.wellorts else [100.]) + [200.] * 11
    qc = QuantumChemistry(par)
    qc.check_qc = Mock(return_value='normal')
    par['bimol'] = 1
    search = Conformers(sp, par, qc)
    sp.confs = search
    search.conf, search.conf_status = 1, [0]
    parent = sp.name if sp.wellorts else sp.name + '_well'
    data = dict(energy=-100. / constants.EVtoHARTREE, zpe=.01,
                frequencies=sp.freq, status='normal')
    for name in (parent, search.get_job_name(0)):
        qc.db.write(atoms, name=name, data=data)
    refined = atoms.copy()
    refined.positions[0, 0] += .002
    data = dict(data, energy=-100.002 / constants.EVtoHARTREE,
                hess=np.eye(3 * sp.natom).tolist())
    source = parent + '_hir_restart_1'
    qc.db.write(refined, name=source, data=data)
    row = list(qc.db.select(name=source))[-1]
    qc.publish_result(row, 'conf/' + sp.name + '_low')
    return sp, par, qc, search, refined


def test_repeated_selection_preserves_complete_refinement(completed):
    sp, par, qc, search, refined = completed
    job = 'conf/' + sp.name + '_low'
    accepted = list(qc.db.select(name=job))[-1]
    for _ in range(2):
        result = search.check_conformers()
        assert result[1] == 'low' and search.selected_job == job
        assert result[3] == pytest.approx(-100.002)
        np.testing.assert_array_equal(result[2], refined.positions)
        assert list(qc.db.select(name=job))[-1].id == accepted.id
        load_calculation_record(sp, qc, search.selected_job)
        assert sp.hess == accepted.data['hess']
    assert search.hir_restart_count == 1


def test_all_failed_samples_do_not_replace_refinement_with_parent(completed):
    sp, par, qc, search, refined = completed
    search.conf_status = [1]
    result = search.check_conformers()
    assert result[3] == pytest.approx(-100.002)
    np.testing.assert_array_equal(result[2], refined.positions)


def test_refinement_with_recovered_native_frequencies_is_retained(completed):
    sp, par, qc, search, refined = completed
    job = 'conf/' + sp.name + '_low'
    row = list(qc.db.select(name=job))[-1]
    data = dict(row.data)
    data['copied_from_job'] += '_freq_recovery_31'
    qc.db.write(refined, name=job, data=data)
    assert search.check_conformers()[3] == pytest.approx(-100.002)
    assert search.hir_restart_count == 1


def test_new_lower_conformer_still_wins(completed):
    sp, par, qc, search, refined = completed
    search.conf, search.conf_status = 2, [0, 0]
    qc.db.write(Atoms(sp.atom, positions=sp.geom), name=search.get_job_name(1),
                data=dict(energy=-100.004 / constants.EVtoHARTREE, zpe=.01,
                          frequencies=sp.freq, status='normal'))
    assert search.check_conformers()[3] == pytest.approx(-100.004)
    assert search.selected_job == search.get_job_name(1)
    # Do not reuse the old refinement job if another HIR restart is necessary.
    assert search.hir_restart_count == 1


@pytest.mark.parametrize('source', ['other_hir_restart_1', 'methanol_high', None])
def test_unrelated_or_ordinary_low_record_is_not_a_refinement(completed, source):
    sp, par, qc, search, refined = completed
    job = 'conf/' + sp.name + '_low'
    row = list(qc.db.select(name=job))[-1]
    data = dict(row.data)
    data.pop('copied_from_job', None)
    if source:
        data['copied_from_job'] = source
    qc.db.write(refined, name=job, data=data)
    assert search.check_conformers()[3] == pytest.approx(-100.)
    assert search.hir_restart_count == 0


def test_resumed_optimization_uses_refinement_without_resubmitting(completed):
    sp, par, qc, search, refined = completed
    opt = Optimize(sp, par, qc, defer_hir=True)
    opt.sconf = opt.scycconf = opt.ssemi_empconf = 1
    qc.qc_opt = Mock(side_effect=AssertionError('Unexpected optimization'))
    qc.qc_opt_ts = Mock(side_effect=AssertionError('Unexpected TS optimization'))
    opt.do_optimization()
    assert opt.restart == 1 and opt.shigh == 1
    assert sp.energy == pytest.approx(-100.002)
    np.testing.assert_array_equal(sp.geom, refined.positions)
    qc.qc_opt.assert_not_called()
    qc.qc_opt_ts.assert_not_called()


def test_l2_stage_is_not_skipped(completed):
    sp, par, qc, search, refined = completed
    par['high_level'] = 1
    opt = Optimize(sp, par, qc, defer_hir=True)
    opt.sconf = opt.scycconf = opt.ssemi_empconf = 1
    qc.qc_opt = Mock()
    qc.qc_opt_ts = Mock()
    qc.check_qc = Mock(side_effect=lambda job, **kwargs: 'running' if job.endswith('_high') else 'normal')
    opt.do_optimization()
    assert opt.restart == 0 and opt.shigh == 0
    (qc.qc_opt_ts if sp.wellorts else qc.qc_opt).assert_called_once()
    assert sp.energy == pytest.approx(-100.002)


def test_ordinary_conformer_selection_does_not_restore_bimolecular_hir(completed):
    sp, par, qc, search, refined = completed
    par['bimol'] = 0
    ordinary = Conformers(sp, par, qc)
    ordinary.conf, ordinary.conf_status = 1, [0]
    ordinary._keep_hir_result = Mock(side_effect=AssertionError(
        'Ordinary conformer selection must not compare a bimolecular refinement'))
    assert ordinary.check_conformers()[3] == pytest.approx(-100.)
    assert ordinary.selected_job == ordinary.get_job_name(0)
    assert ordinary.hir_restart_count == 0


@pytest.mark.parametrize('parent_zpe,sample_zpe,keep_refinement', [
    (.01, .03, False), (.03, .01, True)])
def test_parent_and_refinement_use_their_own_zpe(
        completed, parent_zpe, sample_zpe, keep_refinement):
    sp, par, qc, search, refined = completed
    parent = sp.name if sp.wellorts else sp.name + '_well'
    for job, energy, zpe in ((parent, -100., parent_zpe),
                            (search.get_job_name(0), -99.98, sample_zpe)):
        row = list(qc.db.select(name=job))[-1]
        qc.db.write(row.toatoms(), name=job, data=dict(
            row.data, energy=energy / constants.EVtoHARTREE, zpe=zpe))
    published = 'conf/' + sp.name + '_low'
    row = list(qc.db.select(name=published))[-1]
    qc.db.write(row.toatoms(), name=published, data=dict(row.data, zpe=.02))
    search.check_conformers()
    assert search.selected_job == (published if keep_refinement else parent)


@pytest.mark.parametrize('bimol', [0, 1])
def test_first_rejected_ts_conformer_keeps_selection_change_bimolecular(completed, bimol):
    sp, par, qc, search, refined = completed
    if not sp.wellorts:
        pytest.skip('This comparison uses transition-state frequencies.')
    par['bimol'] = bimol
    search = Conformers(sp, par, qc)
    search.conf, search.conf_status = 3, [0, 0, 0]
    for index, energy, imaginary in ((0, -100.006, -1400.),
                                     (1, -100.004, -1000.),
                                     (2, -100., -1000.)):
        qc.db.write(Atoms(sp.atom, positions=sp.geom), name=search.get_job_name(index),
                    data=dict(energy=energy / constants.EVtoHARTREE, zpe=.01,
                              frequencies=[imaginary] + [200.] * 11,
                              status='normal'))
    search.check_conformers()
    # Bimolecular selection can use the second sample after the first fails
    # its imaginary-mode check. Ordinary selection keeps the base behavior.
    assert search.selected_job == (search.get_job_name(1) if bimol else sp.name)
