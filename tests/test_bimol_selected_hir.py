"""Bimolecular HIR keeps native and embedded selected Hessians consistent."""
import copy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase.build import molecule

from kinbot import constants, frequencies
from kinbot.bimol import load_selected_species_properties
from kinbot.calculation import geometry_reference, load_calculation_record
from kinbot.hindered_rotors import HIR, recover_hir_model
from kinbot.optimize import Optimize
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.stationary_pt import StationaryPoint


def selected_model(tmp_path, monkeypatch, backend='gauss', high_level=0,
                   use_sella=False):
    monkeypatch.chdir(tmp_path)
    atoms = molecule('CH3OH')
    species = StationaryPoint('methanol', 0, 1,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    species.characterize()
    hessian = np.diag(np.linspace(.01, .02, 3 * species.natom))
    raw, expected = frequencies.get_frequencies(species, hessian, species.geom)
    job = 'methanol_high' if high_level else 'methanol_well'
    data = dict(energy=-100. / constants.EVtoHARTREE, zpe=.05,
                frequencies=raw, status='normal')
    if backend == 'fc' or use_sella:
        data['hess'] = hessian.tolist()
    row = SimpleNamespace(name=job, id=1, symbols=list(species.atom),
                          positions=species.geom.copy(), data=data)
    rows = {job: row}
    qc = QuantumChemistry.__new__(QuantumChemistry)
    qc.qc, qc.use_sella = backend, use_sella
    qc.db = SimpleNamespace(select=lambda name: [rows[name]] if name in rows else [])
    qc.check_qc = Mock(return_value='normal')
    qc.qc_freq = Mock(side_effect=AssertionError('Unexpected frequency job'))
    qc.qc_hir = Mock(side_effect=AssertionError('Unexpected HIR job'))
    qc.publish_result = Mock()
    if backend == 'gauss' and not use_sella:
        lower = hessian[np.tril_indices(len(hessian))]
        Path(job + '.fchk').write_text(
            'Cartesian Force Constants R N= ' + str(len(lower)) + '\n' +
            '\n'.join(' '.join(f'{x:.16e}' for x in lower[i:i + 5])
                      for i in range(0, len(lower), 5)) + '\n')
    elif backend == 'qchem' and not use_sella:
        masses = np.repeat([constants.exact_mass[a] for a in species.atom], 3)
        weighted = hessian / np.sqrt(np.outer(masses, masses))
        Path(job + '_freq.out').write_text(' Mass-Weighted Hessian Matrix\n' +
            '\n'.join(' '.join(f'{x:.16e}' for x in row) for row in weighted) +
            '\n Translations and Rotations\n')
    par = dict(high_level=high_level, conformer_search=0, multi_conf_tst=0,
               rotor_scan=1, imagfreq_threshold=50., rigid_hir=False,
               nrotation=12, rotor_0_test=1, plot_hir_profiles=False, bimol=1,
               rotation_restart=3, L3_calc=0)
    qc.par = par
    opt = Optimize.__new__(Optimize)
    opt.species, opt.par, opt.qc = species, par, qc
    opt.name, opt.selected_job = species.name, job
    opt.shigh = opt.shir = 1
    opt.just_high = opt.defer_hir = False
    opt.wait = 0
    opt.restart = 0
    load_calculation_record(species, qc, job)
    hir = HIR(species, qc, par)
    species.hir = hir
    angles = np.arange(12) * 2 * np.pi / 12
    hir.hir_status = [[0] * 12]
    hir.hir_energies = [(-100. + (1. - np.cos(3 * angles)) /
                        constants.AUtoKCAL).tolist()]
    hir.hir_geoms = [[species.geom.copy() for _ in angles]]
    hir.scan_reference = geometry_reference(species)
    hir.scan_reference['dihedrals'] = copy.deepcopy(species.dihed)
    return opt, rows, np.asarray(raw), np.asarray(expected)


def resume_saved_scan(opt, monkeypatch):
    """Submission finds this fixture's already completed scan calculations."""
    opt.defer_hir = True
    hir = opt.species.hir
    hir.generate_hir_geoms = Mock()
    monkeypatch.setattr('kinbot.optimize.HIR', Mock(return_value=hir))


@pytest.mark.parametrize('backend,use_sella',
                         [('gauss', False), ('qchem', False),
                          ('qchem', True), ('fc', False)])
@pytest.mark.parametrize('high_level', [0, 1])
def test_completed_optimization_retains_projection_after_selected_reload(
        tmp_path, monkeypatch, backend, use_sella, high_level):
    opt, rows, raw, expected = selected_model(
        tmp_path, monkeypatch, backend, high_level, use_sella)
    assert opt._ensure_selected_hessian()
    assert recover_hir_model(opt.species, opt.qc, opt.par)
    reference = copy.deepcopy(opt.species.optical_hessian_reference)
    native_hessian = copy.deepcopy(opt.species.hess)
    projection = copy.deepcopy(opt.species.rotor_projection)
    opt.qc.read_qc_hess = Mock(side_effect=AssertionError('Already loaded'))
    assert load_selected_species_properties(opt, opt.qc) == opt.selected_job
    np.testing.assert_allclose(opt.species.freq, raw)
    np.testing.assert_allclose(opt.species.reduced_freqs, expected, rtol=1.e-8)
    assert len(raw) == 12 and len(expected) == 11
    assert opt.species.hess == native_hessian
    assert opt.species.optical_hessian_reference == reference
    assert opt.species.rotor_projection == projection
    opt.qc.qc_freq.assert_not_called()
    opt.qc.qc_hir.assert_not_called()


@pytest.mark.parametrize('backend', ['gauss', 'qchem'])
def test_changed_record_reads_native_hessian_and_repairs_projection(
        tmp_path, monkeypatch, backend):
    opt, rows, raw, expected = selected_model(tmp_path, monkeypatch, backend)
    assert opt._ensure_selected_hessian()
    assert recover_hir_model(opt.species, opt.qc, opt.par)
    rows[opt.selected_job].id = 2
    rows[opt.selected_job].data['energy'] -= .01
    load_selected_species_properties(opt, opt.qc)
    assert opt.species.source_row_id == 2
    assert opt.species.energy == rows[opt.selected_job].data['energy'] * constants.EVtoHARTREE
    assert opt.species.optical_hessian_reference['source_row_id'] == 2
    assert opt.species.hir.scan_reference['source_row_id'] == 2
    np.testing.assert_allclose(opt.species.reduced_freqs, expected, rtol=1.e-8)
    assert opt._projection_source is None
    opt.qc.qc_freq.assert_not_called()
    opt.qc.qc_hir.assert_not_called()


def test_failed_native_read_does_not_leave_raw_modes_with_usable_rotor(
        tmp_path, monkeypatch, caplog):
    opt, rows, raw, expected = selected_model(tmp_path, monkeypatch)
    rows[opt.selected_job].id = 2
    Path(opt.selected_job + '.fchk').unlink()
    load_selected_species_properties(opt, opt.qc)
    np.testing.assert_allclose(opt.species.reduced_freqs, raw)
    assert not opt.species.hir.is_valid_rotor(0)
    assert 'no full Hessian associated' in caplog.text
    opt.qc.qc_freq.assert_not_called()
    opt.qc.qc_hir.assert_not_called()


@pytest.mark.parametrize('backend,use_sella',
                         [('gauss', False), ('qchem', False),
                          ('qchem', True), ('fc', False)])
@pytest.mark.parametrize('high_level', [0, 1])
def test_final_ts_hir_loads_selected_native_or_embedded_hessian(
        tmp_path, monkeypatch, backend, use_sella, high_level):
    opt, rows, raw, expected = selected_model(
        tmp_path, monkeypatch, backend, high_level, use_sella)
    resume_saved_scan(opt, monkeypatch)
    opt.species.wellorts = 1
    opt.species.hir.check_hir = Mock(return_value=True)
    reaction = SimpleNamespace(ts=opt.species, ts_opt=opt)
    generator = ReactionGenerator(opt.species, opt.par, opt.qc, None)
    generator._record_bimol_status = Mock()
    assert generator._advance_bimol_hir(reaction)
    assert reaction.ts_source_job == opt.selected_job
    assert opt.species.hessian_source_job == opt.selected_job
    np.testing.assert_allclose(opt.species.freq, raw)
    np.testing.assert_allclose(opt.species.reduced_freqs, expected, rtol=1.e-8)
    assert generator._record_bimol_status.call_args.kwargs['final_status'] == 'final_hir_valid'
    opt.qc.qc_freq.assert_not_called()
    opt.qc.qc_hir.assert_not_called()


def test_final_ts_waits_for_selected_frequency_recovery(tmp_path, monkeypatch):
    opt, rows, raw, expected = selected_model(tmp_path, monkeypatch, high_level=1)
    native = opt.qc.read_qc_hess(opt.selected_job, opt.species.natom)
    opt.qc.read_qc_hess = Mock(side_effect=[[], [], native])
    recovery = opt.selected_job + '_recovery'
    rows[recovery] = copy.deepcopy(rows[opt.selected_job])
    rows[recovery].name, rows[recovery].id = recovery, 2
    opt.qc.qc_freq = Mock(return_value=recovery)
    opt.qc.check_qc = Mock(side_effect=['running', 'normal'])
    resume_saved_scan(opt, monkeypatch)
    opt.species.hir.check_hir = Mock(return_value=True)
    reaction = SimpleNamespace(ts=opt.species, ts_opt=opt)
    generator = ReactionGenerator(opt.species, opt.par, opt.qc, None)
    original_geometry = opt.species.geom.copy()
    assert not generator._advance_bimol_hir(reaction)
    opt.species.hir.check_hir.assert_not_called()
    assert generator._advance_bimol_hir(reaction)
    assert reaction.ts_source_job == recovery
    assert opt.species.source_row_id == 2
    assert opt.qc.qc_freq.call_args.kwargs['high_level'] == 1
    np.testing.assert_array_equal(opt.species.geom, original_geometry)
    np.testing.assert_allclose(opt.species.reduced_freqs, expected, rtol=1.e-8)


def test_final_ts_failed_frequency_recovery_uses_harmonic_model(
        tmp_path, monkeypatch, caplog):
    opt, rows, raw, expected = selected_model(tmp_path, monkeypatch)
    opt.qc.read_qc_hess = Mock(return_value=[])
    opt.qc.qc_freq = Mock(return_value=opt.selected_job + '_recovery')
    opt.qc.check_qc = Mock(return_value='error')
    resume_saved_scan(opt, monkeypatch)
    reaction = SimpleNamespace(ts=opt.species, ts_opt=opt)
    generator = ReactionGenerator(opt.species, opt.par, opt.qc, None)
    generator._record_bimol_status = Mock()
    assert generator._advance_bimol_hir(reaction)
    assert 'retaining harmonic frequencies' in caplog.text
    assert not opt.species.hir.is_valid_rotor(0)
    np.testing.assert_allclose(opt.species.reduced_freqs, raw)
    assert generator._record_bimol_status.call_args.kwargs['final_status'] == 'final_hir_harmonic'
    opt.qc.qc_hir.assert_not_called()


def test_final_ts_waits_for_replacement_scan(tmp_path, monkeypatch):
    opt, rows, raw, expected = selected_model(tmp_path, monkeypatch)
    hir = opt.species.hir
    hir.scan_input_errors[0] = 'changed rotor definition'
    hir.check_hir = Mock(return_value=False)
    def replacement(*args, **kwargs):
        hir.hir_status = [[-1] * 12]
    hir.generate_hir_geoms = Mock(side_effect=replacement)
    reaction = SimpleNamespace(ts=opt.species, ts_opt=opt)
    generator = ReactionGenerator(opt.species, opt.par, opt.qc, None)
    generator._record_bimol_status = Mock()
    assert not generator._advance_bimol_hir(reaction)
    assert hir.generate_hir_geoms.call_args.kwargs['rotors'] == [0]
    assert generator._record_bimol_status.call_args.kwargs['final_status'] == 'final_hir_running'


@pytest.mark.parametrize('high_level', [0, 1])
@pytest.mark.parametrize('outcome', ['accepted', 'optimizer_failed', 'wrong_product'])
def test_lower_ts_scan_point_uses_shared_refinement_and_checks_new_irc(
        tmp_path, monkeypatch, high_level, outcome):
    """Run the actual Optimize loop; only QC results/submissions are simulated."""
    from kinbot import bimol, reaction_generator, stereochemistry
    opt, rows, raw, expected = selected_model(
        tmp_path, monkeypatch, backend='fc', high_level=high_level)
    species, qc = opt.species, opt.qc
    species.wellorts = 1
    opt.par.update(conformer_search=1, l0_conformer_search=0,
                   calc_aie=0, print_conf=0)
    # The original conformer search is complete and must not run again.
    opt.scycconf = opt.ssemi_empconf = opt.sconf = 1
    opt.skip_conf_check = 0
    species.confs = SimpleNamespace(check_conformers=Mock(
        side_effect=AssertionError('Conformer search repeated')))
    opt.l1_reference_job = opt.selected_job
    row = rows[opt.selected_job]
    row.data['frequencies'][0] = -800.
    load_calculation_record(species, qc, opt.selected_job)
    original = species.geom.copy()
    lower_point = original.copy()
    lower_point[2, 0] += .03
    refined = lower_point.copy()
    refined[2, 0] += .01
    hir = species.hir
    hir.hir_energies[0][6] = -100. - 1.4 / constants.AUtoKCAL
    hir.hir_geoms[0][6] = lower_point.copy()
    hir.check_hir = Mock(return_value=True)
    hir.generate_hir_geoms = Mock()

    new_hir = HIR(species, qc, opt.par)
    def completed_rescan(cart, rigid):
        np.testing.assert_array_equal(cart, refined)
        angles = np.arange(12) * np.pi / 6
        new_hir.hir_status = [[0] * 12]
        new_hir.hir_energies = [(species.energy +
            (1. - np.cos(3 * angles)) / constants.AUtoKCAL).tolist()]
        new_hir.hir_geoms = [[cart.copy() for _ in angles]]
        new_hir.scan_reference = geometry_reference(species)
        new_hir.scan_reference['dihedrals'] = copy.deepcopy(species.dihed)
    new_hir.generate_hir_geoms = Mock(side_effect=completed_rescan)
    new_hir.check_hir = Mock(return_value=True)
    constructor = Mock(side_effect=[hir, new_hir])
    monkeypatch.setattr('kinbot.optimize.HIR', constructor)
    restart_job = opt.log_name(1) if high_level else opt.log_name(0) + '_hir_restart_1'
    statuses = {opt.selected_job: 'normal'}
    qc.invalidate_qc = Mock(side_effect=lambda job: statuses.update({job: 0}))
    qc.qc_opt_ts = Mock(side_effect=lambda *args, **kw:
                        statuses.update({restart_job: 'running'}))
    qc.check_qc = Mock(side_effect=lambda job: statuses.get(job, 0))
    qc.get_qc_geom = Mock(side_effect=lambda job, *args, **kw:
        (0, lower_point.copy() if job == hir.point_job(0, 6) else refined.copy()))
    qc.get_qc_freq = Mock(return_value=(0, row.data['frequencies']))
    monkeypatch.setattr('kinbot.optimize.reader_sella.read_imag_mode',
                        lambda *args: np.arange(3 * species.natom).reshape(-1, 3))
    opt.defer_hir = opt.defer_cleanup = True
    reaction = SimpleNamespace(instance_name='ts_abstraction', ts=species, par=opt.par,
        ts_opt=opt, final_irc_processed=True, final_irc_geometry=original.copy(),
        sequential_pose_search=True, prepare_hir_rotors=Mock())
    reaction.final_irc_name = bimol.final_irc_job_name(reaction)
    old_irc_name = reaction.final_irc_name
    assert f'_final_l{2 if high_level else 1}_' in old_irc_name
    parent = SimpleNamespace(reac_ts_done=[7])
    generator = ReactionGenerator(parent, opt.par, qc, None)
    generator._record_bimol_status = Mock()
    update = Mock()
    monkeypatch.setattr(bimol, 'update_bimol_status', update)
    monkeypatch.setattr(bimol, 'update_ts_topology_from_irc', Mock())
    monkeypatch.setattr(stereochemistry, 'assign_optical_metadata', lambda ts: {})
    monkeypatch.setattr(stereochemistry, 'conformer_symmetry', lambda ts: {})
    irc = Mock()
    irc.check_irc.side_effect = [[0, 0], ['running', 'running'], ['normal', 'normal']]
    irc.irc2stationary_pt.return_value = SimpleNamespace(
        reactant=SimpleNamespace(), product=SimpleNamespace(geom=refined),
        reactant_direction='R', product_direction='F',
        product_matches_expected=outcome != 'wrong_product')
    reaction.actual_bond_changes = {'formed': [[0, 1]], 'broken': [[1, 2]]}
    irc_constructor = Mock(return_value=irc)
    monkeypatch.setattr(reaction_generator, 'IRC', irc_constructor)

    # A measured low point starts a restart; it is not yet an accepted TS.
    generator._advance_irc(0, reaction, 10.)
    assert opt.restart == 1 and opt.shigh == -1
    assert parent.reac_ts_done == [7]
    assert qc.invalidate_qc.call_count == 12 + high_level
    assert reaction.ts.hir is hir
    generator._advance_irc(0, reaction, 10.)
    assert opt.shigh == 0
    assert qc.qc_opt_ts.call_args.kwargs['high_level'] == high_level
    assert qc.check_qc.call_args.args[0] == restart_job
    np.testing.assert_array_equal(qc.qc_opt_ts.call_args.args[1], lower_point)
    irc_constructor.assert_not_called()
    assert constructor.call_count == 1

    if outcome == 'optimizer_failed':
        statuses[restart_job] = 'error'
        generator._advance_irc(0, reaction, 10.)
        assert parent.reac_ts_done == [-999]
        assert not reaction.final_irc_processed
        assert generator._record_bimol_status.call_args.kwargs['final_status'] == 'ts_hir_refinement_failed'
        irc_constructor.assert_not_called()
        return

    # Accepted refinement and its fresh scan share the new geometry/Hessian.
    rows[restart_job] = copy.deepcopy(row)
    rows[restart_job].name, rows[restart_job].id = restart_job, 2
    rows[restart_job].positions = refined.copy()
    rows[restart_job].data['energy'] = -100.003 / constants.EVtoHARTREE
    rows[restart_job].data['hess'] = (np.eye(3 * species.natom) * .02).tolist()
    statuses[restart_job] = 'normal'
    generator._advance_irc(0, reaction, 10.)
    assert opt.shir == opt.shigh == 1
    assert species.hir is new_hir
    assert species.source_job == reaction.ts_source_job == restart_job
    assert species.source_row_id == 2
    assert species.optical_hessian_reference['source_row_id'] == 2
    assert new_hir.scan_reference['source_row_id'] == 2
    assert len(species.reduced_freqs) == len(species.freq) - 1
    new_hir.generate_hir_geoms.assert_called_once()
    assert reaction.final_irc_name != old_irc_name
    assert not reaction.final_irc_processed
    assert parent.reac_ts_done == [7]
    assert generator._record_bimol_status.call_args.kwargs['final_irc_valid'] is False

    # Submit and poll only the new IRC. Do not release product work early.
    generator._advance_irc(0, reaction, 10.)
    generator._advance_irc(0, reaction, 10.)
    assert parent.reac_ts_done == [7]
    irc.do_irc_calculations.assert_called_once()
    generator._advance_irc(0, reaction, 10.)
    assert parent.reac_ts_done == [2]
    np.testing.assert_array_equal(reaction.final_irc_geometry, refined)
    assert reaction.final_irc_processed
    assert reaction.product_matches_expected == (outcome != 'wrong_product')
    reaction.prepare_hir_rotors.assert_called_once()
    assert constructor.call_count == 2
    species.confs.check_conformers.assert_not_called()


@pytest.mark.parametrize('rotor_scan,multi_conf_tst', [(0, 0), (0, 1), (1, 0)])
def test_unchanged_selected_ts_does_not_repeat_irc(
        tmp_path, monkeypatch, rotor_scan, multi_conf_tst):
    opt, rows, raw, expected = selected_model(tmp_path, monkeypatch)
    opt.par.update(rotor_scan=rotor_scan, multi_conf_tst=multi_conf_tst)
    reaction = SimpleNamespace(ts=opt.species, ts_opt=opt,
        final_irc_processed=True, final_irc_geometry=opt.species.geom.copy())
    parent = SimpleNamespace(reac_ts_done=[7])
    generator = ReactionGenerator(parent, opt.par, opt.qc, None)
    irc = Mock(side_effect=AssertionError('Unexpected repeat IRC'))
    monkeypatch.setattr('kinbot.reaction_generator.IRC', irc)
    generator._advance_irc(0, reaction, 10.)
    assert parent.reac_ts_done == [2]
    irc.assert_not_called()


def test_bimol_keeps_calculation_files_after_hir_until_export(tmp_path, monkeypatch):
    opt, rows, raw, expected = selected_model(tmp_path, monkeypatch)
    opt.par['delete_intermediate_files'] = 1
    opt.defer_cleanup = True
    opt.defer_hir = False
    checkpoint = Path(opt.name + '.chk')
    checkpoint.write_text('selected TS needed for final validation')
    opt.delete_files()
    assert checkpoint.read_text() == 'selected TS needed for final validation'
