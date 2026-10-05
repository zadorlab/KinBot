"""Changed HIR inputs recover without relabelling old scans or losing modes."""
import copy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from kinbot import frequencies
from kinbot.calculation import geometry_reference, array_fingerprint
from kinbot.hindered_rotors import recover_hir_model
from kinbot.qc import QuantumChemistry
from kinbot.optimize import Optimize
from kinbot.stereo_routing import StereoRoutingError
from test_thermochemistry_evidence import scanned_point


def model(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('hir').mkdir()
    species, _ = scanned_point()
    species.source_row_id = 1
    species.hess = np.eye(3 * species.natom).tolist()
    species.hessian_source_job = species.source_job
    hir = species.hir
    hir.scan_reference = geometry_reference(species)
    hir.scan_reference['dihedrals'] = copy.deepcopy(species.dihed)
    hir.scan_jobs = [[f'hir/old_{i}' for i in range(12)]]
    qc = SimpleNamespace(qc='fc', hessian_is_massweighted=lambda: False,
        _check_hir_definition=QuantumChemistry._check_hir_definition,
        read_qc_hess=Mock(return_value=species.hess),
        qc_hir=Mock(side_effect=lambda point, geom, rotor, angle, fix, rigid:
                    f'hir/replacement_{rotor}_{angle}'),
        get_qc_geom=Mock(side_effect=lambda *args: (0, species.geom.copy())),
        get_qc_energy=Mock(return_value=(0, species.energy)))
    hir.qc = qc
    species.kinbot_freqs, species.reduced_freqs = frequencies.get_frequencies(
        species, np.asarray(species.hess), species.geom)
    species.optical_hessian_reference = copy.deepcopy(species.rotor_projection['reference'])
    return species, qc, dict(rotor_scan=1, multi_conf_tst=0, rigid_hir=False)


def test_matching_geometry_repairs_source_metadata_without_qc(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.hir.scan_reference['source_job'] = 'old_alias'
    species.hir.scan_reference['source_row_id'] = 99
    assert recover_hir_model(species, qc, par)
    qc.qc_hir.assert_not_called()
    assert species.hir.scan_reference['source_job'] == species.source_job
    assert species.hir.scan_reference['source_row_id'] == 1
    assert len(species.reduced_freqs) == 11
    assert not recover_hir_model(species, qc, par)


def test_healthy_model_is_unchanged(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    projection = copy.deepcopy(species.rotor_projection)
    assert not recover_hir_model(species, qc, par)
    assert species.rotor_projection == projection
    qc.qc_hir.assert_not_called()
    qc.read_qc_hess.assert_not_called()


@pytest.mark.parametrize('failed', [False, True])
def test_changed_input_recomputes_then_preserves_modes_and_stops_retrying(
        tmp_path, monkeypatch, failed):
    species, qc, par = model(tmp_path, monkeypatch)
    old = copy.deepcopy(species.hir.scan_reference)
    Path('hir/old_0.out').write_text('original observation')
    species.geom[2, 0] += .01
    species.optical_hessian_reference.update(geometry_reference(species))
    if failed:
        qc.get_qc_geom.side_effect = lambda *args: (-1, np.zeros_like(species.geom))
    assert recover_hir_model(species, qc, par, allow_qc=True)
    assert qc.qc_hir.call_count == 12
    assert species.hir.recovery_observations[0]['scan_reference'] == old
    assert Path('hir/old_0.out').read_text() == 'original observation'
    assert all(call.args[0].startswith('hir/replacement_')
               for call in qc.get_qc_geom.call_args_list)
    assert len(species.reduced_freqs) == (12 if failed else 11)
    assert species.rotor_projection['internal_rank'] == (0 if failed else 1)
    if failed:
        assert species.reduced_freqs == species.freq
    assert recover_hir_model(species, qc, par) is failed
    assert qc.qc_hir.call_count == 12


def test_missing_hessian_keeps_full_harmonic_spectrum(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.hess = []
    species.rotor_projection['internal_rank'] = 99
    qc.read_qc_hess.return_value = []
    assert recover_hir_model(species, qc, par)
    assert species.reduced_freqs == species.freq
    assert species.rotor_projection['internal_rank'] == 0
    assert not species.hir.is_valid_rotor(0)
    qc.qc_hir.assert_not_called()


def test_projection_recovery_uses_changed_raw_hessian_without_new_scan(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.hess = (2. * np.asarray(species.hess)).tolist()
    species.optical_hessian_reference['hessian_sha256'] = array_fingerprint(species.hess)
    old = list(species.reduced_freqs)
    assert recover_hir_model(species, qc, par)
    np.testing.assert_allclose(species.reduced_freqs, np.sqrt(2.) * np.asarray(old))
    qc.qc_hir.assert_not_called()


def test_invalid_projection_indices_are_rebuilt(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.rotor_projection['rotors'][0]['rotor_index'] = 44
    assert recover_hir_model(species, qc, par)
    assert species.rotor_projection['rotors'][0]['rotor_index'] == 0
    qc.qc_hir.assert_not_called()


def test_unknown_selected_source_falls_back_without_repeated_jobs(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.source_row_id = None
    assert recover_hir_model(species, qc, par)
    assert species.reduced_freqs == species.freq
    assert not species.hir.is_valid_rotor(0)
    recover_hir_model(species, qc, par)
    qc.qc_hir.assert_not_called()
    qc.read_qc_hess.assert_not_called()


def test_configured_cache_rejection_marks_only_this_optimization_failed():
    opt = Optimize.__new__(Optimize)
    opt.name = 'incompatible_product'
    opt._do_optimization = Mock(side_effect=StereoRoutingError('different specified stereoisomer'))
    assert opt.do_optimization() == 0
    assert opt.shigh == -999
    opt._do_optimization.side_effect = ValueError('unrelated implementation error')
    with pytest.raises(ValueError, match='unrelated'):
        opt.do_optimization()


def test_recovery_exception_restores_all_modes_and_preserves_observations(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    old_geometries = copy.deepcopy(species.hir.hir_geoms)
    species.geom[2, 0] += .01
    qc.qc_hir.side_effect = RuntimeError('backend cannot run this replacement')
    assert recover_hir_model(species, qc, par, allow_qc=True)
    assert len(species.reduced_freqs) == len(species.freq)
    assert species.rotor_projection['internal_rank'] == 0
    assert not species.hir.is_valid_rotor(0)
    np.testing.assert_equal(species.hir.recovery_observations[0]['hir_geoms'], old_geometries)
    attempts = qc.qc_hir.call_count
    recover_hir_model(species, qc, par)
    assert qc.qc_hir.call_count == attempts


def test_missing_scan_geometry_keeps_hir_without_optical_qc(tmp_path, monkeypatch, caplog):
    from kinbot.counting_contract import optical_counting
    from kinbot.thermochemistry import hir_evidence
    species, qc, par = model(tmp_path, monkeypatch)
    species.hir.hir_geoms[0][4] = None
    for allow_qc in (False, True):
        assert not recover_hir_model(species, qc, par, allow_qc=allow_qc)
        assert species.hir.is_valid_rotor(0)
        assert len(species.reduced_freqs) == len(species.freq) - 1
        qc.qc_hir.assert_not_called()
    result = optical_counting(species, hir_evidence(species))
    assert result['remaining_multiplier'] == 1.
    assert result['fallback'] == 'unresolved_symmetry'
    assert 'omitting hindered rotors' not in caplog.text


def test_usable_partial_scan_remains_a_rotor(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.hir.hir_status[0][4] = 1
    angles = np.arange(12) * 2. * np.pi / 12
    species.hir.fourier_fit('partial', angles, 0)
    assert not recover_hir_model(species, qc, par)
    assert species.hir.is_valid_rotor(0)
    assert len(species.reduced_freqs) == len(species.freq) - 1
    qc.qc_hir.assert_not_called()


@pytest.mark.parametrize('failure', ['input', 'energy', 'raw_energy', 'replacement'])
def test_one_invalid_rotor_retains_other_partial_scan(tmp_path, monkeypatch, failure):
    from ase.build import molecule
    from kinbot import symmetry
    from kinbot.hindered_rotors import HIR
    from kinbot.stationary_pt import StationaryPoint

    _, qc, par = model(tmp_path, monkeypatch)
    atoms = molecule('CH3CH2OH')
    species = StationaryPoint('ethanol', 0, 1,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    species.characterize()
    symmetry.calculate_symmetry(species)
    assert len(species.dihed) == 2
    species.source_job, species.source_row_id = 'ethanol_well', 1
    species.energy, species.zpe = -100., .05
    qc.get_qc_geom.side_effect = lambda *args: (0, species.geom.copy())
    species.freq = [500.] * 21
    species.hess = np.eye(3 * species.natom).tolist()
    species.hessian_source_job = species.source_job
    hir = species.hir = HIR(species, qc,
        dict(nrotation=12, plot_hir_profiles=False, rotor_0_test=1))
    hir.hir_status = [[0] * 12 for _ in species.dihed]
    hir.hir_energies = [[species.energy] * 12 for _ in species.dihed]
    hir.hir_geoms = [[species.geom.copy() for _ in range(12)] for _ in species.dihed]
    hir.scan_reference = dict(geometry_reference(species), dihedrals=copy.deepcopy(species.dihed))
    species.kinbot_freqs, species.reduced_freqs = frequencies.get_frequencies(
        species, np.asarray(species.hess), species.geom)
    assert len(species.reduced_freqs) == 19
    hir.hir_status[0][4] = 1
    hir.fourier_fit('partial', np.arange(12) * np.pi / 6, 0)
    saved = copy.deepcopy((hir.hir_status[0], hir.hir_energies[0], hir.hir_geoms[0]))
    if failure in ('input', 'replacement'):
        hir.scan_reference['dihedrals'][1][0] = -1
    elif failure == 'raw_energy':
        hir.fourier_fit('complete', np.arange(12) * np.pi / 6, 1)
        hir.hir_raw_energies[1][3] = np.nan
    else:
        hir.hir_energies[1][3] = np.nan
    replacement = failure == 'replacement'
    assert recover_hir_model(species, qc, par, allow_qc=replacement)
    assert hir.is_valid_rotor(0)
    assert hir.is_valid_rotor(1) == replacement
    assert len(species.reduced_freqs) == (19 if replacement else 20)
    assert species.rotor_projection['internal_rank'] == (2 if replacement else 1)
    assert hir.hir_status[0] == saved[0]
    assert hir.hir_energies[0] == saved[1]
    np.testing.assert_array_equal(hir.hir_geoms[0], saved[2])
    assert qc.qc_hir.call_count == (12 if replacement else 0)
    assert all(call.args[2] == 1 for call in qc.qc_hir.call_args_list)
    qc.read_qc_hess.assert_not_called()


def test_optimization_polls_replacement_without_blocking_or_publishing(tmp_path, monkeypatch):
    from kinbot.parameters import Parameters

    species, qc, par = model(tmp_path, monkeypatch)
    Path('input.json').write_text('{"barrier_threshold": 100}')
    par = dict(Parameters('input.json').par, **par)
    par.update(conformer_search=0, high_level=0, L3_calc=0)
    opt = Optimize(species, par, qc, wait=0)
    opt.shigh, opt.shir = 1, 1
    monkeypatch.setattr(opt, '_ensure_selected_hessian', lambda: True)
    species.geom[2, 0] += .01
    species.optical_hessian_reference.update(geometry_reference(species))
    qc.get_qc_geom.side_effect = lambda *args: (1, None)
    publish = Mock()
    monkeypatch.setattr('kinbot.optimize.publish_optimization_result', publish)
    opt.do_optimization()
    assert opt.shir == 0
    assert qc.qc_hir.call_count == 12
    publish.assert_not_called()
    qc.get_qc_geom.side_effect = lambda *args: (0, species.geom.copy())
    opt.do_optimization()
    assert opt.shir == 1
    assert qc.qc_hir.call_count == 12
    publish.assert_called_once()


def test_recovery_does_not_hide_unrelated_programming_errors(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    for error in (AttributeError('unexpected attribute'), IndexError('unexpected index')):
        monkeypatch.setattr(species.hir, 'is_valid_rotor', Mock(side_effect=error))
        with pytest.raises(type(error), match='unexpected'):
            recover_hir_model(species, qc, par)
    assert len(species.reduced_freqs) == 11


@pytest.mark.parametrize('failure', ['input', 'missing_geometry'])
def test_mess_export_does_not_submit_or_wait_for_qc(tmp_path, monkeypatch, failure):
    import json
    from kinbot.mess import MESS
    from kinbot.parameters import Parameters

    species, qc, par = model(tmp_path, monkeypatch)
    par.update(pes=0, high_level=0, me=0, uq=0, epsilon=100., sigma=3., barrier_threshold=100.,
               high_level_method='b3lyp', high_level_basis='6-31G')
    Path('input.json').write_text(json.dumps(par))
    par = Parameters('input.json').par
    species.reac_obj, species.reac_ts_done = [], []
    Path('me').mkdir()
    if failure == 'input':
        species.hir.scan_reference['dihedrals'][0][0] = -1
    else:
        species.hir.hir_geoms[0][4] = None
    for method in ('qc_hir', 'read_qc_hess', 'get_qc_geom', 'get_qc_energy'):
        getattr(qc, method).side_effect = AssertionError('MESS must not access QC')
    monkeypatch.setattr(species.hir, 'check_hir',
                        Mock(side_effect=AssertionError('MESS must not wait for QC')))
    MESS(par, species).write_input(qc)
    output = Path('me/mess_0000.inp').read_text()
    if failure == 'input':
        assert len(species.reduced_freqs) == 12
        assert 'Hindered' not in output
    else:
        assert len(species.reduced_freqs) == 11
        assert 'Hindered' in output
        assert species.mess_optical_counting['remaining_multiplier'] == 1.
        assert species.mess_optical_counting['fallback'] == 'unresolved_symmetry'
