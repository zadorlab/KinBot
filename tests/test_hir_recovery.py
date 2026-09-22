"""Changed HIR inputs recover without relabelling old scans or losing modes."""
import copy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from kinbot import frequencies
from kinbot.calculation import geometry_reference
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
    if failed:
        qc.get_qc_geom.side_effect = lambda *args: (-1, np.zeros_like(species.geom))
    assert recover_hir_model(species, qc, par)
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
    assert recover_hir_model(species, qc, par)
    assert len(species.reduced_freqs) == len(species.freq)
    assert species.rotor_projection['internal_rank'] == 0
    assert not species.hir.is_valid_rotor(0)
    np.testing.assert_equal(species.hir.recovery_observations[0]['hir_geoms'], old_geometries)
    attempts = qc.qc_hir.call_count
    recover_hir_model(species, qc, par)
    assert qc.qc_hir.call_count == attempts


def test_remaining_hir_inconsistency_falls_back_after_recovery(tmp_path, monkeypatch, caplog):
    from kinbot.counting_contract import optical_counting
    from kinbot.thermochemistry import hir_evidence
    from kinbot.mess import MESS
    species, qc, par = model(tmp_path, monkeypatch)
    # The scan/projection references agree, but a supposedly successful scan
    # point has no geometry. This previously survived recovery and stopped MESS.
    species.hir.hir_geoms[0][4] = None
    qc.qc_hir.side_effect = RuntimeError('replacement unavailable')
    assert recover_hir_model(species, qc, par)
    assert species.reduced_freqs == species.freq
    assert species.hir.recovery_observations[0]['hir_status'] == [[0] * 12]
    assert not species.hir.is_valid_rotor(0)
    from kinbot.optical_harmonic import selected_midpoint_diagnostic
    assert selected_midpoint_diagnostic(species, ())['status'] == 'complete'
    result = optical_counting(species, hir_evidence(species))
    assert result['remaining_multiplier'] == 1.
    assert 'omitting hindered rotors' in caplog.text
    writer = MESS(dict(par, pes=1, freq_uq_ref=100., freq_uq_max_exp=2.), species)
    assert 'kept as a harmonic oscillator' in writer.make_rotors(species, 1.)
    block = writer.write_well(species, 0., 1., 0)
    assert 'Frequencies[1/cm]             12' in block
    assert 'Rotor Hindered' not in block and 'Rotor Free' not in block
    assert 'kept as a harmonic oscillator' in block
    assert qc.qc_hir.call_count == 12
    recover_hir_model(species, qc, par)
    assert qc.qc_hir.call_count == 12


def test_missing_observation_is_recovered_before_harmonic_fallback(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.hir.hir_geoms[0][4] = None
    assert recover_hir_model(species, qc, par)
    assert species.hir.is_valid_rotor(0)
    assert len(species.reduced_freqs) == len(species.freq) - 1
    assert qc.qc_hir.call_count == 12
    assert species.hir.recovery_observations[0]['hir_geoms'][0][4] is None


def test_usable_partial_scan_remains_a_rotor(tmp_path, monkeypatch):
    species, qc, par = model(tmp_path, monkeypatch)
    species.hir.hir_status[0][4] = 1
    angles = np.arange(12) * 2. * np.pi / 12
    species.hir.fourier_fit('partial', angles, 0)
    assert not recover_hir_model(species, qc, par)
    assert species.hir.is_valid_rotor(0)
    assert len(species.reduced_freqs) == len(species.freq) - 1
    qc.qc_hir.assert_not_called()
