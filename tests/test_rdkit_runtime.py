"""Required stereo runtime and current-format restarts, without QC jobs."""
import json
import logging
import os
from pathlib import Path
import subprocess
import sys

import pytest
from rdkit import Chem, rdBase

from kinbot import rdkit_config
from kinbot.run_format import ensure_current_run


def test_required_rdkit_and_explicit_perception_in_a_fresh_process():
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    root = str(Path(__file__).resolve().parents[1])
    environment['PYTHONPATH'] = root + os.pathsep + environment.get('PYTHONPATH', '')
    missing = subprocess.run([sys.executable, '-B', '-c',
        "import sys; sys.modules['rdkit'] = None; import kinbot"],
        env=environment, text=True, capture_output=True)
    assert missing.returncode != 0
    assert 'rdkit' in missing.stderr.lower()
    configured = subprocess.run([sys.executable, '-B', '-c',
        'from rdkit import Chem; '
        'Chem.SetUseLegacyStereoPerception(False); '
        'Chem.SetAllowNontetrahedralChirality(False); '
        'import kinbot; '
        'assert Chem.GetUseLegacyStereoPerception(); '
        'assert Chem.GetAllowNontetrahedralChirality()'],
        env=environment, text=True, capture_output=True)
    assert configured.returncode == 0, configured.stderr


def test_minimum_version_is_numeric_and_runtime_is_logged(monkeypatch, caplog):
    installed = rdBase.rdkitVersion
    monkeypatch.setattr(rdBase, 'rdkitVersion', '2025.09.3')
    rdkit_config.configure_rdkit()
    monkeypatch.setattr(rdBase, 'rdkitVersion', '2025.09.2')
    with pytest.raises((ImportError, RuntimeError), match='2025'):
        rdkit_config.configure_rdkit()
    monkeypatch.setattr(rdBase, 'rdkitVersion', installed)
    with caplog.at_level(logging.INFO):
        rdkit_config.log_rdkit(logging.getLogger('KinBot'))
    assert installed in caplog.text
    assert Chem.GetUseLegacyStereoPerception()
    assert Chem.GetAllowNontetrahedralChirality()


def snapshot(directory):
    return {str(path.relative_to(directory)): path.read_bytes()
            for path in directory.rglob('*') if path.is_file()}


def test_fresh_input_and_current_restart_share_one_unchanged_marker(tmp_path):
    (tmp_path / 'input.json').write_text('{"smiles": "CCO"}')
    (tmp_path / 'input.xyz').write_text('1\nstarting geometry\nH 0 0 0\n')
    expected = ensure_current_run(tmp_path, create=True)
    assert expected['schema'] == 'kinbot.calculation.v2'
    assert expected['rdkit_version'] == rdBase.rdkitVersion
    assert expected['use_legacy_stereo_perception'] is True
    assert expected['allow_nontetrahedral_chirality'] is True
    (tmp_path / 'kinbot.db').write_bytes(b'current result placeholder')
    before = snapshot(tmp_path)
    assert ensure_current_run(tmp_path) == expected
    assert ensure_current_run(tmp_path, create=True) == expected
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize('artifact', [
    'kinbot.db', 'kinbot.log', 'summary_123.out', 'me/mess_0000.inp',
    'molpro/123.out', 'orca/123.out',
    '123-s' + 'a'*16 + '_well.out', '123-s' + 'a'*64 + '_well.out',
])
def test_old_calculations_are_rejected_without_changing_files(tmp_path, artifact):
    path = tmp_path / artifact
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b'old calculation\n')
    before = snapshot(tmp_path)
    for create in (False, True):
        with pytest.raises(ValueError):
            ensure_current_run(tmp_path, create=create)
        assert snapshot(tmp_path) == before


@pytest.mark.parametrize(('field', 'value'), [
    ('schema', 'kinbot.calculation.v1'), ('rdkit_version', '2026.03.4'),
    ('use_legacy_stereo_perception', False), ('allow_nontetrahedral_chirality', False),
])
def test_changed_runtime_or_format_cannot_reuse_results(tmp_path, field, value):
    ensure_current_run(tmp_path, create=True)
    marker = tmp_path / '.kinbot_run.json'
    state = json.loads(marker.read_text())
    state[field] = value
    marker.write_text(json.dumps(state))
    (tmp_path / 'kinbot.db').write_bytes(b'original calculation')
    before = snapshot(tmp_path)
    for create in (False, True):
        with pytest.raises(ValueError):
            ensure_current_run(tmp_path, create=create)
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize('program', ['kinbot.kb', 'kinbot.pes'])
def test_command_rejects_old_results_before_logging_or_qc(tmp_path, monkeypatch, program):
    import importlib
    module = importlib.import_module(program)
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'input.json').write_text('{"smiles": "CCO", "barrier_threshold": 100}')
    (tmp_path / 'kinbot.log').write_text('original log\n')
    (tmp_path / 'kinbot.db').write_bytes(b'old calculation')
    before = snapshot(tmp_path)
    monkeypatch.setattr(sys, 'argv', [program, 'input.json'] +
                        (['no-kinbot'] if program.endswith('pes') else []))
    def forbidden(*args, **kwargs):
        raise AssertionError('old calculation reached logging or QC')
    monkeypatch.setattr(module, 'config_log', forbidden)
    if hasattr(module, 'QuantumChemistry'):
        monkeypatch.setattr(module, 'QuantumChemistry', forbidden)
    with pytest.raises(ValueError):
        module.main()
    assert snapshot(tmp_path) == before


def test_unsupported_initial_reactant_stops_before_qc(tmp_path, monkeypatch):
    from kinbot import kb
    from kinbot.stereo_identity import UnsupportedStereochemistry
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'input.json').write_text('{"smiles": "CCO", "barrier_threshold": 100}')
    monkeypatch.setattr(sys, 'argv', ['kinbot', 'input.json'])
    monkeypatch.setattr(kb, 'config_log', lambda *a, **k: logging.getLogger('runtime-test'))
    def unsupported(*args, **kwargs):
        raise UnsupportedStereochemistry('unsupported initial reactant fixture')
    def forbidden(*args, **kwargs):
        raise AssertionError('unsupported initial reactant reached QC')
    monkeypatch.setattr(kb, 'require_supported_identity', unsupported)
    monkeypatch.setattr(kb, 'QuantumChemistry', forbidden)
    with pytest.raises(UnsupportedStereochemistry, match='unsupported initial reactant'):
        kb.main()


@pytest.mark.parametrize('old_worker', [True, False])
def test_pes_checks_worker_format_before_reading_energies(tmp_path, monkeypatch, old_worker):
    from kinbot import pes
    monkeypatch.chdir(tmp_path)
    ensure_current_run(create=True)
    worker = tmp_path / '123'
    worker.mkdir()
    if not old_worker:
        ensure_current_run(worker, create=True)
        marker = worker / '.kinbot_run.json'
        state = json.loads(marker.read_text())
        state['rdkit_version'] = '2026.03.4'
        marker.write_text(json.dumps(state))
    (worker / 'kinbot.db').write_bytes(b'worker calculation')
    before = snapshot(tmp_path)
    def forbidden(*args, **kwargs):
        raise AssertionError('incompatible worker reached energy loading')
    monkeypatch.setattr(pes, 'get_energy', forbidden)
    with pytest.raises(ValueError):
        pes.postprocess({}, ['123'], 'all', [], 30.)
    assert snapshot(tmp_path) == before


def test_pes_submission_rejects_old_worker_before_deleting_or_writing(tmp_path, monkeypatch):
    from kinbot import pes
    monkeypatch.chdir(tmp_path)
    ensure_current_run(create=True)
    worker = tmp_path / '123'
    worker.mkdir()
    for name in ('summary_123.out', 'kinbot_monitor.out', 'kinbot.out', 'kinbot.err'):
        (worker / name).write_text('original worker result\n')
    before = snapshot(tmp_path)
    def forbidden(*args, **kwargs):
        raise AssertionError('old worker reached file removal or job submission')
    monkeypatch.setattr(pes.os, 'system', forbidden)
    monkeypatch.setattr(pes.subprocess, 'Popen', forbidden)
    with pytest.raises(ValueError):
        pes.submit_job('123', {})
    assert snapshot(tmp_path) == before


def test_species_directory_alias_cannot_authorize_a_restart(tmp_path):
    target = tmp_path / 'original'
    ensure_current_run(target, create=True)
    (target / 'kinbot.db').write_bytes(b'original calculation')
    alias = tmp_path / '123'
    alias.symlink_to(target, target_is_directory=True)
    before = snapshot(target)
    for create in (False, True):
        with pytest.raises(ValueError):
            ensure_current_run(alias, create=create)
    assert alias.is_symlink()
    assert snapshot(target) == before
    with pytest.raises(ValueError):
        ensure_current_run(tmp_path, create=True)
    assert not (tmp_path / '.kinbot_run.json').exists()
