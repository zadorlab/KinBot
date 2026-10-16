"""Numerical trial wiring/provenance only; these tests do not run Molpro."""
import json
from pathlib import Path
import re
from types import SimpleNamespace

import pytest

from kinbot.execution import file_sha256
from kinbot.molpro import _vrc_multireference_method, validate_mr_correlation_options
from kinbot import molpro_reconnaissance as recon
from test_mrci_simultaneous import model
from test_molpro_confirmation_recovery import fixture
from test_molpro_confirmation_package import PACKAGE


@pytest.mark.parametrize('bad', [{}, [], {'pair': 0.1}, {'pairs': True},
    {'singles': '0.1'}, {'internal': -0.1}, {'pairs': float('nan')},
    {'pairs': float('inf')}])
def test_reject_invalid_preconditioning(bad):
    with pytest.raises(ValueError, match='denominator_shifts'):
        validate_mr_correlation_options({'denominator_shifts': bad})


def test_method_specific_control_and_unchanged_default():
    state = model()
    plain, _ = _vrc_multireference_method('mrci-f12+q(5,4)', 25, 1, state_model=state)
    assert not re.search(r'shift[psi]=', plain)
    state['correlation_options']['denominator_shifts'] = {'pairs': 0.1, 'singles': 0.1, 'internal': 0.1}
    trial, _ = _vrc_multireference_method('mrci-f12+q(5,4)', 25, 1, state_model=state)
    assert trial.replace(',shiftp=0.1,shifts=0.1,shifti=0.1', '') == plain
    state['correlation_options']['state_policy'] = 'independent'
    with pytest.raises(ValueError, match='not a CASPT2 level shift'):
        _vrc_multireference_method('caspt2-f12(5,4)', 25, 1, state_model=state)
    # Partial maps and zero shifts are explicit overrides; missing keys stay absent.
    state['correlation_options']['denominator_shifts'] = {'pairs': 0}
    trial, _ = _vrc_multireference_method('mrci+q(5,4)', 25, 1, state_model=state)
    assert 'shiftp=0' in trial and 'shifts=' not in trial and 'shifti=' not in trial


def stopped_multistate(tmp_path, monkeypatch):
    adapter, original, source, _, review = fixture(tmp_path, monkeypatch)
    record = adapter.prepare(original, source, PACKAGE/'multistate_resolution_plan.json', review)
    for system in record['systems']:
        p = Path(system['manifest'])
        manifest = json.loads(p.read_text())
        for job in manifest['jobs']:
            d, name = Path(job['directory']), job['name']
            if system['title'].startswith('ch3ss_ss') and job['anchor'] != 'far':
                job['status'] = 'incomplete' if job['anchor'] == 'center' else 'prepared'
                if job['anchor'] == 'center':
                    # Collection can remove terminal_failure; execution remains authoritative.
                    (d/'execution.json').write_text(json.dumps({'status': 'failed', 'returncode': 255}))
                    (d/(name+'.out')).write_text('? dsygv failed\n? The problem occurs in cif12solve\n')
                continue
            out, wfu = d/(name+'.out'), d/(name+'.wfu')
            out.write_text('Synthetic completed point fixture')
            wfu.write_bytes(b'Synthetic CAS(7,5) checkpoint fixture')
            (d/'execution.json').write_text(json.dumps({'status': 'complete', 'returncode': 0}))
            job.update(status='complete', output_sha256=file_sha256(out),
                       wfu=str(wfu), wfu_sha256=file_sha256(wfu))
        p.write_text(json.dumps(manifest))
    record['status'] = 'failed'
    (source/'confirmation.json').write_text(json.dumps(record))
    # _audit intentionally writes a manifest. Ensure the new adapter isolates
    # that write; native-audit validation itself is tested in reconnaissance tests.
    def audit_copy(path):
        assert source not in Path(path).parents
        data = json.loads(Path(path).read_text())
        data['fixture_audit_write'] = True
        Path(path).write_text(json.dumps(data))
        return data
    monkeypatch.setattr(recon, '_audit', audit_copy)
    return adapter, source, tmp_path/'solver_retry'


def test_trial_preserves_five_results_and_only_submits_two_points(tmp_path, monkeypatch):
    adapter, source, target = stopped_multistate(tmp_path, monkeypatch)
    before = {p: file_sha256(p) for p in source.rglob('*') if p.is_file()}
    plan = PACKAGE/'mrci_preconditioner_plan.json'
    record = adapter.prepare(source, target, plan, None)
    assert len(record['systems']) == 1 and record['systems'][0]['jobs'] == 2
    assert record['execution']['max_jobs'] == 1 and record['execution']['memory_mb'] == 104000
    cfg = json.loads(Path(record['systems'][0]['configuration']).read_text())
    assert cfg['probe_definition']['orbital_partition'] == {
        'source_active_electrons': 7, 'source_active_orbitals': 5,
        'close_lowest': 0, 'discard_highest': 0, 'orbital_swaps': []}
    snapshot = json.loads((target/'source_completion_snapshot.json').read_text())
    assert snapshot['production_approved'] is False
    assert sum(j['status'] == 'complete' for s in snapshot['systems'] for j in s['jobs']) == 5
    m = json.loads(Path(record['systems'][0]['manifest']).read_text())
    assert [j['anchor'] for j in m['jobs']] == ['center', 'tilted']
    assert not m['jobs'][0].get('orbital_predecessor')
    assert m['jobs'][1]['orbital_predecessor'] == m['jobs'][0]['name']
    seed = next(iter(cfg['initial_records'].values()))
    assert file_sha256(seed['path']) == seed['sha256']
    assert all(file_sha256(p) == digest for p, digest in before.items())
    assert adapter.prepare(source, target, plan, None) == record

    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *a, **k: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *a: None)
    commands = []
    def launch(command, **kwargs):
        native = (kwargs['cwd']/command[-1]).read_text()
        blocks = re.findall(r'\{mrci-f12[^{}]*\}', native)
        assert len(blocks) == 2
        for block in blocks:
            for token in ('shiftp=0.1,shifts=0.1,shifti=0.1', 'cabs_singles=3,noswap',
                          'state,6,1,2,3,4,5,6;', 'refstate,6,1,2,3,4,5,6;',
                          'maxdav=20', 'maxvi=60', 'thrvar=1.d-10', 'cluster=0'):
                assert token in block
        assert len(re.findall(r'^reconc\d+=', native, re.M)) == 12
        for slot in range(1,7):
            assert native.count(f'KBVRCSEL=energd({slot})') == 2
        commands.append(command)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', launch)
    for job in m['jobs']:
        job['retry_solver_options'] = m['solver_options']
        if job.get('orbital_predecessor'):
            (Path(job['directory'])/(job['name']+'.wfu')).write_bytes(b'Synthetic predecessor')
        recon._execute(job)
    assert len(commands) == 2


@pytest.mark.parametrize('failure', ['completion', 'checkpoint', 'live', 'different_error'])
def test_trial_refuses_invalid_source_without_creating_destination(tmp_path, monkeypatch, failure):
    adapter, source, target = stopped_multistate(tmp_path, monkeypatch)
    original_audit = recon._audit
    if failure in ('completion', 'checkpoint'):
        def bad_audit(path):
            result = original_audit(path)
            job = next(j for j in result['jobs'] if j['status'] == 'complete')
            job['status' if failure == 'completion' else 'wfu_sha256'] = 'changed'
            return result
        monkeypatch.setattr(recon, '_audit', bad_audit)
    elif failure == 'live':
        monkeypatch.setattr(adapter.run_cleanup, 'inspect', lambda p: {'jobs': ['live'], 'processes': []})
    else:
        (source/'ch3ss_ss_atz_mrci_simultaneous_cabs3/atz__center/atz__center.out').write_text('Different error')
    with pytest.raises((ValueError, RuntimeError)):
        adapter.prepare(source, target, PACKAGE/'mrci_preconditioner_plan.json', None)
    assert not target.exists()
