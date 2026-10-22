"""Native-diagnostic VPT2 recovery is shared by every graph owner."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot.anl import dispatch, live_control
from kinbot import run_recover, run_control
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def graph(tmp_path, *, symmetry='Symmetry=None', freq='Freq=Anharmonic', route='#p', linear=False):
    spec = dispatch_spec()
    task = deepcopy(next(t for t in spec['tasks'] if t['id'] == 'gaussian_vpt2'))
    task.update(id='custom_anharmonic')
    task.pop('depends_on')
    task['input_template'] = task['input_template'].replace('NoSymm', symmetry).replace('Freq=Anharmonic', freq).replace('#p', route)
    spec['tasks'] = [spec['tasks'][0], task, molpro_task('live', 'rhf\n', geometry_from='initial'),
                     molpro_task('done', 'rhf\n', geometry_from='initial')]
    if linear:
        spec['molecule'] = dict(symbols=['H', 'C', 'N'], positions=[[0,0,-1.06],[0,0,0],[0,0,1.15]], charge=0, multiplicity=1)
    root = dispatch.prepare(spec, tmp_path/'graph')
    _mock_execution(root, 'l2_geometry', geometry='final.xyz')
    dispatch.advance(root)
    return root


def fail(root, text=None, *, legacy_linear=False, status='failed'):
    d = root/'tasks/custom_anharmonic'
    state = dispatch._load(root)[2]
    entry = state['tasks']['custom_anharmonic']
    entry.update(status=status, job_id='42')
    if legacy_linear:
        p = d/'vpt2.com'
        p.write_text(p.read_text().replace('NoSymm', 'Symmetry=None'))
        entry['input_sha256'] = file_sha256(p)
    text = text or ('ERROR: Inconsistency found in framework group definition:\n'
                    'New: CS - Old:\nError termination via Lnk1e in l717.exe\n')
    (d/'vpt2.log').write_text(text)
    atomic_json(d/'execution.json', dict(schema=1, task_id='custom_anharmonic',
        geometry_sha256=entry['geometry_sha256'], status='failed', slurm_job_id='42',
        failure_phase='execution', error='Program exited with status 1.'))
    atomic_json(root/'state.json', state)
    return d


@pytest.mark.parametrize('group', ['CS', 'C1', 'C2V', 'TD'])
@pytest.mark.parametrize('syntax', ['Symmetry=None', 'Symm=(None)'])
def test_common_dispatch_recovers_any_task_and_native_group(tmp_path, monkeypatch, group, syntax):
    root = graph(tmp_path, symmetry=syntax, freq='Frequencies=(Anharmonic,VibRot)')
    d = fail(root, f'ERROR: Inconsistency found in framework group definition:\nNew: {group} - Old:\nError termination via Lnk1e in l717.exe\n')
    geometry = (d/'geometry.xyz').read_bytes()
    _mock_execution(root, 'done')
    state = dispatch._load(root)[2]
    state['tasks']['live'].update(status='submitted', job_id='43')
    atomic_json(root/'state.json', state)
    hashes = {p:file_sha256(p) for p in (root/'tasks/live').iterdir() if p.is_file()}
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job == '43')
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['attempt'] == 2
    assert state['tasks']['custom_anharmonic']['status'] == 'staged'
    assert f'Symmetry=(PG={group})' in (d/'vpt2.com').read_text()
    assert 'Frequencies=(Anharmonic,VibRot)' in (d/'vpt2.com').read_text()
    assert (d/'geometry.xyz').read_bytes() == geometry
    assert state['tasks']['done']['status'] == 'complete'
    assert state['tasks']['live']['job_id'] == '43'
    assert hashes == {p:file_sha256(p) for p in hashes}
    # Repeating the native failure cannot create an unbounded retry loop.
    fail(root)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert len(list((root/'attempts').glob('live_*'))) == 1


@pytest.mark.parametrize('route', ['#', '#P', '#n', '#t'])
def test_compound_frequency_without_symmetry_gets_route_keyword(tmp_path, monkeypatch, route):
    root = graph(tmp_path, symmetry='', freq='Freq=(Anharmonic,VibRot)', route=route)
    d = fail(root, 'ERROR: Inconsistency found in framework group definition:\nNew: CS - Old: C1\nError termination via Lnk1e\n')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    dispatch.advance(root)
    assert 'Symmetry=None' in (d/'vpt2.com').read_text()
    assert 'Freq=(Anharmonic,VibRot)' in (d/'vpt2.com').read_text()


@pytest.mark.parametrize('stage', ['failed', 'staged'])
def test_linear_native_fallback_reaches_common_dispatcher(tmp_path, monkeypatch, stage):
    root = graph(tmp_path, linear=True, freq='Freq=(Anharmonic,VibRot)')
    d = fail(root, 'Molecule is non-linear but NTrRo=5\nError termination via Lnk1e in l717.exe\n', legacy_linear=True, status=stage)
    if stage == 'staged':
        (d/'execution.json').unlink()
        state = dispatch._load(root)[2]
        state['tasks']['custom_anharmonic'].pop('job_id')
        atomic_json(root/'state.json', state)
    geometry = (d/'geometry.xyz').read_bytes()
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'NoSymm' in (d/'vpt2.com').read_text()
    assert (d/'geometry.xyz').read_bytes() == geometry
    assert 'gaussian_vpt2_linear_framework_v1' in state['native_interface_repairs']['custom_anharmonic']


@pytest.mark.parametrize('diagnostic', ['Error: segmentation violation', 'SCF convergence failure', 'Molecule is non-linear but NTrRo=5'])
def test_unmatched_or_nonlinear_failure_is_not_relabelled(tmp_path, monkeypatch, diagnostic):
    root = graph(tmp_path)
    fail(root, diagnostic+'\nError termination via Lnk1e\n')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert not (root/'.live_control_transaction.json').exists()


def test_recovery_resumes_interrupted_manifest_transaction(tmp_path, monkeypatch):
    root = graph(tmp_path); fail(root)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    original = run_recover.atomic_json
    def interrupt(path, value):
        if Path(path) == root/'state.json':
            raise OSError('interrupt')
        return original(path, value)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='interrupt'):
        dispatch.advance(root)
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['attempt'] == 2
    assert 'Symmetry=(PG=CS)' in (root/'tasks/custom_anharmonic/vpt2.com').read_text()
    assert len(list((root/'attempts').glob('live_*'))) == 1


@pytest.mark.parametrize('block', ['active', 'legacy', 'unknown'])
def test_allocation_fences_prevent_recovery_race(tmp_path, monkeypatch, block):
    root = graph(tmp_path); d = fail(root)
    before = (d/'vpt2.com').read_bytes()
    if block == 'legacy':
        state = dispatch._load(root)[2]
        state['tasks']['live'].update(status='submitted', job_id='43')
        atomic_json(root/'state.json', state)
        p = root/'tasks/live/task.json'; record=json.loads(p.read_text())
        record.pop('live_control_schema'); atomic_json(p, record)
        state['tasks']['live']['task_sha256']=file_sha256(p)
        atomic_json(root/'state.json', state)
    def active(job):
        if block == 'unknown':
            raise dispatch.SchedulerUnavailableError('queue unavailable')
        return block == 'active' or job == '43'
    monkeypatch.setattr(dispatch, '_job_active', active)
    if block == 'unknown':
        with pytest.raises(dispatch.SchedulerUnavailableError):
            dispatch.advance(root)
    else:
        dispatch.advance(root)
    assert (d/'vpt2.com').read_bytes() == before
    assert not (root/'.live_control_transaction.json').exists()


def test_reconcile_is_previewed_non_submitting_and_respects_live_owner(tmp_path, monkeypatch):
    import fcntl
    root = graph(tmp_path); fail(root)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    before = (root/'state.json').read_bytes()
    assert run_control.reconcile(tmp_path)['graphs'][0]['status'] == 'reviewed'
    assert (root/'state.json').read_bytes() == before
    with (root/'drive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert run_control.reconcile(tmp_path, apply=True)['graphs'][0]['status'] == 'owned_by_controller'
        assert (root/'state.json').read_bytes() == before
    def no_submit(*a, **kw):
        raise AssertionError('No scheduler mutation permitted')
    monkeypatch.setattr(dispatch.subprocess, 'run', no_submit)
    report = run_control.reconcile(tmp_path, ['graph'], apply=True)['graphs'][0]
    assert report['status'] == 'reconciled'
    assert report['changes']['custom_anharmonic']['after'] == 'staged'
    assert run_control.reconcile(tmp_path, ['no_match'], apply=True)['graphs'] == []


@pytest.mark.parametrize('frequency', ['Freq=(Anharmonic,VibRot)', 'Frequencies=(VibRot,Anharmonic)', 'Freq(Anharmonic)', 'Freq=Anharm'])
def test_output_reader_accepts_same_compound_options_as_input(frequency):
    from tests.test_anl_results import _GAUSSIAN_VPT2
    from kinbot.anl.results import parse_gaussian_vpt2
    expected = parse_gaussian_vpt2(_GAUSSIAN_VPT2, method='B3LYP', basis='cc-pVTZ')
    result = parse_gaussian_vpt2(_GAUSSIAN_VPT2.replace('Freq=Anharmonic', frequency), method='B3LYP', basis='cc-pVTZ')
    assert result == expected
    with pytest.raises(ValueError, match='requested Freq=Anharmonic'):
        parse_gaussian_vpt2(_GAUSSIAN_VPT2.replace('Freq=Anharmonic', 'Freq=Harmonic'), method='B3LYP', basis='cc-pVTZ')


def test_general_framework_recovery_ladder_is_bounded(tmp_path, monkeypatch):
    root = graph(tmp_path, symmetry='Symmetry=PG=C1')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    d = fail(root, 'ERROR: Inconsistency found in framework group definition:\nNew: TD - Old: T\nError termination via Lnk1e\n')
    dispatch.advance(root)
    assert 'Symmetry=None' in (d/'vpt2.com').read_text()
    fail(root)
    dispatch.advance(root)
    assert 'Symmetry=(PG=CS)' in (d/'vpt2.com').read_text()
    fail(root)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'failed'
    assert state['tasks']['custom_anharmonic']['attempt'] == 3


def test_unrecognized_user_symmetry_requires_manual_review_without_abandoning_siblings(tmp_path, monkeypatch):
    root = graph(tmp_path, symmetry='Symmetry=(Loose,Follow)')
    d = fail(root)
    _mock_execution(root, 'done')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'failed'
    assert state['tasks']['done']['status'] == 'complete'
    assert not (root/'.live_control_transaction.json').exists()


def test_reconcile_isolates_one_graph_error(tmp_path, monkeypatch):
    bad = graph(tmp_path/'bad'); fail(bad)
    good = graph(tmp_path/'good'); fail(good)
    original = dispatch.advance
    def advance(path, **kwargs):
        if path == bad:
            raise RuntimeError('bad provenance')
        assert kwargs == {'submit': False}
        return original(path, **kwargs)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    monkeypatch.setattr(dispatch, 'advance', advance)
    records = run_control.reconcile(tmp_path, apply=True)['graphs']
    assert [(r['graph'], r['status']) for r in records] == [('bad/graph', 'blocked'), ('good/graph', 'reconciled')]


def test_mismatched_allocation_cannot_trigger_recovery(tmp_path, monkeypatch):
    root = graph(tmp_path); d = fail(root)
    p = d/'execution.json'; outcome = json.loads(p.read_text())
    outcome['slurm_job_id'] = '99'; atomic_json(p, outcome)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    with pytest.raises(RuntimeError, match='different allocation'):
        dispatch.advance(root)
    assert not (root/'.live_control_transaction.json').exists()
