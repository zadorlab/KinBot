"""Outage inspection must not run dispatch, alter sources or trust checkpoints."""
import hashlib
import json
from pathlib import Path
import tarfile

import pytest

from kinbot import checkpoint_preserve as cp, run_outage_audit as audit


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def snapshot(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file() and not p.is_symlink()}


def tree(tmp_path):
    root = tmp_path/'run'
    graph = root/'alpha/123/composite'
    task = graph/'tasks/geom'
    put(graph/'workflow.json', {'tasks': [{'id': 'geom', 'kind': 'ase_optimize'}]})
    put(graph/'state.json', {'tasks': {'geom': {'status': 'submitted', 'job_id': '42'}}})
    put(task/'task.json', {'geometry_sha256': 'staged-provenance'})
    put(task/'native.runtime.json', {'scratch_host': 'node1',
        'scratch_directory': '/scratch/user/kinbot-molpro-old',
        'native_work_directory': '/scratch/user/kinbot-molpro-old/work'})
    (task/'optimization.traj').write_bytes(b'not executed or parsed by capture')
    (task/'native.wfu').write_bytes(b'checkpoint-only-metadata')
    (task/'native.out').write_text('old\n'*30000+'last line\n')
    (root/'kinbot.db').write_bytes(b'database-not-opened')
    (root/'sample.pkl').write_bytes(b'never unpickled')
    return root, task


def test_capture_is_read_only_scoped_and_repeatable(tmp_path, monkeypatch):
    root, task = tree(tmp_path)
    outside = tmp_path/'secret'; outside.write_text('not captured')
    (root/'escape.json').symlink_to(outside)
    (root/'escape_dir').symlink_to(tmp_path/'other', target_is_directory=True)
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda root: [{'pid': 12, 'state': 'T'}])
    before = snapshot(root)
    first = audit.capture(root, tmp_path/'export', query_scheduler=False)
    second = audit.capture(root, tmp_path/'export', query_scheduler=False)
    assert first['directory'] != second['directory']
    assert snapshot(root) == before
    with tarfile.open(first['archive']) as archive:
        report = json.load(archive.extractfile('report.json'))
        files = {f['path']: f for f in report['files']}
        assert all(not n.endswith(('native.wfu', 'kinbot.db', 'sample.pkl', 'escape.json')) for n in archive.getnames())
        out = files[str((task/'native.out').relative_to(root))]
        assert out['offset'] > 0 and out['kind'] == 'tail'
        for item in files.values():
            if 'member' in item:
                data = archive.extractfile(item['member']).read()
                assert hashlib.sha256(data).hexdigest() == item['sha256']
        assert report['graphs'][0]['tasks'][0]['acceptance'] == 'not_revalidated_by_capture'
        assert report['scratch_locations'][0]['existence'] == 'not_checked_on_node'
    assert first['checkpoint_backup'] is False


def test_scheduler_unknown_is_not_absent(tmp_path, monkeypatch):
    def failed(argv):
        return {'argv': argv, 'status': 'unknown', 'error': 'slurmctld unavailable'}
    monkeypatch.setattr(audit, 'command', failed)
    result = audit.scheduler(tmp_path, ['42'], '2026-10-01')
    assert result['queue']['status'] == 'unknown'
    assert result['accounting'][0]['status'] == 'unknown'
    assert '--duplicates' in result['accounting'][0]['argv']


def test_requeued_active_job_is_reported_without_recovery(tmp_path, monkeypatch):
    root, task = tree(tmp_path)
    seen = []
    def command(argv):
        seen.append(argv)
        text = f'42|RUNNING|node1|04:00|{task}\n' if argv[0] == 'squeue' else 'Restarts=2'
        return {'argv': argv, 'status': 'observed', 'returncode': 0, 'stdout': text, 'stderr': ''}
    monkeypatch.setattr(audit, 'command', command)
    result = audit.capture(root, tmp_path/'exports', jobs=['42'])
    report = json.loads((Path(result['directory'])/'report.json').read_text())
    assert report['graphs'][0]['tasks'][0]['observed_allocations'][0]['job_id'] == '42'
    assert {c[0] for c in seen} == {'squeue', 'sacct', 'scontrol'}
    assert not result['restart_applied']


def test_budget_and_invalid_json_are_explicit(tmp_path):
    root, _ = tree(tmp_path)
    (root/'broken.json').write_text('{')
    result = audit.capture(root, tmp_path/'out', query_scheduler=False, logs_mb=0)
    report = json.loads((Path(result['directory'])/'report.json').read_text())
    assert report['omitted_content_count'] >= 1
    assert any(x.get('path') == 'broken.json' for x in report['issues'])


@pytest.mark.parametrize('bad', ['1;scancel', '-1', '42,43'])
def test_bad_job_ids_rejected_before_any_io(tmp_path, bad):
    root = tmp_path/'run'; root.mkdir()
    with pytest.raises(ValueError):
        audit.capture(root, tmp_path/'out', jobs=[bad])
    assert not (tmp_path/'out').exists()


def test_output_inside_run_refused(tmp_path):
    root, _ = tree(tmp_path)
    with pytest.raises(ValueError):
        audit.capture(root, root/'export', query_scheduler=False)


def scratch_tree(tmp_path, monkeypatch):
    monkeypatch.setattr(cp.socket, 'gethostname', lambda: 'node1')
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda root: [{'pid': 2, 'state': 'T'}])
    run = tmp_path/'run'; run.mkdir()
    scratch = tmp_path/'scratch/kinbot-mrcc-old/work'; scratch.mkdir(parents=True)
    (scratch/'fort.16').write_bytes(b'amplitudes')
    (scratch/'integrals').write_bytes(b'not-a-checkpoint')
    (scratch/'job.out').write_text('native evidence')
    (scratch/'job.wfu').write_bytes(b'molpro records')
    secret = tmp_path/'secret'; secret.write_text('do not copy')
    (scratch/'escape.chk').symlink_to(secret)
    other = scratch.parent.parent/'kinbot-molpro-orphan'
    other.mkdir(); (other/'previous.chk').write_bytes(b'unknown attempt')
    report = tmp_path/'report.json'
    put(report, {'root': str(run), 'scratch_locations': [dict(source='runtime.json',
        scratch_host='node1.site', scratch_directory=str(scratch.parent))]})
    return run, scratch.parent, report


def test_checkpoint_copy_hashes_preserves_sources_and_marks_unattributed(tmp_path, monkeypatch):
    run, scratch, report = scratch_tree(tmp_path, monkeypatch)
    before = snapshot(scratch.parent)
    result = cp.preserve(report, tmp_path/'preserved', include_unattributed=True)
    saved = json.loads(Path(result['receipt']).read_text())
    assert snapshot(scratch.parent) == before
    assert len(saved['unreferenced_scratch']) == 1
    assert {Path(f['source']).name for f in saved['files']} == {'fort.16', 'job.out', 'job.wfu', 'previous.chk'}
    for item in saved['files']:
        assert item['verified_copy'] and item['source_stable']
        assert cp.file_hash(item['target']) == item['sha256']
    assert next(f for f in saved['files'] if f['source'].endswith('previous.chk'))['attribution'] == 'unknown'
    assert not result['cancellation_authorized_by_receipt']
    assert not list(run.iterdir())


def test_checkpoint_budget_does_not_claim_complete_backup(tmp_path, monkeypatch):
    _, _, report = scratch_tree(tmp_path, monkeypatch)
    result = cp.preserve(report, tmp_path/'preserved', max_gib=1e-12)
    saved = json.loads(Path(result['receipt']).read_text())
    assert result['issues'] and result['files'] == 0
    assert saved['status'] == 'review_required'


def test_default_leaves_unattributed_and_wrong_host_uncopied(tmp_path, monkeypatch):
    _, _, report = scratch_tree(tmp_path, monkeypatch)
    result = cp.preserve(report, tmp_path/'preserved')
    saved = json.loads(Path(result['receipt']).read_text())
    assert all(f['attribution'] == 'recorded' for f in saved['files'])
    monkeypatch.setattr(cp.socket, 'gethostname', lambda: 'node2')
    result = cp.preserve(report, tmp_path/'preserved')
    assert result['files'] == 0


@pytest.mark.parametrize('state', ['RUNNING', 'STOPPED'])
def test_copy_orchestration_preserves_frozen_allocations_and_repeats(tmp_path, monkeypatch, state):
    root = tmp_path/'run'; root.mkdir()
    report, freeze = tmp_path/'report.json', tmp_path/'freeze.json'
    put(report, {'root': str(root), 'scratch_locations': []})
    put(freeze, {'root': str(root), 'actions': [dict(action='freeze_running', returncode=0, jobs=['41', '42'])]})
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda root: [{'state': 'T'}])
    jobs = [dict(job_id=j, state=state, reason=n, directory=str(root/'work'))
            for j,n in [('41','node1'), ('42','node1'), ('43','node2')]]
    monkeypatch.setattr('kinbot.run_cleanup._jobs', lambda root: jobs)
    calls = []
    def execute(argv, **kwargs):
        from types import SimpleNamespace
        calls.append(argv)
        assert argv[0] == 'srun' and '--jobid=41' in argv
        assert '--overlap' in argv and '--include-unattributed' in argv
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(cp.subprocess, 'run', execute)
    first = cp.allocations(report, freeze, tmp_path/'backup')
    second = cp.allocations(report, freeze, tmp_path/'backup')
    assert first['receipt'] != second['receipt'] and len(calls) == 2
    assert first['omitted_jobs'] == 1 and first['steps'] == 1
    assert first['cancels_jobs'] is False and len(jobs) == 3


@pytest.mark.parametrize('problem', ['running_owner', 'scheduler', 'different_root'])
def test_copy_orchestration_requires_frozen_scoped_ownership(tmp_path, monkeypatch, problem):
    root = tmp_path/'run'; root.mkdir()
    report, freeze = tmp_path/'report.json', tmp_path/'freeze.json'
    put(report, {'root': str(root)})
    put(freeze, {'root': str(root if problem != 'different_root' else tmp_path), 'actions': []})
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda root: [{'state': 'S' if problem == 'running_owner' else 'T'}])
    def unavailable(root):
        raise RuntimeError('scheduler unavailable')
    monkeypatch.setattr('kinbot.run_cleanup._jobs', unavailable)
    def forbidden(*a, **kw):
        pytest.fail('Should not start a copy step')
    monkeypatch.setattr(cp.subprocess, 'run', forbidden)
    with pytest.raises((RuntimeError, ValueError)):
        cp.allocations(report, freeze, tmp_path/'backup')
