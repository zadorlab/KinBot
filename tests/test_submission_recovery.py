"""Real pool admission/dispatch with simulated Slurm; no native QC runs."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from kinbot import run_policy
from kinbot.anl import dispatch, submission_recovery
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_scratch_recovery import failed


@pytest.fixture
def pool(tmp_path, monkeypatch):
    root = tmp_path/'pool'; root.mkdir()
    run_policy.update(root, max_jobs=3, fair_submission=False)
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task('bad', 'rhf\n', geometry_from='initial')]
    graph = dispatch.prepare(spec, root/'graph')
    jobs, calls = {}, []
    def slurm(command, **kwargs):
        if command[0] == 'sbatch':
            calls.append((command, kwargs))
            job = str(100+len(calls)); jobs[job] = str(kwargs['cwd'])
            return SimpleNamespace(returncode=0, stdout=job, stderr='')
        if command[0] == 'squeue':
            form = command[command.index('-o')+1]
            if form == '%i':
                selected = jobs if '-j' not in command else {
                    k:v for k,v in jobs.items() if k == command[command.index('-j')+1]}
                output = '\n'.join(selected)
            elif form == '%i|%D|%Z':
                output = '\n'.join(f'{k}|1|{v}' for k,v in jobs.items())
            elif form == '%i|%T|%R|%Z':
                output = '\n'.join(f'{k}|RUNNING|node8|{v}' for k,v in jobs.items())
            else:
                pytest.fail(str(command))
            return SimpleNamespace(returncode=0, stdout=output, stderr='')
        if command[:3] == ['scontrol', 'show', 'job']:
            return SimpleNamespace(returncode=1, stdout='', stderr='slurm_load_jobs error: Invalid job id specified')
        pytest.fail(str(command))
    monkeypatch.setattr(subprocess, 'run', slurm)
    monkeypatch.setattr(run_policy, 'slurm_job_record', lambda j: {})
    monkeypatch.setenv('KINBOT_QC_SCRATCH_RETRIES', '2')
    return SimpleNamespace(root=root, graph=graph, jobs=jobs, calls=calls, slurm=slurm)


def archived_failure(pool, *, legacy_submitting=False):
    dispatch.advance(pool.graph, submit=True)
    failed(pool.graph, job='101', node='badnode', legacy=False)
    pool.jobs.clear()  # Completed before the next admission poll; Slurm purged it.
    state = dispatch.advance(pool.graph)
    assert state['tasks']['bad']['attempt'] == 2
    if legacy_submitting:
        state['tasks']['bad']['status'] = 'submitting'
        atomic_json(pool.graph/'state.json', state)
    return state


def test_archived_terminal_execution_releases_known_claim_once(pool):
    state = archived_failure(pool)
    archived = Path(state['tasks']['bad']['previous_attempt'])
    preserved = {p:p.read_bytes() for p in archived.rglob('*') if p.is_file()}
    state = dispatch.advance(pool.graph, submit=True)
    assert state['tasks']['bad']['job_id'] == '102'
    assert state['tasks']['bad']['attempt'] == 2
    proofs = list((pool.root/'.kinbot_policy_reconciled').glob('*.json'))
    assert len(proofs) == 1
    proof = json.loads(proofs[0].read_text())
    assert proof['job_id'] == '101' and proof['execution_status'] == 'failed'
    assert proof['execution_sha256'] == file_sha256(archived/'execution.json')
    dispatch.advance(pool.graph, submit=True)
    assert len(pool.calls) == 2
    assert all(p.read_bytes() == b for p,b in preserved.items())


def test_current_terminal_execution_can_reconcile_without_retry(pool):
    dispatch.advance(pool.graph, submit=True)
    failed(pool.graph, job='101', legacy=False)
    pool.jobs.clear()
    ledger = pool.root/'.kinbot_policy_admissions.json'
    claims = json.loads(ledger.read_text())
    assert run_policy._reconcile(pool.root, {}, claims) == 0
    assert not claims
    assert json.loads(ledger.read_text())  # A read-only inspection changes nothing.
    assert not (pool.root/'.kinbot_policy_reconciled').exists()


def test_legacy_pre_admission_state_uses_old_claim_and_archived_execution(pool):
    state = archived_failure(pool, legacy_submitting=True)
    original_spec = (pool.graph/'workflow.json').read_bytes()
    # A pool poll must not remove the proof before the owner repairs its state.
    claims = json.loads((pool.root/'.kinbot_policy_admissions.json').read_text())
    assert run_policy._reconcile(pool.root, {}, claims) == 1
    repaired = submission_recovery.reconcile_pre_admission(pool.graph)
    assert repaired['tasks']['bad']['status'] == 'staged'
    assert repaired['tasks']['bad']['attempt'] == 2
    assert repaired['scratch_node_retries'] == state['scratch_node_retries']
    assert submission_recovery.reconcile_pre_admission(pool.graph) == repaired
    result = dispatch.advance(pool.graph, submit=True)
    assert result['tasks']['bad']['job_id'] == '102' and len(pool.calls) == 2
    assert (pool.graph/'workflow.json').read_bytes() == original_spec


@pytest.mark.parametrize('fault', ['unknown_claim','live_old','live_new','new_receipt','new_claim',
                                  'worker_artifact','wrong_execution','wrong_archive','receipt_changed',
                                  'new_submission_marker','scheduler_offline'])
def test_legacy_repair_never_waives_uncertainty(pool, monkeypatch, fault):
    state = archived_failure(pool, legacy_submitting=True)
    entry = state['tasks']['bad']; directory = pool.graph/'tasks/bad'
    ledger = pool.root/'.kinbot_policy_admissions.json'
    claims = json.loads(ledger.read_text()); token = next(iter(claims)); claim = claims[token]
    if fault == 'unknown_claim':
        claim.pop('job_id'); claim['status'] = 'submitting'; atomic_json(ledger, claims)
    elif fault in ('live_old','live_new'):
        pool.jobs['101' if fault == 'live_old' else '999'] = str(directory)
    elif fault == 'new_receipt':
        atomic_json(pool.root/'.kinbot_policy_submissions/999.json', dict(claim, job_id='999'))
    elif fault == 'new_claim':
        claims['new'] = dict(claim, job_id='999'); atomic_json(ledger, claims)
    elif fault == 'worker_artifact':
        (directory/'slurm.stdout').write_text('')
    elif fault == 'wrong_execution':
        path = Path(entry['previous_attempt'])/'execution.json'
        atomic_json(path, dict(json.loads(path.read_text()), slurm_job_id='999'))
    elif fault == 'wrong_archive':
        (Path(entry['previous_attempt'])/'task.json').write_text('{}')
    elif fault == 'receipt_changed':
        atomic_json(pool.root/'.kinbot_policy_submissions/101.json', dict(claim, directory='/elsewhere'))
    elif fault == 'new_submission_marker':
        entry['submission_started_utc'] = '2026-01-01T00:00:00+00:00'
        atomic_json(pool.graph/'state.json', state)
    else:
        def offline(command, **kwargs):
            if command[0] == 'scontrol':
                return SimpleNamespace(returncode=1, stdout='', stderr='controller unavailable')
            return pool.slurm(command, **kwargs)
        monkeypatch.setattr(subprocess, 'run', offline)
    before = ledger.read_bytes()
    assert submission_recovery.reconcile_pre_admission(pool.graph)['tasks']['bad']['status'] == 'submitting'
    assert ledger.read_bytes() == before and len(pool.calls) == 1
    assert not (pool.root/'.kinbot_policy_reconciled').exists()


def test_admission_rejection_does_not_mark_task_submitting(pool):
    directory = pool.graph/'tasks/bad'
    atomic_json(pool.root/'.kinbot_policy_admissions.json', {'unknown':dict(
        directory=str(directory),command=['sbatch','--parsable','job.slurm'],status='submitting')})
    with pytest.raises(RuntimeError, match='unconfirmed'):
        dispatch.advance(pool.graph, submit=True)
    entry = dispatch._load(pool.graph)[2]['tasks']['bad']
    assert entry['status'] == 'staged' and 'submission_started_utc' not in entry
    assert not pool.calls


def test_sbatch_timeout_remains_fenced(pool, monkeypatch):
    def slurm(command, **kwargs):
        if command[0] == 'sbatch':
            raise subprocess.TimeoutExpired(command, 30)
        return pool.slurm(command, **kwargs)
    monkeypatch.setattr(subprocess, 'run', slurm)
    with pytest.raises(RuntimeError, match='timed out'):
        dispatch.advance(pool.graph, submit=True)
    state = dispatch._load(pool.graph)[2]
    assert state['tasks']['bad']['status'] == 'submitting'
    assert state['tasks']['bad']['submission_started_utc']
    assert submission_recovery.reconcile_pre_admission(pool.graph) == state


def test_crash_between_state_repair_and_claim_retirement_is_idempotent(pool, monkeypatch):
    archived_failure(pool, legacy_submitting=True)
    original = submission_recovery.atomic_json
    def interrupt(path, value):
        if Path(path).name == '.kinbot_policy_admissions.json':
            raise OSError('interrupted')
        original(path, value)
    monkeypatch.setattr(submission_recovery, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='interrupted'):
        submission_recovery.reconcile_pre_admission(pool.graph)
    monkeypatch.setattr(submission_recovery, 'atomic_json', original)
    assert dispatch._load(pool.graph)[2]['tasks']['bad']['status'] == 'staged'
    assert dispatch.advance(pool.graph, submit=True)['tasks']['bad']['job_id'] == '102'
    assert len(pool.calls) == 2
