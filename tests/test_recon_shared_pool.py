"""Review and production allocations share admission, including uncertain starts."""
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from kinbot import requested_orbital_recon as review, molpro_reconnaissance as recon
from kinbot import run_policy as pool
from kinbot.execution import atomic_json, file_sha256
from test_requested_orbital_recon import request_tree, RECIPE


@pytest.fixture
def shared(request_tree, monkeypatch):
    request, case, _ = request_tree
    root = case.parent.parent
    policy = root/'.kinbot_resource_policy.json'
    atomic_json(policy, dict(schema=1, max_jobs=23, excluded_nodes=['badnode'],
                            fair_submission=False, revision=0))
    before = policy.read_bytes()
    output = root/'orbital_recon'/'screen'
    review.prepare(request, case, RECIPE, output, shared_pool=True)
    assert policy.read_bytes() == before
    manifest = output/'native/reconnaissance.json'
    recon._audit(manifest, live_job_ids=[])
    jobs = {str(800+i): 1 for i in range(23)}
    monkeypatch.setattr(pool, '_queue', lambda root: jobs.copy())
    monkeypatch.setattr(pool, 'slurm_job_record', lambda ident: {})
    from kinbot import run_walltime
    monkeypatch.setattr(run_walltime, 'partition_walltime', lambda p: '31-00:00:00')
    return root, output, manifest, jobs


def submit(manifest):
    return recon._submit(manifest, max_jobs=3, startup_checked=True, live_job_ids=[],
                         audited_manifest=json.loads(manifest.read_text()))


def test_full_pool_has_no_false_submitting_then_admits_once(shared, monkeypatch):
    root, output, manifest, jobs = shared
    calls = []
    def launch(command, **kwargs):
        calls.append(command)
        assert json.loads(manifest.read_text())['jobs'][0]['submission_status'] == 'submitting'
        jobs['900'] = 1
        return SimpleNamespace(returncode=0, stdout='900\n', stderr='')
    monkeypatch.setattr(pool.subprocess, 'run', launch)
    before = (root/'.kinbot_resource_policy.json').read_bytes()
    state = submit(manifest)
    assert not calls and not any(j.get('submission_status') == 'submitting' for j in state['jobs'])
    jobs.pop('800')
    state = submit(manifest)
    assert state['jobs'][0]['job_id'] == '900'
    assert len(calls) == 1 and len(jobs) == 23
    assert '--exclude=badnode' in calls[0] and '--time=31-00:00:00' in calls[0]
    submit(manifest)
    assert len(calls) == 1
    assert (root/'.kinbot_resource_policy.json').read_bytes() == before
    assert (root/'.kinbot_policy_submissions/900.json').exists()


def test_timeout_keeps_both_journals_and_never_resubmits(shared, monkeypatch):
    root, _, manifest, jobs = shared
    jobs.pop('800')
    def timeout(command, **kwargs):
        raise subprocess.TimeoutExpired(command, 30)
    monkeypatch.setattr(pool.subprocess, 'run', timeout)
    with pytest.raises(RuntimeError, match='uncertain'):
        submit(manifest)
    claims = json.loads((root/'.kinbot_policy_admissions.json').read_text())
    assert len(claims) == 1 and next(iter(claims.values()))['status'] == 'submitting'
    with pytest.raises(RuntimeError, match='Interrupted Slurm submission'):
        submit(manifest)


def test_failed_exec_is_retryable_and_pause_changes_no_job(shared, monkeypatch):
    root, _, manifest, jobs = shared
    jobs.pop('800')
    def missing(*args, **kwargs):
        raise FileNotFoundError('sbatch')
    monkeypatch.setattr(pool.subprocess, 'run', missing)
    with pytest.raises(FileNotFoundError):
        submit(manifest)
    assert json.loads(manifest.read_text())['jobs'][0]['submission_status'] == 'retry_ready'
    assert json.loads((root/'.kinbot_policy_admissions.json').read_text()) == {}
    pool.update(root, max_jobs=0)
    submit(manifest)  # No new execution despite the prepared inventory.
    assert len(jobs) == 22


def test_policy_and_native_exclusions_union_without_mutating_policy(shared, monkeypatch):
    root, _, manifest, jobs = shared
    jobs.pop('800')
    state = json.loads(manifest.read_text()); state['slurm_excluded_nodes'] = ['failednode']
    atomic_json(manifest, state)
    pool.update(root, allow=['badnode'], exclude=['newbad'])
    before = (root/'.kinbot_resource_policy.json').read_bytes()
    calls = []
    def launch(command, **kwargs):
        calls.append(command); jobs['900'] = 1
        return SimpleNamespace(returncode=0, stdout='900\n', stderr='')
    monkeypatch.setattr(pool.subprocess, 'run', launch)
    submit(manifest)
    assert '--exclude=failednode,newbad' in calls[0]
    assert (root/'.kinbot_resource_policy.json').read_bytes() == before


def test_shared_review_cannot_escape_policy_or_bypass_limit(shared, monkeypatch):
    root, output, manifest, jobs = shared
    with pytest.raises(ValueError, match='live pool limit'):
        review.run(output, queue_all=True)
    state = json.loads(manifest.read_text()); state['shared_policy_file'] = '/wrong/policy.json'
    atomic_json(manifest, state)
    monkeypatch.setattr(pool.subprocess, 'run', lambda *a, **k: pytest.fail('submitted'))
    with pytest.raises(ValueError, match='policy/path changed'):
        submit(manifest)


@pytest.mark.parametrize('archived', [False, True])
def test_purged_terminal_attempt_releases_claim_with_exact_provenance(shared, monkeypatch, archived):
    from kinbot.anl import submission_recovery
    root, _, manifest, jobs = shared
    state = json.loads(manifest.read_text()); job = state['jobs'][0]
    directory = Path(job['directory']); job.update(job_id='900', submission_status='submitted')
    command = ['sbatch', '--export=ALL', '--parsable', job['name']+'.sbatch']
    claim = dict(job_id='900', status='submitted', directory=str(directory), command=command)
    receipt = root/'.kinbot_policy_submissions/900.json'; atomic_json(receipt, claim)
    attempt = directory
    if archived:
        attempt = directory/'attempts/1'; attempt.mkdir(parents=True)
        inp = directory/(job['name']+'.inp'); (attempt/inp.name).write_bytes(inp.read_bytes())
        job['attempt_history'] = [dict(job_id='900', archive=str(attempt), input_sha256=job['input_sha256'])]
        job['job_id'] = '901'; job['submission_status'] = 'submitted'
    atomic_json(attempt/'execution.json', dict(status='failed', slurm_job_id='900',
                                             input_sha256=job['input_sha256']))
    atomic_json(manifest, state)
    monkeypatch.setattr(submission_recovery, '_purged', lambda ident: 'Invalid job id specified')
    assert submission_recovery.terminal_evidence(root, claim)['job_id'] == '900'
    claims = {'token': claim}
    assert pool._reconcile(root, jobs, claims, archive=True) == 23 and claims == {}
    assert (root/'.kinbot_policy_reconciled/token.json').exists()
    atomic_json(attempt/'execution.json', dict(status='failed', slurm_job_id='999',
                                             input_sha256=job['input_sha256']))
    assert submission_recovery.terminal_evidence(root, claim) is None


def test_execution_records_allocation_even_before_native_staging(shared, monkeypatch):
    _, _, manifest, _ = shared
    job = json.loads(manifest.read_text())['jobs'][0]
    monkeypatch.setenv('SLURM_JOB_ID', '900')
    monkeypatch.setattr(recon, '_execute', lambda j: (_ for _ in ()).throw(ValueError('bad seed')))
    with pytest.raises(ValueError, match='bad seed'):
        recon.execute(manifest, job['name'])
    receipt = json.loads((Path(job['directory'])/'execution.json').read_text())
    assert receipt['slurm_job_id'] == '900' and receipt['status'] == 'failed'
