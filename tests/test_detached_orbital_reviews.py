"""Simulated detached review ownership; no native scientific approval."""
from types import SimpleNamespace

import pytest

from kinbot import requested_orbital_recon as review, run_policy, run_walltime
from kinbot.execution import atomic_json
from test_state_support_panel import (
    support_tree, support_plan, combined, spectrum, continuation_tree,
    resolution_tree, resolution_plan, tree, followup_plan, plan, original_tree, snapshot,
)


def test_detached_full_panel_run_resume_and_source_preservation(support_tree, support_plan,
                                                               tmp_path, monkeypatch):
    source, _ = support_tree
    pool = tmp_path/'review_pool'; pool.mkdir()
    run_policy.update(pool, max_jobs=3, exclude=['node9'], fair_submission=False)
    output = pool/'panel'
    before = snapshot(source)
    options = dict(shared_pool=True, resource_root=pool)
    assert review.batch(source, support_plan/'screens.json', output,
                        preflight=True, **options)['tasks'] == 52
    assert snapshot(source) == before and not output.exists()
    assert review.batch(source, support_plan/'screens.json', output, **options)['tasks'] == 52
    review.batch(source, support_plan/'screens.json', output, resume=True, **options)
    for system in review._saved(output)['systems']:
        state = review._read(system['manifest'])
        assert state['shared_policy_file'] == str(pool/run_policy.NAME)
    launched, running, done = [], {}, set()
    # The unrelated allocation must survive and must not consume this pool.
    unrelated = {'555': 1}
    def queue(root):
        assert root == pool
        return running.copy()
    monkeypatch.setattr(run_policy, '_queue', queue)
    monkeypatch.setattr(run_policy, 'slurm_job_record', lambda j:
        {'JobId': j, 'JobState': 'COMPLETED'} if j in done else {})
    monkeypatch.setattr(run_walltime, 'partition_walltime', lambda p: '31-00:00:00')
    def scheduler(command, **kwargs):
        if command[0] == 'bash':
            return SimpleNamespace(returncode=0, stdout='/mock/molpro\n', stderr='')
        if command[0] == 'squeue':
            return SimpleNamespace(returncode=0, stdout='\n'.join({**running, **unrelated}), stderr='')
        assert command[0] == 'sbatch' and len(running) < 3
        assert str(kwargs['cwd']).startswith(str(output)+'/')
        assert '--exclude=node9' in command
        job = str(1000+len(launched)); launched.append(str(kwargs['cwd'])); running[job] = 1
        return SimpleNamespace(returncode=0, stdout=job+'\n', stderr='')
    monkeypatch.setattr(review.subprocess, 'run', scheduler)
    def audit(path, **kwargs):
        state = review._read(path)
        for job in state['jobs']:
            job['status'] = 'complete' if job.get('job_id') in done else 'incomplete'
        state['status'] = ('complete_review_required' if all(j['status'] == 'complete'
                           for j in state['jobs']) else 'incomplete')
        atomic_json(path, state)
        return state
    monkeypatch.setattr(review.recon, '_audit', audit)
    monkeypatch.setattr(review.recon, '_recover_slurm_startup_failures', lambda *a: None)
    monkeypatch.setattr(review.recon, '_retry_allocation_failures', lambda p, m, *a: m)
    def finish(_):
        done.update(running); running.clear()
    monkeypatch.setattr(review.recon.time, 'sleep', finish)
    assert review.run(output)['status'] == 'complete_review_required'
    assert len(launched) == len(set(launched)) == 52
    review.run(output)
    assert len(launched) == 52 and unrelated == {'555': 1}
    assert snapshot(source) == before
    assert review._read(pool/'.kinbot_policy_admissions.json') == {}


def test_detached_pool_cannot_alias_source_or_move_on_resume(support_tree, support_plan, tmp_path):
    source, _ = support_tree
    with pytest.raises(ValueError, match='disjoint'):
        review.batch(source, support_plan/'screens.json', source/'panel',
                     shared_pool=True, resource_root=source)
    pool = tmp_path/'pool'; pool.mkdir()
    run_policy.update(pool, max_jobs=3)
    output = pool/'panel'
    review.batch(source, support_plan/'screens.json', output, shared_pool=True, resource_root=pool)
    with pytest.raises(ValueError, match='ownership differs'):
        review.batch(source, support_plan/'screens.json', output, shared_pool=True, resume=True)


def test_complete_request_inventory_is_read_only(support_tree):
    from kinbot.orbital_review_inventory import inventory
    source, _ = support_tree
    before = snapshot(source)
    report = inventory(source)
    paths = list(source.glob('*/[0-9]*/vrctst/orbital_review/*/request.json'))
    assert len(report['requests']) == len(paths) > 7
    assert sum(report['counts'].values()) == len(paths)
    assert all(r['sha256'] and 'installed_channels' in r for r in report['requests'])
    assert any(r['endpoints'] for r in report['requests'])
    assert snapshot(source) == before
