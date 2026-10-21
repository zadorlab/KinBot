"""Walltime changes are scoped scheduler requests, never QC restarts."""
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from kinbot import run_walltime as wt


@pytest.mark.parametrize('text,expected', [('1-02:03:04', 93784), ('23:30:00', 84600),
    ('3:20', 200), ('30', 1800), ('2-03', 183600), ('UNLIMITED', float('inf'))])
def test_slurm_duration_formats(text, expected):
    assert wt.seconds(text) == expected


@pytest.mark.parametrize('text', ['0:00:00', '24:60:00', '1', '-1:00:00', 'UNLIMITED', None])
def test_invalid_policy_duration(text):
    with pytest.raises(ValueError):
        wt.validate_minimum(text)


def test_submission_caps_partition_and_preserves_environment_and_native_script(tmp_path, monkeypatch):
    script = tmp_path/'job.slurm'
    text = '#!/bin/bash\n#SBATCH --time=12:00:00\n#SBATCH --partition=week\nrun_qc\n'
    script.write_text(text)
    monkeypatch.setattr(wt, 'partition_limits', lambda **kw: {'week': 604800, 'day': 86400})
    assert wt.submission_override(script, '7-00:00:00', env={})['effective'] == '7-00:00:00'
    override = wt.submission_override(script, '7-00:00:00', env={'SBATCH_PARTITION': 'day'})
    assert override['effective'] == '1-00:00:00' and override['capped_to_partition']
    assert wt.submission_override(script, '7-00:00:00', env={'SBATCH_TIMELIMIT': '10-00:00:00'}) == {}
    assert wt.submission_override(script, '7-00:00:00', env={'SBATCH_TIMELIMIT': '0'}) == {}
    assert script.read_text() == text


def test_missing_partition_fails_and_implicit_time_uses_site_default(tmp_path, monkeypatch):
    script = tmp_path/'job.slurm'; script.write_text('#SBATCH --time=1:00:00\n')
    monkeypatch.setattr(wt, 'partition_limits', lambda **kw: {'week': 604800})
    with pytest.raises(RuntimeError, match='partition'):
        wt.submission_override(script, '7-00:00:00', env={})
    script.write_text('echo hello\n#SBATCH --time=1:00:00\n')
    monkeypatch.setattr(wt, 'partition_defaults', lambda: {'': 10*86400})
    assert wt.submission_override(script, '7-00:00:00', env={}) == {}
    monkeypatch.setattr(wt, 'partition_defaults', lambda: {'': 86400})
    monkeypatch.setattr(wt, 'partition_limits', lambda **kw: {'': 31*86400})
    assert wt.submission_override(script, '7-00:00:00', env={})['effective'] == '7-00:00:00'


@pytest.fixture
def scheduler(tmp_path, monkeypatch):
    root = tmp_path/'run'; root.mkdir()
    calls = []
    queue = [{'id': '101', 'state': 'RUNNING', 'limit': '1-00:00:00', 'partition': 'long',
              'directory': str(root/'case'), 'elapsed': '18:00:00'},
             {'id': '102', 'state': 'PENDING', 'limit': '12:00:00', 'partition': 'long',
              'directory': str(root/'case'), 'elapsed': '0:00'},
             {'id': '103', 'state': 'RUNNING', 'limit': '10-00:00:00', 'partition': 'long',
              'directory': str(root/'case'), 'elapsed': '18:00:00'},
             {'id': '999', 'state': 'RUNNING', 'limit': '12:00:00', 'partition': 'long',
              'directory': str(tmp_path/'other'), 'elapsed': '18:00:00'}]
    state = SimpleNamespace(root=root, queue=queue, calls=calls, deny=False,
                            timeout=False, drop=False, queries=0, before_update=None)
    def command(args):
        calls.append(args)
        if args[0] == 'sinfo':
            return SimpleNamespace(stdout='long*|31-00:00:00\nshort|1-00:00:00\n')
        if args[0] == 'squeue':
            state.queries += 1
            if state.queries == 2 and state.before_update:
                state.before_update(queue)
            return SimpleNamespace(stdout='\n'.join('|'.join(j[k] for k in
                ('id', 'state', 'elapsed', 'limit', 'partition', 'directory')) for j in queue))
        assert args[:2] == ['scontrol', 'update']
        ids = args[2].split('=', 1)[1].split(',')
        assert set(ids) <= {'101', '102'}
        assert args[3] == 'TimeLimit=7-00:00:00'
        if state.deny:
            # Slurm may apply a pending update while refusing a running one.
            queue[1]['limit'] = '7-00:00:00'
            raise subprocess.CalledProcessError(1, args, stderr='Access/permission denied')
        for j in queue:
            if j['id'] in ids:
                j['limit'] = args[3].split('=', 1)[1]
        if state.drop:
            queue[:] = [j for j in queue if j['id'] != '101']
        if state.timeout:
            raise subprocess.TimeoutExpired(args, 20)
        return SimpleNamespace(stdout='')
    monkeypatch.setattr(wt, 'command', command)
    return state


def test_dry_run_is_scoped_readonly_and_apply_verifies_without_cancelling(scheduler):
    before = deepcopy(scheduler.queue)
    report = wt.extend(scheduler.root, '7-00:00:00')
    assert scheduler.queue == before
    assert {j['id'] for j in report['selected']} == {'101', '102'}
    assert all(c[0] in ('sinfo', 'squeue') for c in scheduler.calls)
    report = wt.extend(scheduler.root, '7-00:00:00', apply=True)
    assert {j['id'] for j in report['verified']} == {'101', '102'}
    assert report['unresolved'] == []
    assert len([c for c in scheduler.calls if c[0] == 'scontrol']) == 1
    assert scheduler.queue[2:] == before[2:]


def test_partial_permission_error_and_vanished_jobs_are_not_success(scheduler):
    scheduler.deny = True
    report = wt.extend(scheduler.root, '7-00:00:00', apply=True)
    assert [j['id'] for j in report['verified']] == ['102']
    assert [j['id'] for j in report['unresolved']] == ['101']
    assert 'denied' in report['request_error']
    scheduler.deny = False; scheduler.drop = True
    report = wt.extend(scheduler.root, '7-00:00:00', apply=True)
    assert report['unresolved'][0]['after'] is None


def test_timeout_is_reconciled_not_repeated(scheduler):
    scheduler.timeout = True
    report = wt.extend(scheduler.root, '7-00:00:00', apply=True)
    assert len(report['verified']) == 2 and report['request_error']
    assert len([c for c in scheduler.calls if c[0] == 'scontrol']) == 1


def test_partition_cap_and_fresh_limit_increase_are_respected(scheduler):
    scheduler.queue[0]['partition'] = 'short'
    scheduler.before_update = lambda queue: queue[1].update(limit='14-00:00:00')
    report = wt.extend(scheduler.root, '7-00:00:00', apply=True)
    assert [j['id'] for j in report['skipped']] == ['101']
    assert not report['command']
    assert not any(c[0] == 'scontrol' for c in scheduler.calls)


def test_aggregate_array_and_missing_workdir_fail_closed(scheduler):
    scheduler.queue[0]['id'] = '101_[1-8]'
    with pytest.raises(ValueError, match='aggregate'):
        wt.extend(scheduler.root, '7-00:00:00', apply=True)
    scheduler.queue[0]['id'] = '101'
    scheduler.queue[0]['directory'] = '(null)'
    with pytest.raises(ValueError, match='work directory'):
        wt.extend(scheduler.root, '7-00:00:00', apply=True)
    assert not any(c[0] == 'scontrol' for c in scheduler.calls)


def test_not_established_limit_is_reported_not_modified(scheduler):
    scheduler.queue[1]['limit'] = 'NOT_SET'
    report = wt.extend(scheduler.root, '7-00:00:00', apply=True)
    assert [j['id'] for j in report['selected']] == ['101']
    assert report['skipped'][0]['id'] == '102'
    assert wt.duration(report['skipped'][0]['remaining_seconds']) == 'unknown'


def test_progress_compares_content_not_mirrored_mtime_and_reports_curvature(scheduler, monkeypatch):
    directory = scheduler.root/'case'; directory.mkdir()
    log = directory/'energy.live.log'
    log.write_text('Iteration 16 CC energy: -437.46024626\nNorm of residual: 0.00000383\nWall time [min]: 1057.278\n')
    curv = directory/'minimum_curvature_0.json'
    curv.write_text(json.dumps({'status': 'running', 'probes': [1, 2, 3],
                               'geometry_angstrom': [[0, 0, 0], [1, 0, 0]]}))
    snapshot = {p.name: p.read_bytes() for p in directory.iterdir()}
    monkeypatch.setattr(wt.time, 'sleep', lambda _: os.utime(log, None))
    counters = iter([({'101.batch': ['01:00:00', '1G', '100', '200']}, None),
                     ({'101.batch': ['01:00:20', '1G', '101', '200']}, None)])
    monkeypatch.setattr(wt, 'counters', lambda _: next(counters))
    report = wt.progress(scheduler.root)
    job = report['jobs'][0]
    assert job['native_changed'] is False
    assert job['cpu_changed'] is True and job['io_changed'] is True
    assert job['remaining_seconds'] == 6*3600
    assert job['evidence']['curvature']['total'] == 12
    assert job['evidence']['curvature']['completed'] == 3
    assert any('residual' in row for row in job['evidence']['native_rows'])
    assert {p.name: p.read_bytes() for p in directory.iterdir()} == snapshot
    assert all(c[0] == 'squeue' for c in scheduler.calls)


def test_cli_cannot_overwrite_run_or_mix_progress_with_apply(scheduler):
    with pytest.raises(SystemExit):
        wt.main([str(scheduler.root), '--progress', '--apply'])
    with pytest.raises(SystemExit):
        wt.main([str(scheduler.root), '--output', str(scheduler.root/'evidence.json')])
