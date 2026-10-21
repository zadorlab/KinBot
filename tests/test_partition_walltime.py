"""Partition defaults are resolved at the scheduler boundary, never guessed."""
from types import SimpleNamespace

import pytest

from kinbot.anl import site
from kinbot import qc, run_walltime


def entries():
    return [dict(name='short', default=True, cores=96, memory_mb=126000, seconds=1800),
            dict(name='month', default=False, cores=96, memory_mb=126000, seconds=31*86400)]


@pytest.mark.parametrize('requested,expected', [(None, '31-00:00:00'),
                                              ('partition_max', '31-00:00:00'),
                                              ('02:00:00', '02:00:00')])
def test_fixed_resources_use_selected_partition_max_unless_explicit(monkeypatch, requested, expected):
    monkeypatch.setattr(site, '_partitions', entries)
    resources = dict(cores=8, memory_mb=64000, partition='month')
    if requested is not None:
        resources['walltime'] = requested
    spec = {'tasks': [dict(id='energy', resources=resources)]}
    site.assign_partitions(spec)
    assert resources['partition'] == 'month'
    assert resources['walltime'] == expected
    assert (resources['cores'], resources['memory_mb']) == (8, 64000)


def test_unlimited_partition_is_valid_slurm_zero(monkeypatch):
    unlimited = entries()[1]
    unlimited['seconds'] = float('inf')
    monkeypatch.setattr(site, '_partitions', lambda: [unlimited])
    task = dict(id='energy', resources=dict(cores=8, memory_mb=64000, partition='month'))
    site.assign_partitions({'tasks': [task]})
    assert task['resources']['walltime'] == '00:00:00'


def test_discovery_failure_cannot_substitute_an_arbitrary_time(monkeypatch):
    monkeypatch.setattr(site, '_partitions', lambda: [])
    with pytest.raises(RuntimeError, match='maximum walltime require sinfo'):
        site.assign_partitions({'tasks': [dict(id='energy', resources=dict(
            cores=8, memory_mb=64000, partition='month'))]})


@pytest.fixture
def scheduler(monkeypatch):
    monkeypatch.setattr(qc, '_SLURM_WALLTIME_CACHE', {})
    calls = []
    def run(args, **kwargs):
        calls.append(args)
        return SimpleNamespace(returncode=0, stdout='31-00:00:00\n')
    monkeypatch.setattr(qc.subprocess, 'run', run)
    return calls


def test_ordinary_default_max_and_explicit_shorter_override(scheduler):
    assert qc._ordinary_slurm_walltime({'queue_name': 'month'}) == '#SBATCH --time=31-00:00:00'
    assert scheduler == [['sinfo', '-h', '-p', 'month', '-o', '%l']]
    assert qc._ordinary_slurm_walltime({'queue_name': 'month', 'queue_walltime': '00:30:00'}) == '#SBATCH --time=00:30:00'
    assert len(scheduler) == 1


def test_ordinary_discovery_unknown_is_not_cached_as_success(monkeypatch, scheduler):
    from kinbot.anl.dispatch import SchedulerUnavailableError
    monkeypatch.setattr(qc.subprocess, 'run', lambda *a, **k: SimpleNamespace(
        returncode=0, stdout='31-00:00:00\nNOT_SET\n'))
    with pytest.raises(SchedulerUnavailableError, match='no job was submitted'):
        qc._ordinary_slurm_walltime({'queue_name': 'month'})
    assert not qc._SLURM_WALLTIME_CACHE


def test_ordinary_refreshes_expired_partition_limit(scheduler, monkeypatch):
    monkeypatch.setattr(qc.time, 'monotonic', lambda: 100)
    qc._SLURM_WALLTIME_CACHE['month'] = (1, '#SBATCH --time=7-00:00:00')
    assert qc._ordinary_slurm_walltime({'queue_name': 'month'}) == '#SBATCH --time=31-00:00:00'
    assert len(scheduler) == 1


def test_rotd_default_is_independent_of_ordinary_override(scheduler, tmp_path, monkeypatch):
    from kinbot import pes
    monkeypatch.setattr(pes, 'render_site_setup', lambda *a: '# mock site\n')
    path = pes._write_rotdpy_slurm_template(dict(queuing='slurm', queue_name='month',
        queue_walltime='00:30:00', rotdpy_walltime='partition_max'), tmp_path)
    from pathlib import Path
    assert '#SBATCH --time=31-00:00:00' in Path(path).read_text()


def test_recon_resolves_max_at_submit_and_preserves_environment_override(monkeypatch):
    from kinbot.molpro_reconnaissance import _slurm_submission_options
    monkeypatch.delenv('SBATCH_PARTITION', raising=False)
    monkeypatch.delenv('SBATCH_TIMELIMIT', raising=False)
    monkeypatch.setattr(run_walltime, 'partition_limits', lambda **kw: {'month': 31*86400})
    batch = '#!/bin/bash\n#SBATCH -p month\n# KINBOT_WALLTIME=partition_max\nset -e\n'
    assert '--time=31-00:00:00' in _slurm_submission_options(batch, [])
    monkeypatch.setenv('SBATCH_TIMELIMIT', '02:00:00')
    assert not any(x.startswith('--time') for x in _slurm_submission_options(batch, []))


def test_recon_explicit_time_needs_no_discovery(monkeypatch):
    from kinbot.molpro_reconnaissance import _slurm_submission_options
    def unexpected(**kw):
        pytest.fail('Explicit time must not query partitions')
    monkeypatch.setattr(run_walltime, 'partition_limits', unexpected)
    assert not any(x.startswith('--time') for x in _slurm_submission_options(
        '#!/bin/bash\n#SBATCH -p month\n#SBATCH --time=01:00:00\nset -e\n', []))
