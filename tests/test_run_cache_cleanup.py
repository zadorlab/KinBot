import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import run_cache_cleanup as cleanup
from kinbot.execution import atomic_json


@pytest.fixture
def cache(tmp_path, monkeypatch):
    home = tmp_path/'home'
    home.mkdir()
    archive = home/'runs/r6.stopped'
    archive.mkdir(parents=True)
    new = home/'runs/r7'
    new.mkdir()
    wfu = home/'wfu'
    wfu.mkdir()
    (wfu/'old.wfu').write_text('disposable native file')
    (archive/'approved.wfu').write_text('retained r6 result')
    (new/'seed.wfu').write_text('reviewed seed')
    (home/'.cache').mkdir()
    (home/'.cache/package').write_text('not KinBot scratch')
    monkeypatch.setattr(cleanup.Path, 'home', classmethod(lambda cls: home))
    monkeypatch.setattr(cleanup, '_mounts', lambda: [])
    checks = []
    monkeypatch.setattr(cleanup, 'assert_idle', lambda paths: checks.append(list(paths)))
    options = {'empty_home_wfu': True, 'remove_verified_home_scratch': True}
    return SimpleNamespace(home=home, archive=archive, new=new, wfu=wfu, options=options, checks=checks)


def test_only_authorized_wfu_and_provenance_matched_scratch_removed(cache):
    c = cache
    root = c.home/'.kinbot_scratch'
    for name in ('kinbot-molpro-old', 'unknown', 'kinbot-molpro-retained', 'kinbot-molpro-otherhost'):
        p = root/name
        p.mkdir(parents=True)
        (p/'data').write_text(name)
    host = cleanup.socket.gethostname()
    for name, props in [('old', {}), ('retained', {'retain_scratch': True}),
                        ('otherhost', {'scratch_host': 'other-compute-node'})]:
        atomic_json(c.archive/(name+'.runtime.json'), dict(scratch_directory=str(root/('kinbot-molpro-'+name)),
                    scratch_host=host, **{k: v for k, v in props.items() if k != 'scratch_host'})
                    | ({'scratch_host': props['scratch_host']} if 'scratch_host' in props else {}))
    report = cleanup.cleanup(c.archive, c.new, c.options)
    assert report['status'] == 'complete'
    assert list(c.wfu.iterdir()) == []
    assert not (root/'kinbot-molpro-old').exists()
    assert (root/'unknown/data').exists()
    assert (root/'kinbot-molpro-retained/data').exists()
    assert (root/'kinbot-molpro-otherhost/data').exists()
    assert len(report['skipped']) == 3
    assert (c.archive/'approved.wfu').read_text() == 'retained r6 result'
    assert (c.new/'seed.wfu').read_text() == 'reviewed seed'
    assert (c.home/'.cache/package').exists()
    assert len(c.checks) == 3
    # Repeat is a no-op even if other work later recreated cache contents.
    (c.wfu/'new.wfu').write_text('new work')
    cleanup.cleanup(c.archive, c.new, c.options)
    assert (c.wfu/'new.wfu').read_text() == 'new work'
    assert len(c.checks) == 3


def test_queue_or_local_guard_failure_changes_no_cache(cache, monkeypatch):
    def busy(paths):
        raise RuntimeError('another allocation')
    monkeypatch.setattr(cleanup, 'assert_idle', busy)
    with pytest.raises(RuntimeError, match='another allocation'):
        cleanup.cleanup(cache.archive, cache.new, cache.options)
    assert (cache.wfu/'old.wfu').exists()


def test_new_job_after_move_leaves_retired_data_and_retries_safely(cache, monkeypatch):
    checks = []
    def guard(paths):
        checks.append(1)
        if len(checks) > 1:
            raise RuntimeError('new job appeared')
    monkeypatch.setattr(cleanup, 'assert_idle', guard)
    with pytest.raises(RuntimeError, match='new job appeared'):
        cleanup.cleanup(cache.archive, cache.new, cache.options)
    journal = json.loads((cache.new/'cache_cleanup.json').read_text())
    assert (Path(journal['entries'][0]['retired'])/'old.wfu').exists()
    (cache.wfu/'new.wfu').write_text('keep')
    monkeypatch.setattr(cleanup, 'assert_idle', lambda paths: None)
    cleanup.cleanup(cache.archive, cache.new, cache.options)
    assert (cache.wfu/'new.wfu').read_text() == 'keep'


@pytest.mark.parametrize('boundary', ['rename', 'rmtree'])
def test_interrupted_cleanup_preserves_new_contents_and_resumes_retired_tree(cache, monkeypatch, boundary):
    c = cache
    rename, remove = Path.rename, cleanup.shutil.rmtree
    def interrupt_rename(path, dest):
        result = rename(path, dest)
        if path == c.wfu:
            raise RuntimeError('connection lost after move')
        return result
    def interrupt_remove(path, *args, **kwargs):
        remove(path, *args, **kwargs)
        raise RuntimeError('connection lost after delete')
    interrupt_remove.avoids_symlink_attacks = True
    if boundary == 'rename':
        monkeypatch.setattr(Path, 'rename', interrupt_rename)
    else:
        monkeypatch.setattr(cleanup.shutil, 'rmtree', interrupt_remove)
    with pytest.raises(RuntimeError, match='connection lost'):
        cleanup.cleanup(c.archive, c.new, c.options)
    c.wfu.mkdir(exist_ok=True)
    (c.wfu/'new.wfu').write_text('new work')
    monkeypatch.setattr(Path, 'rename', rename)
    monkeypatch.setattr(cleanup.shutil, 'rmtree', remove)
    report = cleanup.cleanup(c.archive, c.new, c.options)
    assert report['status'] == 'complete'
    assert (c.wfu/'new.wfu').read_text() == 'new work'
    assert not list(c.home.glob('.kinbot-retired-cache-*'))


def test_root_symlink_is_refused_but_internal_symlink_never_followed(cache):
    c = cache
    c.wfu.rename(c.home/'saved-wfu')
    c.wfu.symlink_to(c.home/'saved-wfu', target_is_directory=True)
    with pytest.raises(RuntimeError, match='non-symlink'):
        cleanup.cleanup(c.archive, c.new, c.options)
    c.wfu.unlink()
    (c.home/'saved-wfu').rename(c.wfu)
    (c.wfu/'link').symlink_to(c.archive, target_is_directory=True)
    cleanup.cleanup(c.archive, c.new, c.options)
    assert (c.archive/'approved.wfu').exists()


def test_mount_and_protected_overlap_refused_before_deletion(cache, monkeypatch):
    c = cache
    monkeypatch.setattr(cleanup, '_mounts', lambda: [c.wfu/'nested'])
    with pytest.raises(RuntimeError, match='mount point'):
        cleanup.cleanup(c.archive, c.new, c.options)
    assert (c.wfu/'old.wfu').exists()
    monkeypatch.setattr(cleanup, '_mounts', lambda: [])
    with pytest.raises(RuntimeError, match='protected'):
        cleanup.cleanup(c.archive, c.new, c.options, protected=[c.wfu/'old.wfu'])


def test_account_name_and_all_allocations_guard(monkeypatch):
    monkeypatch.setattr(cleanup, 'account_identity', lambda: {'name': 'researcher', 'uid': 1001})
    def queue(args, **kw):
        assert args == ['squeue', '-h', '-u', 'researcher', '-o', '%i|%T']
        assert kw['check'] and kw['timeout'] == 20
        return SimpleNamespace(stdout='17|COMPLETING\n18|PENDING\n')
    monkeypatch.setattr(cleanup.subprocess, 'run', queue)
    with pytest.raises(RuntimeError, match='ALL account allocations'):
        cleanup.assert_idle([])


def test_idle_guard_detects_local_qc_and_open_cache(tmp_path, monkeypatch):
    proc = tmp_path/'proc'
    process = proc/'12345'
    process.mkdir(parents=True)
    (process/'fd').mkdir()
    (process/'cwd').symlink_to(tmp_path, target_is_directory=True)
    (process/'cmdline').write_bytes(b'/opt/molpro/bin/molpro\0job.inp\0')
    original = cleanup.Path
    monkeypatch.setattr(cleanup, 'Path', lambda p: proc if str(p) == '/proc' else original(p))
    monkeypatch.setattr(cleanup.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout=''))
    with pytest.raises(RuntimeError, match='local QC/controller'):
        cleanup.assert_idle([])
    (process/'cmdline').write_bytes(b'cat\0cache.wfu\0')
    (process/'fd/8').symlink_to(tmp_path/'wfu/file')
    with pytest.raises(RuntimeError, match='uses'):
        cleanup.assert_idle([tmp_path/'wfu'])
    (process/'fd/8').unlink()
    cleanup.assert_idle([tmp_path/'wfu'])


def test_publication_failure_preserves_home_scratch(cache):
    c = cache
    scratch = c.home/'.kinbot_scratch/kinbot-molpro-failed'
    scratch.mkdir(parents=True)
    atomic_json(c.archive/'calculation.runtime.json', dict(scratch_directory=str(scratch),
        scratch_host=cleanup.socket.gethostname()))
    atomic_json(c.archive/'execution.json', {'error': 'QC result publication failed; local files retained'})
    cleanup.cleanup(c.archive, c.new, c.options)
    assert scratch.exists()
