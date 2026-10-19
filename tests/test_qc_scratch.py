"""Native processes must never spill integrals into HOME/result repositories."""
import errno
import json
import os
import shlex
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import pytest

from kinbot import qc_scratch as scratch
from kinbot.anl.runtime import _check_scratch_io, cleanup_qc_runtime, qc_runtime_environment


def test_actual_account_not_environment_username(monkeypatch):
    monkeypatch.setenv('USER', 'wrong-user')
    monkeypatch.setattr(scratch.pwd, 'getpwuid', lambda uid: SimpleNamespace(pw_name='alice', pw_dir='/home/alice'))
    assert scratch.account_scratch() == Path('/scratch/alice')


@pytest.fixture
def missing_passwd(monkeypatch):
    monkeypatch.setattr(scratch.os, 'getuid', lambda: 155382)
    def unavailable(uid):
        raise KeyError(f'getpwuid(): uid not found: {uid}')
    monkeypatch.setattr(scratch.pwd, 'getpwuid', unavailable)
    return unavailable


@pytest.mark.parametrize('identity,source', [
    ({'USER': 'alice'}, 'USER'),
    ({'LOGNAME': 'alice'}, 'LOGNAME'),
    ({'USER': 'different', 'SLURM_JOB_UID': '155382', 'SLURM_JOB_USER': 'alice'}, 'slurm_uid'),
    ({'KINBOT_QC_ACCOUNT': json.dumps({'uid': 155382, 'name': 'alice', 'home': '/home/alice'})}, 'inherited_uid'),
])
def test_missing_passwd_uses_exported_account_name(missing_passwd, identity, source):
    assert scratch.account_scratch(identity) == Path('/scratch/alice')
    assert scratch.account_identity(identity)['source'] == source


@pytest.mark.parametrize('identity', [
    {}, {'USER': '../alice'}, {'USER': '/home/alice'}, {'USER': ''},
    {'SLURM_JOB_UID': '999', 'SLURM_JOB_USER': 'alice'},
    {'KINBOT_QC_ACCOUNT': '{broken'},
    {'KINBOT_QC_ACCOUNT': json.dumps({'uid': 999, 'name': 'alice'})},
])
def test_missing_passwd_without_usable_name_fails_clearly(missing_passwd, identity):
    with pytest.raises(RuntimeError, match='Cannot resolve QC account'):
        scratch.account_scratch(identity)


def test_controller_carries_name_when_compute_node_has_no_nss(monkeypatch):
    monkeypatch.setattr(scratch.pwd, 'getpwuid', lambda uid:
                        SimpleNamespace(pw_name='alice', pw_dir='/home/alice'))
    monkeypatch.setenv('KINBOT_QC_ACCOUNT', 'old')
    scratch.configure({})
    record = os.environ['KINBOT_QC_ACCOUNT']
    def unavailable(uid):
        raise KeyError(uid)
    monkeypatch.setattr(scratch.pwd, 'getpwuid', unavailable)
    assert scratch.account_scratch({'KINBOT_QC_ACCOUNT': record, 'USER': 'wrong'}) == Path('/scratch/alice')


@pytest.mark.parametrize('backend', ['molpro', 'gaussian', 'cfour', 'mrcc', 'qchem', 'orca'])
def test_all_backends_allocate_with_user_and_missing_passwd(tmp_path, monkeypatch, missing_passwd, backend):
    preferred = tmp_path/'scratch/alice'
    preferred.parent.mkdir()
    original = scratch.account_scratch
    monkeypatch.setattr(scratch, 'account_scratch', lambda env=None:
                        tmp_path/'scratch'/original(env).name)
    monkeypatch.setattr(Path, 'home', lambda: pytest.fail('must not require a passwd-based HOME lookup'))
    home = tmp_path/'home'
    path, info = scratch.allocate({'USER': 'alice', 'HOME': str(home),
        'KINBOT_QC_SCRATCH': str(home/'old-cache')}, backend, 1, probe=_check_scratch_io)
    try:
        assert path.parent == preferred
        assert info['scratch_source'] == 'account_scratch'
        assert info['scratch_account']['name'] == 'alice'
        assert not home.exists()
    finally:
        cleanup_qc_runtime(info)


def test_missing_passwd_does_not_allow_home_fallback(tmp_path, monkeypatch, missing_passwd):
    home = tmp_path/'home'
    monkeypatch.setattr(scratch, 'candidates', lambda *a: iter([('override', home/'cache')]))
    with pytest.raises(OSError, match='Home directories cannot'):
        scratch.allocate({'USER': 'alice', 'HOME': str(home)}, 'molpro', 1, probe=_check_scratch_io)
    assert not home.exists()


@pytest.mark.parametrize('backend', ['molpro', 'gaussian', 'cfour', 'mrcc', 'qchem', 'orca'])
def test_account_scratch_precedes_even_explicit_legacy_override(tmp_path, monkeypatch, backend):
    preferred = tmp_path / 'scratch' / 'alice'
    preferred.parent.mkdir()
    monkeypatch.setattr(scratch, 'account_scratch', lambda env=None: preferred)
    home = tmp_path / 'home'
    env = {'HOME': str(home), 'KINBOT_MOLPRO_SCRATCH': str(home/'cache'),
           'KINBOT_QC_SCRATCH': str(tmp_path/'other'), 'SLURM_NNODES': '1'}
    path, info = scratch.allocate(env, backend, 1, probe=_check_scratch_io)
    try:
        assert path.parent == preferred
        assert info['scratch_source'] == 'account_scratch'
        assert info['scratch_io_checked'] and info['scratch_host']
        assert not home.exists()
        assert not (tmp_path/'other').exists()
    finally:
        cleanup_qc_runtime(info)


def test_shared_filesystem_and_home_override_are_rejected(tmp_path, monkeypatch):
    home, nfs, local = (tmp_path/n for n in ('home', 'network', 'local'))
    monkeypatch.setattr(scratch, 'candidates', lambda *a: iter([
        ('override', home/'cache'), ('site', nfs), ('local', local)]))
    real = scratch.filesystem
    def mount(path):
        if path == nfs:
            raise OSError('network filesystem nfs4')
        return real(path)
    monkeypatch.setattr(scratch, 'filesystem', mount)
    path, info = scratch.allocate({'HOME': str(home)}, 'molpro', 1, probe=_check_scratch_io)
    try:
        assert path.parent == local
        assert len(info['scratch_rejected_candidates']) == 2
        assert not home.exists() and not nfs.exists()
    finally:
        cleanup_qc_runtime(info)


def test_linux_mount_detection_uses_deepest_mount_including_symlinks(tmp_path, monkeypatch):
    mount = tmp_path / 'shared space'
    mount.mkdir()
    link = tmp_path / 'alias'
    link.symlink_to(mount)
    native_read = Path.read_text
    escaped_mount = str(mount).replace(' ', r'\040')
    records = ('1 0 8:1 / / rw - ext4 /dev/sda rw\n'
               f'2 1 0:2 / {escaped_mount} rw - nfs4 server:/export rw\n')
    monkeypatch.setattr(scratch.sys, 'platform', 'linux')
    monkeypatch.setattr(Path, 'read_text', lambda p, *a, **kw:
                        records if str(p) == '/proc/self/mountinfo' else native_read(p, *a, **kw))
    with pytest.raises(OSError, match='nfs4'):
        scratch.filesystem(link/'child')
    assert scratch.filesystem(tmp_path)['filesystem'] == 'ext4'


def test_io_error_keeps_retryable_errno_and_no_home_fallback(tmp_path, monkeypatch):
    monkeypatch.setattr(scratch, 'candidates', lambda *a: iter([('site', tmp_path/'local')]))
    def broken(path):
        raise OSError(errno.EIO, 'Input/output error')
    with pytest.raises(OSError) as caught:
        scratch.allocate({}, 'molpro', 1, probe=broken)
    assert caught.value.errno == errno.EIO
    assert not list((tmp_path/'local').iterdir())


def test_multinode_is_rejected_before_local_allocation():
    with pytest.raises(ValueError, match='single-node'):
        scratch.allocate({'SLURM_NNODES': '2'}, 'molpro', 1, probe=_check_scratch_io)


@pytest.mark.parametrize('backend,deck,checkpoint', [
    ('molpro', 'file,2,seed.wfu\nput,xml,job.xml\n', 'seed.wfu'),
    ('gaussian', '%chk=seed\n#p hf/sto-3g\n', 'seed.chk')])
def test_workspace_stages_restart_and_publishes_results_only(tmp_path, backend, deck, checkpoint):
    result, local = tmp_path/'results', tmp_path/'scratch'
    result.mkdir(); local.mkdir()
    inp = 'job.inp' if backend == 'molpro' else 'job.com'
    (result/inp).write_text(deck)
    (result/checkpoint).write_text('restart')
    info = {'scratch_directory': str(local)}
    with scratch.workspace(result, info, backend=backend, inputs=[inp], label='job') as work:
        assert work == local/'work' and (work/checkpoint).read_text() == 'restart'
        (work/checkpoint).write_text('updated checkpoint')
        (work/'job.out').write_text('native result')
        (work/'huge.integrals').write_bytes(b'intermediate')
        (work/'job.xml').write_text('orbitals')
        assert not (result/'job.out').exists()
    assert (result/'job.out').read_text() == 'native result'
    assert (result/checkpoint).read_text() == 'updated checkpoint'
    assert not (result/'huge.integrals').exists()
    assert json.loads((result/'job.runtime.json').read_text())['native_work_directory'] == str(work)
    cleanup_qc_runtime(info)
    assert not local.exists()


def test_copyback_failure_retains_local_results_even_after_cleanup(tmp_path, monkeypatch):
    result, local = tmp_path/'results', tmp_path/'scratch'
    result.mkdir(); local.mkdir()
    (result/'job.inp').write_text('rhf\n')
    info = {'scratch_directory': str(local)}
    copy = scratch.shutil.copy2
    def fail_output(src, dest):
        if Path(src).name == 'job.out':
            raise OSError(errno.EIO, 'home filesystem unavailable')
        return copy(src, dest)
    monkeypatch.setattr(scratch.shutil, 'copy2', fail_output)
    with pytest.raises(RuntimeError, match='local files retained'):
        with scratch.workspace(result, info, backend='molpro', inputs=['job.inp']) as work:
            (work/'job.out').write_text('expensive result')
    cleanup_qc_runtime(info)
    assert (local/'work/job.out').read_text() == 'expensive result'


@pytest.mark.parametrize('backend,text', [
    ('molpro', 'file,2,/home/user/wfu/old.wfu\n'),
    ('molpro', 'file,2,local.wfu;file,3,/home/user/seed.wfu\n'),
    ('molpro', 'put,xml,../orbital.xml\n'),
    ('gaussian', '%chk=/home/user/a.chk\n'),
    ('gaussian', '%rwf=/home/user/a.rwf\n')])
def test_absolute_native_output_and_seed_paths_cannot_bypass_locality(backend, text):
    with pytest.raises(ValueError, match='relative basename'):
        scratch.references(text, backend)


def test_molpro_repository_flags_are_explicit_and_unique():
    args = scratch.molpro_command(['molpro', '-W', '/home/wfu', '-d/old', '-I', '/shared',
                                  '-n', '8', 'point.inp'], {'scratch_directory': '/scratch/a'}, '/scratch/a/work')
    assert args == ['molpro', '-d', '/scratch/a', '-W', '/scratch/a/work',
                    '-I', '/scratch/a', '-n', '8', 'point.inp']


def test_cfour_library_failure_allocates_no_scratch(tmp_path, monkeypatch):
    monkeypatch.setattr('kinbot.anl.runtime.shutil.which', lambda *a, **k: None)
    def unexpected(*a, **kw):
        pytest.fail('scratch allocated before missing executable was checked')
    monkeypatch.setattr(scratch, 'allocate', unexpected)
    with pytest.raises(RuntimeError, match='executable unavailable'):
        qc_runtime_environment('xcfour', 'cfour', work_directory=tmp_path)


def test_gaussian_executes_locally_and_preserves_checkpoint_between_calls(tmp_path, monkeypatch):
    from kinbot.ase_modules.calculators.gaussian import Gaussian
    result, local = tmp_path/'results', tmp_path/'local'
    result.mkdir()
    monkeypatch.setattr(scratch, 'candidates', lambda *a: iter([('test', local)]))
    (result/'job.com').write_text('%chk=job\n#p HF/STO-3G\n')
    (result/'job.chk').write_text('seed')
    native = tmp_path/'native.py'
    native.write_text('from pathlib import Path\nimport os,sys\n'
        'assert Path.cwd().parent == Path(os.environ["GAUSS_SCRDIR"])\n'
        'assert "%chk=job" in sys.stdin.read()\n'
        'p=Path("job.chk"); p.write_text(p.read_text()+"+updated")\n'
        'Path("huge.rwf").write_bytes(b"scratch")\n'
        'print("Normal termination of Gaussian 16")\n')
    calculator = SimpleNamespace(prefix='job', directory=result,
        command=f'{shlex.quote(sys.executable)} {shlex.quote(str(native))} < PREFIX.com > PREFIX.log')
    Gaussian.execute(calculator)
    Gaussian.execute(calculator)
    assert (result/'job.chk').read_text() == 'seed+updated+updated'
    assert 'Normal termination' in (result/'job.log').read_text()
    assert not (result/'huge.rwf').exists() and not list(local.iterdir())


@pytest.mark.parametrize('backend', ['cfour', 'mrcc'])
def test_external_native_backend_runs_locally_and_retains_alternative_output_names(tmp_path, monkeypatch, backend):
    from kinbot.anl.dispatch import _run_external
    result, local = tmp_path/'results', tmp_path/'local'
    result.mkdir()
    monkeypatch.setattr(scratch, 'candidates', lambda *a: iter([('test', local)]))
    (result/'MINP').write_text('native input')
    (result/'geometry.xyz').write_text('1\nfixture\nH 0 0 0\n')
    native = tmp_path/'native'
    native.write_text(f'#!{sys.executable}\nfrom pathlib import Path\nimport os\n'
        'assert Path.cwd().parent == Path(os.environ["TMPDIR"])\n'
        'assert Path("MINP").read_text() == "native input"\n'
        'Path("fort.20").write_text("completed result")\n'
        'Path("huge.integrals").write_bytes(b"scratch")\n')
    native.chmod(0o755)
    task = {'id': 'test', 'backend': backend, 'kind': 'external',
        'resources': {'memory_mb': 8192, 'cores': 2}, 'command': [str(native)],
        'input_name': 'MINP', 'required_outputs': ['energy.dat'],
        'output_aliases': {'energy.dat': ['fort.20']}}
    details = _run_external(result, {'task': task})
    assert details['returncode'] == 0
    assert (result/'energy.dat').read_text() == 'completed result'
    assert not (result/'huge.integrals').exists() and not list(local.iterdir())
