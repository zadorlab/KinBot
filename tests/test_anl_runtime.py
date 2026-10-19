"""Native QC runtime checks are scoped to the launched program."""

from pathlib import Path
import errno
import shutil
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from kinbot.anl.runtime import cleanup_qc_runtime, qc_runtime_environment


def test_cfour_finds_only_missing_fortran_library_for_child():
    with TemporaryDirectory() as temporary:
        library = Path(temporary) / 'libgfortran.so.4'
        library.write_bytes(b'fixture')

        def missing(_, env):
            return [] if env.get('LD_PRELOAD') == str(library) else ['libgfortran.so.4']

        original = {'PATH': '/site/bin', 'LD_LIBRARY_PATH': '/site/current/lib'}
        with patch('kinbot.anl.runtime.shutil.which', return_value='/site/bin/xcfour'), \
                patch('kinbot.anl.runtime._missing_libraries', side_effect=missing), \
                patch('kinbot.anl.runtime._runtime_candidates', return_value=[library]):
            child, provenance = qc_runtime_environment('xcfour', 'cfour', original)
        assert original == {'PATH': '/site/bin', 'LD_LIBRARY_PATH': '/site/current/lib'}
        assert child['LD_LIBRARY_PATH'] == '/site/current/lib'
        assert child['LD_PRELOAD'] == str(library)
        assert provenance == {'runtime_library': str(library.resolve())}


def test_cfour_uses_module_runtime_when_available():
    with patch('kinbot.anl.runtime.shutil.which', return_value='/site/bin/xcfour'), \
            patch('kinbot.anl.runtime._missing_libraries', return_value=[]):
        child, provenance = qc_runtime_environment('xcfour', 'cfour', {'PATH': '/site/bin'})
    assert child == {'PATH': '/site/bin'}
    assert provenance == {}


def test_cfour_rejects_unresolved_runtime_before_launch():
    with TemporaryDirectory() as temporary:
        library = Path(temporary) / 'libgfortran.so.4'
        library.write_bytes(b'fixture')
        with patch('kinbot.anl.runtime.shutil.which', return_value='/site/bin/xcfour'), \
                patch('kinbot.anl.runtime._missing_libraries',
                      return_value=['libgfortran.so.4']), \
                patch('kinbot.anl.runtime._runtime_candidates', return_value=[library]):
            with pytest.raises(RuntimeError, match='CFOUR needs libgfortran.so.4'):
                qc_runtime_environment('xcfour', 'cfour', {'PATH': '/site/bin'})


def test_single_node_molpro_defaults_to_shared_memory_mpi():
    original = {'PATH': '/site/bin', 'SLURM_NNODES': '1'}
    child, provenance = qc_runtime_environment('molpro', 'molpro', original)
    assert original == {'PATH': '/site/bin', 'SLURM_NNODES': '1'}
    assert child['I_MPI_FABRICS'] == 'shm'
    assert provenance == {'mpi_fabrics': 'shm'}


def test_molpro_preserves_explicit_fabric_and_multinode_environment():
    child, provenance = qc_runtime_environment(
        'molpro', 'molpro',
        {'SLURM_NNODES': '1', 'I_MPI_FABRICS': 'ofi'})
    assert child['I_MPI_FABRICS'] == 'ofi'
    assert provenance == {'mpi_fabrics': 'ofi'}

    child, provenance = qc_runtime_environment(
        'molpro', 'molpro', {'SLURM_NNODES': '2'})
    assert 'I_MPI_FABRICS' not in child
    assert provenance == {}


def test_molpro_never_falls_back_to_home_when_tmp_is_too_small():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        home = root / 'home'
        work = root / 'work'
        tmp = root / 'tmp'
        for path in (home, work, tmp):
            path.mkdir()

        def disk_usage(path):
            free_mb = 100 if Path(path).resolve() == tmp.resolve() else 10000
            unit = 1024 * 1024
            return shutil._ntuple_diskusage(
                20000 * unit, (20000 - free_mb) * unit, free_mb * unit)

        env = {'HOME': str(home), 'TMPDIR': str(tmp), 'SLURM_NNODES': '1'}
        with patch('kinbot.anl.runtime.shutil.disk_usage', side_effect=disk_usage):
            child, provenance = qc_runtime_environment(
                'molpro', 'molpro', env, work_directory=work,
                scratch_min_mb=1000)
        scratch = Path(child['TMPDIR'])
        assert scratch.is_dir()
        assert not scratch.is_relative_to(home.resolve())
        assert not (home / '.cache').exists()
        assert provenance['scratch_source'] == 'system_tmp'
        assert provenance['scratch_minimum_mb'] == 1000
        cleanup_qc_runtime(provenance)
        assert not scratch.exists()


def test_environment_can_raise_but_not_lower_a_molpro_task_disk_requirement(tmp_path):
    from kinbot.anl.runtime import _molpro_scratch
    unit = 1024 * 1024
    def capacity(root):
        return shutil._ntuple_diskusage(10000*unit, 5000*unit, 5000*unit)
    environment = {'KINBOT_MOLPRO_SCRATCH': str(tmp_path),
                   'KINBOT_MOLPRO_SCRATCH_MIN_MB': '6000'}
    with patch('kinbot.anl.runtime.shutil.disk_usage', side_effect=capacity):
        with pytest.raises(RuntimeError, match='6000 MB'):
            _molpro_scratch(environment, tmp_path, 1000)
        environment['KINBOT_MOLPRO_SCRATCH_MIN_MB'] = '1'
        with pytest.raises(RuntimeError, match='7000 MB'):
            _molpro_scratch(environment, tmp_path, 7000)
        environment['KINBOT_MOLPRO_SCRATCH_MIN_MB'] = '1.5'
        with pytest.raises(ValueError, match='positive integer'):
            _molpro_scratch(environment, tmp_path, 1000)


@pytest.mark.parametrize('failure', ['directory_creation', 'parent_permissions', 'fsync', 'readback'])
def test_molpro_rejects_broken_scratch_and_uses_only_configured_fallback(tmp_path, monkeypatch, failure):
    from kinbot.anl import runtime
    primary, fallback = tmp_path / 'primary', tmp_path / 'fallback'
    env = {'KINBOT_MOLPRO_SCRATCH': str(primary),
           'KINBOT_MOLPRO_SCRATCH_FALLBACKS': str(fallback), 'SLURM_NNODES': '1'}
    monkeypatch.setattr(runtime.shutil, 'disk_usage', lambda _: shutil._ntuple_diskusage(
        20000 * 1024**2, 0, 20000 * 1024**2))
    original_mkdtemp = runtime.tempfile.mkdtemp
    original_check = runtime._check_scratch_io
    original_mkdir = Path.mkdir
    def mkdir(path, *args, **kwargs):
        if failure == 'parent_permissions' and path == primary:
            raise PermissionError(errno.EACCES, 'Permission denied', str(path))
        return original_mkdir(path, *args, **kwargs)
    def create(*args, **kwargs):
        if failure == 'directory_creation' and Path(kwargs['dir']) == primary:
            raise OSError(errno.EIO, 'Input/output error')
        return original_mkdtemp(*args, **kwargs)
    def check(directory):
        if directory.parent == primary:
            if failure == 'fsync':
                with patch.object(runtime.os, 'fsync', side_effect=OSError(errno.EIO, 'Input/output error')):
                    return original_check(directory)
            if failure == 'readback':
                with patch.object(Path, 'read_bytes', return_value=b'corrupt data'):
                    return original_check(directory)
        return original_check(directory)
    monkeypatch.setattr(runtime.tempfile, 'mkdtemp', create)
    monkeypatch.setattr(runtime, '_check_scratch_io', check)
    monkeypatch.setattr(Path, 'mkdir', mkdir)
    child, provenance = qc_runtime_environment('molpro', 'molpro', env,
        work_directory=tmp_path, scratch_min_mb=1000)
    try:
        scratch = Path(child['TMPDIR'])
        assert scratch.parent == fallback
        assert list(scratch.iterdir()) == []
        assert not primary.exists() or list(primary.iterdir()) == []
        assert provenance['scratch_io_checked'] is True
        assert provenance['scratch_host']
        assert len(provenance['scratch_rejected_candidates']) == 1
        assert str(primary) in provenance['scratch_rejected_candidates'][0]
        assert 'TMPDIR' not in env
    finally:
        cleanup_qc_runtime(provenance)


def test_explicit_scratch_failure_does_not_silently_choose_another_filesystem(tmp_path, monkeypatch):
    from kinbot.anl import runtime
    primary = tmp_path / 'primary'
    monkeypatch.setattr(runtime.shutil, 'disk_usage', lambda _: shutil._ntuple_diskusage(
        20000 * 1024**2, 0, 20000 * 1024**2))
    def failed(directory):
        raise OSError(errno.EIO, 'Input/output error')
    monkeypatch.setattr(runtime, '_check_scratch_io', failed)
    with pytest.raises(RuntimeError, match='Input/output error'):
        runtime._molpro_scratch({'KINBOT_MOLPRO_SCRATCH': str(primary),
                               'HOME': str(tmp_path / 'home')}, tmp_path, 1000)
    assert list(primary.iterdir()) == []
    assert not (tmp_path / 'home').exists()


def test_scratch_fallback_cannot_weaken_the_capacity_floor(tmp_path, monkeypatch):
    from kinbot.anl import runtime
    primary, fallback = tmp_path / 'primary', tmp_path / 'fallback'
    monkeypatch.setattr(runtime.shutil, 'disk_usage', lambda _: shutil._ntuple_diskusage(
        20000 * 1024**2, 19900 * 1024**2, 100 * 1024**2))
    with pytest.raises(RuntimeError, match='1000 MB'):
        runtime._molpro_scratch({'KINBOT_MOLPRO_SCRATCH': str(primary),
            'KINBOT_MOLPRO_SCRATCH_FALLBACKS': str(fallback)}, tmp_path, 1000)
    assert list(primary.iterdir()) == list(fallback.iterdir()) == []


def test_empty_scratch_fallback_path_is_rejected(tmp_path):
    from kinbot.anl.runtime import _molpro_scratch
    with pytest.raises(ValueError, match='empty path'):
        _molpro_scratch({'KINBOT_MOLPRO_SCRATCH': str(tmp_path),
                        'KINBOT_MOLPRO_SCRATCH_FALLBACKS': ':/another'}, tmp_path, 1000)


def test_failure_capture_retains_rank_evidence_before_scratch_cleanup(tmp_path):
    import json
    from kinbot.anl.runtime import preserve_qc_failure
    scratch = tmp_path/'scratch'
    work = tmp_path/'work'
    scratch.mkdir(); work.mkdir()
    (scratch/'rank.log.2').write_text('x'*70000 + '\nError writing 1024 words\n')
    (scratch/'huge.wfu').write_bytes(b'not a diagnostic')
    (scratch/'outside.log').symlink_to(work/'secret')
    (work/'secret').write_text('unrelated')
    (work/'step.stderr').write_text('MPI_Abort\n')
    provenance = {'scratch_directory': str(scratch)}
    path = preserve_qc_failure(provenance, work, 'step', 255)
    cleanup_qc_runtime(provenance)
    data = json.loads(path.read_text())
    assert not scratch.exists()
    assert data['explicit_io_failure'] is True
    assert set(data['files']) == {'scratch/rank.log.2', 'work/step.stderr'}
    assert len(data['files']['scratch/rank.log.2']['tail']) == 65536
    assert data['returncode'] == 255


def test_failure_capture_does_not_call_mpi_abort_an_io_error(tmp_path):
    import json
    from kinbot.anl.runtime import preserve_qc_failure
    (tmp_path/'step.stderr').write_text('application called MPI_Abort; KILLED BY SIGNAL: 9')
    path = preserve_qc_failure({}, tmp_path, 'step', 255)
    assert json.loads(path.read_text())['explicit_io_failure'] is False
