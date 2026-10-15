"""Real batch-body retries, checkpoint rollback and controller handoff."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pytest

from kinbot.slurm_retry import io_failure, retry_limit, run, wrap_script


@pytest.mark.parametrize('text', [
    'mktemp: cannot create directory: Input/output error',
    'OSError: [Errno 5] Input/output error',
    'sqlite3.OperationalError: disk I/O error',
    'ERROR WRITING 32768 WORDS AT OFFSET 569671680. TO FILE21',
    '? I/O error. Perhaps full disk?\n? The problem occurs in writew',
    'No space left on device', 'Stale file handle',
    'Erroneous write. Write -1 instead of 4096.',
])
def test_explicit_io_failures(text):
    assert io_failure(text)


@pytest.mark.parametrize('text', ['SCF did not converge', 'Permission denied',
    'String overflow inserting', 'I/O error handling enabled', 'Basis not found'])
def test_chemistry_and_other_failures_do_not_retry(text):
    assert not io_failure(text)


def test_retry_limit_is_bounded_and_can_be_disabled():
    assert retry_limit({}) == 2
    assert retry_limit({'KINBOT_SLURM_IO_RETRIES': '0'}) == 0
    for value in ('-1', '11', '1.5', 'yes'):
        with pytest.raises(ValueError, match='0 to 10'):
            retry_limit({'KINBOT_SLURM_IO_RETRIES': value})


def test_directory_retry_does_not_copy_or_rollback_archived_attempts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('SLURM_JOB_ID', '43')
    old = tmp_path / 'attempts/1/old.wfu'
    old.parent.mkdir(parents=True)
    old.write_bytes(b'previous failed native evidence')
    Path('job.inp').write_text('immutable input')
    command = [sys.executable, '-c', 'from pathlib import Path; Path("job.out").write_text("done")']
    assert run(command, identity='job', directories=[tmp_path]) == 0
    record = next(tmp_path.glob('.kinbot_slurm_io_retries/*/43/execution.json'))
    assert set(json.loads(record.read_text())['initial_artifacts']) == {'job.inp'}
    assert old.read_bytes() == b'previous failed native evidence'
    assert not (record.parent / 'initial/attempts').exists()


def test_requeued_allocation_archives_partial_products_and_restores_original_seed(tmp_path, monkeypatch):
    import hashlib
    import fcntl
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('SLURM_JOB_ID', '42')
    monkeypatch.setenv('SLURM_RESTART_COUNT', '0')
    Path('job.inp').write_text('file,2,kb_seed.wfu\n')
    Path('kb_seed.wfu').write_text('valid seed')
    command = [sys.executable, '-c',
               'from pathlib import Path; assert Path("kb_seed.wfu").read_text()=="valid seed"; '
               'assert not Path("job.out").exists(); Path("job.out").write_text("complete")']
    def interrupted(*args, **kwargs):
        Path('kb_seed.wfu').write_text('partial corrupt seed')
        Path('job.out').write_text('partial output')
        raise KeyboardInterrupt
    real_run = subprocess.run
    monkeypatch.setattr('kinbot.slurm_retry.subprocess.run', interrupted)
    with pytest.raises(KeyboardInterrupt):
        run(command, identity='job', prefixes=['job'])
    token = hashlib.sha256(b'job').hexdigest()[:20]
    folder = Path('.kinbot_slurm_io_retries') / token / '42'
    original_record = (folder / 'execution.json').read_bytes()
    monkeypatch.setattr('kinbot.slurm_retry.subprocess.run', real_run)
    monkeypatch.setenv('SLURM_RESTART_COUNT', '1')
    with (folder / 'lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):
            run(command, identity='job', prefixes=['job'])
        assert Path('job.out').read_text() == 'partial output'
    assert run(command, identity='job', prefixes=['job']) == 0
    assert (folder / 'execution.json').read_bytes() == original_record
    assert (folder / 'restart-1/interrupted/kb_seed.wfu').read_text() == 'partial corrupt seed'
    assert (folder / 'restart-1/interrupted/job.out').read_text() == 'partial output'
    assert Path('job.out').read_text() == 'complete'
    assert json.loads((folder / 'restart-1/execution.json').read_text())['random_seed'] == json.loads(original_record)['random_seed']
    with pytest.raises(RuntimeError, match='duplicate execution'):
        run(command, identity='job', prefixes=['job'])
    monkeypatch.setenv('SLURM_RESTART_COUNT', '2')
    (folder / 'initial/kb_seed.wfu').write_text('tampered')
    with pytest.raises(RuntimeError, match='initial artifacts changed'):
        run(command, identity='job', prefixes=['job'])


@pytest.mark.parametrize('checkpoint', ['job.wfu', 'kb_1234.wfu'])
def test_batch_retries_io_and_restores_exact_checkpoint_without_touching_other_jobs(tmp_path, monkeypatch, checkpoint):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('SLURM_JOB_ID', '42')
    Path(checkpoint).write_bytes(b'verified starting orbitals')
    Path('job.inp').write_text(f'file,2,{checkpoint}\n')
    Path('job-other.out').write_text('another active job')
    Path('counter').write_text('0')
    program = (
        'from pathlib import Path\nimport sys\n'
        'counter=Path("counter"); n=int(counter.read_text())+1; counter.write_text(str(n))\n'
        'assert Path("job.wfu").read_bytes()==b"verified starting orbitals"\n'
        'if n == 1:\n'
        ' Path("job.wfu").write_bytes(b"partial corrupt orbitals")\n'
        ' Path("job.out").write_text("ERROR WRITING 32768 WORDS AT OFFSET 0. TO FILE21")\n'
        ' sys.exit(1)\n'
        'assert not Path("job.out").exists()\n'
        'Path("job.out").write_text("normal termination")\n')
    program = program.replace('job.wfu', checkpoint)
    assert run([sys.executable, '-c', program], identity='job', prefixes=['job'],
               exit_marker='job.exitcode') == 0
    journal = next(Path('.kinbot_slurm_io_retries').glob('*/42/execution.json'))
    record = json.loads(journal.read_text())
    assert record['status'] == 'complete' and len(record['attempts']) == 2
    import hashlib
    assert record['initial_artifacts'][checkpoint] == hashlib.sha256(b'verified starting orbitals').hexdigest()
    assert record['attempts'][0]['io_failure'] is True
    assert Path('job.exitcode').read_text() == '42 0\n'
    assert Path('job-other.out').read_text() == 'another active job'
    assert (journal.parent / 'attempt-1/artifacts' / checkpoint).read_bytes() == b'partial corrupt orbitals'


@pytest.mark.parametrize('message,limit,attempts', [
    ('Input/output error', '2', 3), ('SCF failed', '2', 1), ('Input/output error', '0', 1)])
def test_retries_are_bounded_and_final_failure_remains_visible(tmp_path, monkeypatch, message, limit, attempts):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('KINBOT_SLURM_IO_RETRIES', limit)
    command = [sys.executable, '-c', f'import sys; sys.stderr.write({message!r}); sys.exit(9)']
    assert run(command, identity='test', exit_marker='test.exitcode') == 9
    record = json.loads(next(Path('.kinbot_slurm_io_retries').glob('*/*/execution.json')).read_text())
    assert record['status'] == 'failed' and len(record['attempts']) == attempts
    assert Path('test.exitcode').read_text().endswith(' 9\n')


def test_real_slurm_guard_keeps_directives_and_does_not_recurse(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    original = '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks=8\n' + shlex.join(
        [sys.executable, '-c', 'from pathlib import Path; Path("called").write_text("once")']) + '\n'
    wrapped = wrap_script(original, identity='sample_{surf_id}', template=True)
    rendered = wrapped.format(surf_id=7)
    assert rendered.startswith('#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks=8\n')
    assert wrap_script(rendered, identity='sample_7') == rendered
    path = Path('sample.sh')
    path.write_text(rendered)
    subprocess.run(['bash', '-n', str(path)], check=True)
    subprocess.run(['bash', str(path)], check=True, timeout=20)
    assert Path('called').read_text() == 'once'


def test_failed_qc_row_can_trigger_retry_despite_zero_python_exit(tmp_path, monkeypatch):
    from ase import Atoms
    from ase.db import connect
    monkeypatch.chdir(tmp_path)
    connect('kinbot.db').write(Atoms('H'), name='job', data={'status': 'error'})
    command = [sys.executable, '-c', 'from pathlib import Path; Path("job.log").write_text("Input/output error")']
    assert run(command, identity='job', prefixes=['job'], qc_job='job') == 1
    assert Path('kinbot.db').is_file()
    record = json.loads(next(Path('.kinbot_slurm_io_retries').glob('*/*/execution.json')).read_text())
    assert len(record['attempts']) == 3


def test_current_normal_qc_pickle_takes_precedence_over_historical_error_row(tmp_path, monkeypatch):
    from ase import Atoms
    from ase.db import connect
    monkeypatch.chdir(tmp_path)
    connect('kinbot.db').write(Atoms('H'), name='job', data={'status': 'error'})
    command = [sys.executable, '-c',
               'import pickle; pickle.dump({"data":{"status":"normal"}}, open("job.pkl", "wb"))']
    assert run(command, identity='job', prefixes=['job'], qc_job='job') == 0
    record = json.loads(next(Path('.kinbot_slurm_io_retries').glob('*/*/execution.json')).read_text())
    assert len(record['attempts']) == 1


def test_ordinary_qc_batch_defers_marker_until_io_retry_finishes(tmp_path, monkeypatch):
    from kinbot import kb_path
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('SLURM_JOB_ID', '731')
    Path('counter').write_text('0')
    Path('job.py').write_text(
        'from pathlib import Path\nimport pickle\n'
        'assert not Path("job.exitcode").exists()\n'
        'p=Path("counter"); n=int(p.read_text())+1; p.write_text(str(n))\n'
        'Path("job.log").write_text("Input/output error" if n == 1 else "normal termination")\n'
        'pickle.dump({"data":{"status":"error" if n == 1 else "normal"}}, open("job.pkl", "wb"))\n')
    body = (Path(kb_path) / 'tpl/slurm_python.tpl').read_text().format(
        name='job', python_executable=shlex.quote(sys.executable), python_file='job.py', arguments='')
    Path('job.sbatch').write_text(wrap_script('#!/bin/bash\n#SBATCH --ntasks=1\n' + body,
        identity='job', prefixes=['job'], exit_marker='job.exitcode', qc_job='job'))
    subprocess.run(['bash', 'job.sbatch'], check=True, timeout=20)
    assert Path('counter').read_text() == '2'
    assert Path('job.exitcode').read_text() == '731 0\n'


def test_scratch_setup_exception_is_classified_without_native_output(tmp_path, monkeypatch):
    from kinbot.anl import runtime
    import errno
    monkeypatch.setenv('KINBOT_MOLPRO_SCRATCH', str(tmp_path))
    monkeypatch.setattr(runtime, '_check_scratch_io', lambda _: (_ for _ in ()).throw(OSError(errno.EIO, 'broken')))
    with pytest.raises(runtime.MolproScratchError) as error:
        runtime.qc_runtime_environment('molpro', 'molpro', work_directory=tmp_path, scratch_min_mb=1)
    assert error.value.retryable_io


def test_retry_prefers_another_allowed_scratch_root_and_records_provenance(tmp_path, monkeypatch):
    from kinbot.anl.runtime import qc_runtime_environment, cleanup_qc_runtime
    monkeypatch.chdir(tmp_path)
    initial, alternate = tmp_path / 'first', tmp_path / 'second'
    monkeypatch.setenv('KINBOT_MOLPRO_SCRATCH', str(initial))
    monkeypatch.setenv('KINBOT_MOLPRO_SCRATCH_FALLBACKS', str(alternate))
    roots = []
    def native(command, **kwargs):
        env, provenance = qc_runtime_environment('molpro', 'molpro', kwargs['env'],
            work_directory=tmp_path, scratch_min_mb=1)
        roots.append(provenance['scratch_root'])
        cleanup_qc_runtime(provenance)
        if len(roots) == 1:
            kwargs['stderr'].write('mktemp: Input/output error\n')
        from types import SimpleNamespace
        return SimpleNamespace(returncode=1 if len(roots) == 1 else 0)
    monkeypatch.setattr('kinbot.slurm_retry.subprocess.run', native)
    assert run(['native'], identity='test') == 0
    assert roots == [str(initial), str(alternate)]


class _AtomicFlux:
    draws = []
    fail = True

    def run(self, length, face, sample):
        import numpy as np
        self.draws.append(float(np.random.random()))
        if self.fail:
            raise OSError(5, 'Input/output error')


def test_rotd_retry_preserves_running_pickle_and_random_sequence(tmp_path, monkeypatch):
    from enum import IntEnum
    import pickle
    from types import SimpleNamespace
    from kinbot import rotdpy_worker
    class Tag(IntEnum):
        FLUX_TAG = 1
        SURF_TAG = 2
        STOP_TAG = 3
    monkeypatch.setitem(sys.modules, 'rotd_py.system', SimpleNamespace(FluxTag=Tag))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('KINBOT_SLURM_IO_RETRY_SEED', '1234')
    monkeypatch.setattr(_AtomicFlux, 'draws', [])
    monkeypatch.setattr(_AtomicFlux, 'fail', True)
    path = Path('sample.pkl')
    original = pickle.dumps([1, _AtomicFlux(), 0, 0, 10, 0, 'RUNNING'])
    path.write_bytes(original)
    with pytest.raises(OSError, match='Input/output error'):
        rotdpy_worker.run('sample')
    assert path.read_bytes() == original
    monkeypatch.setattr(_AtomicFlux, 'fail', False)
    rotdpy_worker.run('sample')
    assert _AtomicFlux.draws[0] == _AtomicFlux.draws[1]
    assert pickle.loads(path.read_bytes())[-1] == 'COMPLETED'
    assert list(tmp_path.glob('*.tmp')) == []


def test_rotd_failed_completion_write_never_publishes_partial_pickle(tmp_path, monkeypatch):
    import pickle
    from enum import IntEnum
    from types import SimpleNamespace
    from kinbot import rotdpy_worker
    class Tag(IntEnum):
        FLUX_TAG = 1
        SURF_TAG = 2
        STOP_TAG = 3
    monkeypatch.setitem(sys.modules, 'rotd_py.system', SimpleNamespace(FluxTag=Tag))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(_AtomicFlux, 'fail', False)
    path = Path('sample.pkl')
    original = pickle.dumps([1, _AtomicFlux(), 0, 0, 10, 0, 'RUNNING'])
    path.write_bytes(original)
    def failed_dump(payload, stream):
        stream.write(b'partial')
        raise OSError(5, 'Input/output error')
    monkeypatch.setattr(rotdpy_worker.pickle, 'dump', failed_dump)
    with pytest.raises(OSError):
        rotdpy_worker.run('sample')
    assert path.read_bytes() == original
    assert list(tmp_path.glob('*.tmp')) == []


def test_existing_rotd_template_acquires_guard_and_atomic_worker(tmp_path):
    from kinbot.rotdpy import _guard_slurm_template
    path = tmp_path / 'qu.tpl'
    path.write_text('#!/bin/bash\n#SBATCH --ntasks={procs}\n'
                    '/usr/bin/python surf{surf_id}_face{face_id}_samp{samp_id}.py\n')
    _guard_slurm_template(path)
    once = path.read_text()
    _guard_slurm_template(path)
    assert once == path.read_text()
    rendered = once.format(procs=8, surf_id=0, face_id=1, samp_id=2)
    assert '-m kinbot.slurm_retry run' in rendered
    assert '/usr/bin/python -m kinbot.rotdpy_worker surf0_face1_samp2' in rendered
    subprocess.run(['bash', '-n'], input=rendered, text=True, check=True)
