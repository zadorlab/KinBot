"""Shared, bounded I/O retries for KinBot-generated Slurm batch scripts.

Retries stay within the allocation, so they cannot bypass queue admission or
node caps. Only explicitly owned artifacts are rolled back; shared databases
and unrelated jobs are never moved. Each attempt retains its native evidence.
"""
from argparse import ArgumentParser
import fcntl
import hashlib
import json
import os
from pathlib import Path
import pickle
import re
import secrets
import shlex
import shutil
import socket
import stat
import sqlite3
import subprocess
import sys

from kinbot.execution import atomic_json, file_sha256


_MARKER = '# KinBot shared Slurm I/O retry guard v1'
_BINARY = {'.wfu', '.chk', '.fchk', '.pkl', '.db', '.sqlite', '.traj', '.npy', '.npz'}


def retry_limit(env=None):
    value = (os.environ if env is None else env).get('KINBOT_SLURM_IO_RETRIES', '2')
    if not re.fullmatch(r'[0-9]+', str(value)) or not 0 <= int(value) <= 10:
        raise ValueError('KINBOT_SLURM_IO_RETRIES must be an integer from 0 to 10.')
    return int(value)


def io_failure(text):
    """Match actual OS/native I/O diagnostics, not an arbitrary failed job."""
    lower = text.lower()
    if any(marker in lower for marker in (
            'input/output error', 'stale file handle', 'no space left on device',
            'disk quota exceeded', 'sqlite3.operationalerror: disk i/o error')):
        return True
    if re.search(r'(?im)^\s*error (?:writing|reading)\s+\d+\s+words\b', text):
        return True
    if re.search(r'(?im)^\s*erroneous (?:write|read)\.', text):
        return True
    return bool(re.search(r'(?im)^\s*\?\s*i/o error\.', text)
                and re.search(r'(?i)problem occurs in (?:writew|readw)\b', text))


def wrap_script(contents, *, identity, prefixes=(), directories=(), exit_marker=None,
                qc_job=None, python=None, template=False):
    """Keep all batch directives and the original body, adding one shared guard.

    ``template`` preserves a later format pass used by ROTD_py. No custom
    chemistry, resource directive or site setup is replaced.
    """
    if _MARKER in contents:
        return contents
    lines = contents.splitlines(keepends=True)
    index = 0
    while index < len(lines) and (not lines[index].strip() or lines[index].lstrip().startswith('#')):
        index += 1
    command = [python or sys.executable, '-m', 'kinbot.slurm_retry', 'run', '--identity', identity]
    token = hashlib.sha256(identity.encode()).hexdigest()[:20]
    command += ['--guard-token', token]
    for path in prefixes:
        command += ['--prefix', str(path)]
    for path in directories:
        command += ['--directory', str(path)]
    if exit_marker:
        command += ['--exit-marker', str(exit_marker)]
    if qc_job:
        command += ['--qc-job', qc_job]
    guard = (_MARKER + '\n'
             f'if [ "${{KINBOT_SLURM_RETRY_CHILD:-}}" != {shlex.quote(token)} ]; then\n'
             f'    exec {shlex.join(command)} -- bash "$0" "$@"\n'
             'fi\n')
    if template:
        # Escape only shell parameter braces; identity/scope placeholders
        # deliberately remain for ROTD's surf/face/sample format pass.
        guard = guard.replace('${KINBOT_SLURM_RETRY_CHILD:-}', '${{KINBOT_SLURM_RETRY_CHILD:-}}')
    return ''.join(lines[:index]) + guard + ''.join(lines[index:])


def _artifacts(prefixes, directories, journal):
    scheduler_streams = set()
    for descriptor in (1, 2):
        metadata = os.fstat(descriptor)
        if stat.S_ISREG(metadata.st_mode):
            scheduler_streams.add((metadata.st_dev, metadata.st_ino))
    paths = set()
    for directory in directories:
        directory = Path(directory).resolve()
        paths.update(path for path in directory.rglob('*')
                     if path.is_file() and not path.is_symlink() and journal not in path.parents)
    for prefix in prefixes:
        prefix = Path(prefix).absolute()
        if prefix.parent.is_dir():
            paths.update(path for path in prefix.parent.iterdir()
                         if path.is_file() and not path.is_symlink()
                         and (path.name.startswith(prefix.name + '.')
                              or re.match(re.escape(prefix.name) + r'_(?:step|basis)_\d+(?:[._]|$)', path.name)
                              or path.name in (prefix.name + '_sella.log', prefix.name + '_freq.out')))
    # Molpro FILE names must be short: hashed checkpoints need not share the
    # job prefix. Follow only local FILE,2 WFU basenames in this job's inputs
    # so an I/O retry restores its seeds instead of reusing corrupted orbitals.
    for inp in tuple(paths):
        if inp.suffix != '.inp':
            continue
        names = re.findall(r'(?im)^\s*file\s*,\s*2\s*,\s*([a-z0-9][a-z0-9_.-]*\.wfu)\s*$',
                           inp.read_text(errors='replace'))
        for name in names:
            checkpoint = inp.parent / name
            if len(name) <= 31 and checkpoint.is_file() and not checkpoint.is_symlink():
                paths.add(checkpoint)
    # A shared calculation database is transactional and belongs to the
    # controller. Native/job-specific checkpoints are covered by the scopes.
    return {path.resolve() for path in paths
            if path.name not in ('kinbot.db', 'kinbot.db-wal', 'kinbot.db-shm', 'kinbot.db-journal')
            and journal not in path.resolve().parents
            and (path.stat().st_dev, path.stat().st_ino) not in scheduler_streams}


def _text(path):
    if path.suffix in _BINARY:
        return ''
    try:
        with path.open('rb') as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - 262144))
            return stream.read().decode(errors='replace')
    except OSError:
        return ''


def _qc_failed(job):
    if not job:
        return False
    pending = Path(job + '.pkl')
    if pending.is_file():
        try:
            with pending.open('rb') as stream:
                payload = pickle.load(stream)
            return payload['data'].get('status') == 'error'
        except (EOFError, pickle.UnpicklingError, KeyError, TypeError):
            return True
    if not Path('kinbot.db').is_file():
        return False
    from ase.db import connect
    rows = list(connect('kinbot.db').select(name=job))
    return bool(rows and rows[-1].data.get('status') == 'error')


def run(command, *, identity, prefixes=(), directories=(), exit_marker=None, qc_job=None,
        guard_token=None):
    limit = retry_limit()
    root = Path.cwd().resolve()
    journal = root / '.kinbot_slurm_io_retries'
    token = hashlib.sha256(identity.encode()).hexdigest()[:20]
    allocation = re.sub(r'[^A-Za-z0-9_.-]', '_', os.environ.get('SLURM_JOB_ID', f'local-{os.getpid()}'))
    allocation_folder = journal / token / allocation
    allocation_folder.mkdir(parents=True, exist_ok=True)
    restart = os.environ.get('SLURM_RESTART_COUNT', '0')
    if not re.fullmatch(r'[0-9]+', restart):
        raise ValueError('SLURM_RESTART_COUNT must be a nonnegative integer.')
    restart = int(restart)
    folder = allocation_folder if not restart else allocation_folder / f'restart-{restart}'
    original_seed = None
    # A Slurm requeue keeps its job ID. Serialize all incarnations of that
    # allocation, while retaining a separate journal for each restart count.
    with (allocation_folder / 'lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        folder.mkdir(exist_ok=True)
        if (folder / 'execution.json').exists():
            raise RuntimeError('This Slurm allocation already has a retry journal; refusing duplicate execution.')
        if restart and (allocation_folder / 'execution.json').is_file():
            original = json.loads((allocation_folder / 'execution.json').read_text())
            if original['identity'] != identity or original['command'] != list(command):
                raise RuntimeError('Slurm requeue command differs from the original allocation.')
            original_seed = original.get('random_seed')
            restored = {}
            for relative, digest in original['initial_artifacts'].items():
                target = (root / relative).resolve()
                saved = allocation_folder / 'initial' / relative
                if root not in target.parents or not saved.is_file() or file_sha256(saved) != digest:
                    raise RuntimeError('Slurm requeue initial artifacts changed; refusing restart.')
                restored[target] = saved
            # Preserve interrupted products, then reset to the pre-allocation
            # inputs/checkpoints. A partial WFU is never a new starting guess.
            for current in _artifacts(prefixes, directories, journal):
                if root not in current.parents:
                    raise ValueError('Slurm retry artifact scopes must lie inside the job working tree.')
                if current in restored and file_sha256(current) == original['initial_artifacts'][str(current.relative_to(root))]:
                    continue
                destination = folder / 'interrupted' / current.relative_to(root)
                destination.parent.mkdir(parents=True, exist_ok=True)
                current.rename(destination)
            for target, saved in restored.items():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(saved, target)
        seeds = _artifacts(prefixes, directories, journal)
        # Scopes must be inside the submission working tree; never roll back
        # external data, a library, or another run's absolute checkpoint.
        if any(root not in path.parents for path in seeds):
            raise ValueError('Slurm retry artifact scopes must lie inside the job working tree.')
        baseline = {}
        for path in seeds:
            relative = path.relative_to(root)
            saved = folder / 'initial' / relative
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, saved)
            baseline[path] = {'saved': str(saved), 'sha256': file_sha256(saved)}
        env = dict(os.environ, KINBOT_SLURM_RETRY_CHILD=guard_token or token,
                   KINBOT_SLURM_RETRY_ACTIVE='1')
        if original_seed is not None:
            env['KINBOT_SLURM_IO_RETRY_SEED'] = original_seed
        else:
            env.setdefault('KINBOT_SLURM_IO_RETRY_SEED', str(secrets.randbits(32)))
        record = {'schema': 1, 'identity': identity, 'command': list(command),
                  'host': socket.gethostname(), 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
                  'slurm_restart_count': restart,
                  'initial_artifacts': {str(path.relative_to(root)): seed['sha256']
                                        for path, seed in baseline.items()},
                  'retry_limit': limit, 'random_seed': env['KINBOT_SLURM_IO_RETRY_SEED'],
                  'status': 'running', 'attempts': []}
        record_path = folder / 'execution.json'
        avoided = []
        for number in range(limit + 1):
            attempt = folder / f'attempt-{number + 1}'
            attempt.mkdir()
            env['KINBOT_SLURM_IO_ATTEMPT_DIRECTORY'] = str(attempt)
            stdout_path, stderr_path = attempt / 'stdout.txt', attempt / 'stderr.txt'
            atomic_json(record_path, record)
            with stdout_path.open('w') as stdout, stderr_path.open('w') as stderr:
                try:
                    result = subprocess.run(command, env=env, stdout=stdout, stderr=stderr, check=False)
                    returncode = result.returncode
                except OSError as error:
                    stderr.write(f'{type(error).__name__}: {error}\n')
                    returncode = 127 if isinstance(error, FileNotFoundError) else 1
            probe_error = ''
            try:
                failed = bool(returncode or _qc_failed(qc_job))
            except (OSError, sqlite3.DatabaseError) as error:
                failed = True
                probe_error = f'{type(error).__name__}: {error}'
            current = _artifacts(prefixes, directories, journal)
            changed = {path for path in current if path not in baseline
                       or file_sha256(path) != baseline[path]['sha256']}
            detail = _text(stderr_path) + '\n' + _text(stdout_path) + '\n' + probe_error
            detail += '\n' + '\n'.join(_text(path) for path in changed)
            retryable = failed and io_failure(detail)
            # Native Python workers can supply precise scratch-setup evidence.
            for path in changed:
                if path.name == 'execution.json':
                    try:
                        native = json.loads(path.read_text())
                    except (OSError, ValueError):
                        continue
                    retryable |= failed and native.get('retryable_io') is True
                    scratch = native.get('runtime', {}).get('scratch_root')
                    if scratch and scratch not in avoided:
                        avoided.append(scratch)
            for provenance in attempt.glob('scratch-*.json'):
                scratch = json.loads(provenance.read_text()).get('scratch_root')
                if scratch and scratch not in avoided:
                    avoided.append(scratch)
            entry = {'number': number + 1, 'returncode': returncode,
                     'failed_result': failed, 'io_failure': bool(retryable),
                     'stdout': str(stdout_path), 'stderr': str(stderr_path)}
            record['attempts'].append(entry)
            retrying = retryable and number < limit
            if not retrying:
                record['status'] = 'failed' if failed else 'complete'
                atomic_json(record_path, record)
                if len(record['attempts']) > 1:
                    print(f'{identity}: {record["status"]} after {len(record["attempts"]) - 1} I/O '
                          f'retry/retries; evidence: {record_path}', file=sys.stderr, flush=True)
                for source, target in ((stdout_path, sys.stdout), (stderr_path, sys.stderr)):
                    tail = _text(source)
                    if tail:
                        target.write(tail)
                code = returncode if returncode else (1 if failed else 0)
                if exit_marker:
                    marker = Path(exit_marker)
                    marker.parent.mkdir(parents=True, exist_ok=True)
                    temporary = marker.with_name(marker.name + f'.{os.getpid()}.tmp')
                    temporary.write_text(f'{os.environ.get("SLURM_JOB_ID", "unknown")} {code}\n')
                    os.replace(temporary, marker)
                return code
            # Archive partial results, then restore the exact pre-job seeds.
            # In particular, a WFU/CHK written during failed I/O is not reused.
            for path in changed:
                destination = attempt / 'artifacts' / path.relative_to(root)
                destination.parent.mkdir(parents=True, exist_ok=True)
                path.rename(destination)
            for path, seed in baseline.items():
                saved = Path(seed['saved'])
                if file_sha256(saved) != seed['sha256']:
                    raise RuntimeError('Slurm retry initial checkpoint changed; refusing retry.')
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(saved, path)
            env['KINBOT_MOLPRO_AVOID_SCRATCH_ROOTS'] = json.dumps(avoided)
            atomic_json(record_path, record)


def main(argv=None):
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['run'])
    parser.add_argument('--identity', required=True)
    parser.add_argument('--prefix', action='append', default=[])
    parser.add_argument('--directory', action='append', default=[])
    parser.add_argument('--exit-marker')
    parser.add_argument('--qc-job')
    parser.add_argument('--guard-token')
    parser.add_argument('command', nargs='*')
    # argparse's REMAINDER would swallow the guard's own options. Split the
    # explicit command boundary instead, preserving shell argument quoting.
    arguments = list(sys.argv[1:] if argv is None else argv)
    boundary = arguments.index('--')
    args = parser.parse_args(arguments[:boundary])
    return run(arguments[boundary + 1:], identity=args.identity, prefixes=args.prefix,
               directories=args.directory, exit_marker=args.exit_marker, qc_job=args.qc_job,
               guard_token=args.guard_token)


if __name__ == '__main__':
    raise SystemExit(main())
