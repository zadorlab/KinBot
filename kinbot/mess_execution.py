"""Submit MESS calculations and distinguish completion from solver success."""
import logging
from pathlib import Path
import re
import shlex
import subprocess
import time

from kinbot import constants

logger = logging.getLogger('KinBot')
_RESULT_GRACE_SECONDS = 60.
_QUERY_TIMEOUT_SECONDS = 30.
_MAX_QUERY_FAILURES = 3


class MESSExecutionError(RuntimeError):
    """MESS did not produce a successful rate calculation."""


class _ResultPending(MESSExecutionError):
    """A completed queued job's result files are not yet visible."""


def _run(command, *, timeout=None):
    return subprocess.run(command, capture_output=True, text=True, check=False,
                          timeout=timeout)


def _submit(queue, script):
    result = _run([constants.qsubmit[queue], str(script)])
    if result.returncode:
        raise MESSExecutionError(f'MESS submission failed: {result.stderr.strip()}')
    if queue == 'slurm':
        match = re.search(r'\bSubmitted batch job (\d+)\b', result.stdout)
    else:
        match = re.fullmatch(r'\s*(\d+(?:\.[\w.-]+)?)\s*', result.stdout)
    if not match:
        raise MESSExecutionError(f'Cannot read MESS job ID from {result.stdout!r}.')
    return match.group(1)


def _active(queue, pid):
    command = (['squeue', '--noheader', '--jobs', pid, '--format=%T']
               if queue == 'slurm' else ['qstat', '-f', pid])
    try:
        result = _run(command, timeout=_QUERY_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired as error:
        raise MESSExecutionError(f'Cannot query MESS job {pid}: scheduler query timed out.') from error
    except OSError as error:
        raise MESSExecutionError(f'Cannot query MESS job {pid}: {error}') from error
    if queue == 'slurm':
        if result.returncode:
            if re.search(r'invalid job id', result.stderr, re.IGNORECASE):
                return False
            raise MESSExecutionError(f'Cannot query MESS job {pid}: {result.stderr.strip()}')
        states = result.stdout.split()
        terminal = {'COMPLETED', 'CANCELLED', 'FAILED', 'TIMEOUT', 'NODE_FAIL',
                    'OUT_OF_MEMORY', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}
        return any(state not in terminal for state in states)
    if result.returncode:
        if re.search(r'unknown job|unknown job id|job has finished|invalid job id',
                     result.stderr, re.IGNORECASE):
            return False
        raise MESSExecutionError(f'Cannot query MESS job {pid}: {result.stderr.strip()}')
    state = re.search(r'\bjob_state\s*=\s*(\w+)', result.stdout)
    if not state:
        raise MESSExecutionError(f'Cannot read state of MESS job {pid}.')
    return state.group(1) not in ('C', 'F')


def _stamp(path):
    try:
        stat = path.stat()
    except FileNotFoundError:
        return None
    return stat.st_mtime_ns, stat.st_size


def _verify(index, previous_output, exit_code=None):
    index = f'{index:04d}' if isinstance(index, int) else index
    stem = Path('me') / f'mess_{index}'
    if exit_code is None:
        status = stem.with_suffix('.exitcode')
        try:
            exit_code = int(status.read_text().strip())
        except (OSError, ValueError) as error:
            raise _ResultPending(
                f'MESS calculation {index} left the queue without a solver exit result; '
                'check the scheduler output for cancellation or failure.') from error
    if exit_code:
        raise MESSExecutionError(
            f'MESS calculation {index} failed with exit code {exit_code}; '
            'reaction generation may have completed, but rates were not obtained. '
            'Inspect the MESS log and scheduler output.')
    output = stem.with_suffix('.out')
    try:
        current_output = _stamp(output)
    except OSError as error:
        raise _ResultPending(f'MESS calculation {index}: cannot read rate output at {output}.') from error
    if current_output is None or not current_output[1] or current_output == previous_output:
        raise _ResultPending(
            f'MESS calculation {index} exited successfully but produced no new, '
            f'nonempty rate output at {output}.')
    logger.info('MESS calculation %s completed successfully: %s', index, output)


def run_mess(writer):
    """Write submission scripts, then run and verify every ready network/UQ calculation.

    Queue absence means only that execution ended. The script's exit result
    and newly written rate output establish whether MESS itself succeeded.
    """
    queue = writer.par['queuing']
    if queue not in ('local', 'slurm', 'pbs'):
        raise MESSExecutionError(f'Unsupported MESS queue: {queue}')
    scripts = []
    if hasattr(writer, 'mess_jobs'):
        indices = [job['stem'].removeprefix('mess_') for job in writer.mess_jobs
                   if job['status'] == 'ready']
    else:
        indices = [f'{index:04d}' for index in range(writer.par['uq_n'])]
    for index in indices:
        extension = '.sh' if queue == 'local' else constants.qext[queue]
        script = Path('me') / f'run_mess_{index}{extension}'
        if queue == 'local':
            script.write_text(f'#!/bin/sh\ncd me || exit 1\nexec mess mess_{index}.inp\n')
        else:
            writer.write_submitscript(str(script), index)
        scripts.append(script)
    command = 'sh' if queue == 'local' else constants.qsubmit[queue]
    batch = Path('batch_me.sub')
    batch.write_text(''.join(f'{command} {shlex.quote(str(script))}\n' for script in scripts))
    batch.chmod(0o700)
    if not writer.par['run_me']:
        return 0

    active = {}
    ended = {}
    query_failures = {}
    failures = []
    limit = max(1, int(writer.par['uq_max_runs']))

    def poll():
        time.sleep(5)
        for pid, (index, previous) in list(active.items()):
            if pid not in ended:
                try:
                    running = _active(queue, pid)
                except MESSExecutionError as error:
                    count = query_failures.get(pid, 0) + 1
                    query_failures[pid] = count
                    if count < _MAX_QUERY_FAILURES:
                        logger.warning('%s Retrying scheduler query (%s/%s).',
                                       error, count, _MAX_QUERY_FAILURES)
                        continue
                    raise MESSExecutionError(
                        f'{error} Monitoring stopped after {count} consecutive query failures; '
                        f'submitted jobs may still be active: {", ".join(active)}.') from error
                query_failures.pop(pid, None)
                if running:
                    continue
                ended[pid] = time.monotonic()
            try:
                _verify(index, previous)
            except _ResultPending as error:
                if time.monotonic() - ended[pid] < _RESULT_GRACE_SECONDS:
                    continue
                failures.append(f'{error} Result files did not become available within '
                                f'{_RESULT_GRACE_SECONDS:g} seconds after queue completion.')
            except MESSExecutionError as error:
                failures.append(str(error))
            del active[pid]
            ended.pop(pid, None)

    for index, script in zip(indices, scripts):
        output = Path('me') / f'mess_{index}.out'
        previous = _stamp(output)
        if queue == 'local':
            with open(f'me/mess_{index}.stdout', 'w') as stdout, \
                    open(f'me/mess_{index}.err', 'w') as stderr:
                result = subprocess.run(['sh', str(script)], stdout=stdout, stderr=stderr,
                                        check=False)
            _verify(index, previous, result.returncode)
        else:
            while len(active) >= limit:
                poll()
            if failures:
                break
            (Path('me') / f'mess_{index}.exitcode').unlink(missing_ok=True)
            try:
                pid = _submit(queue, script)
            except MESSExecutionError as error:
                failures.append(str(error))
                break
            active[pid] = (index, previous)
    while active:
        poll()
    if failures:
        raise MESSExecutionError('\n'.join(failures))
    return 0
