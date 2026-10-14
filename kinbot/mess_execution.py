"""Submit MESS calculations and distinguish completion from solver success."""
import hashlib
import json
import logging
import os
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
_SCHEDULER_UNAVAILABLE = re.compile(
    r'(?:unable to contact (?:slurm )?controller|'
    r'unable to establish controller connection|'
    r'(?:scheduler )?controller unavailable|'
    r'slurm_load_jobs error|socket timed out|connection refused|'
    r'resource temporarily unavailable|network is unreachable|'
    r'host is down|no route to host|'
    r'communication connection failure|controller not responding|'
    r'scheduler query timed out)', re.IGNORECASE)


class MESSExecutionError(RuntimeError):
    """MESS did not produce a successful rate calculation."""


class _ResultPending(MESSExecutionError):
    """A completed queued job's result files are not yet visible."""


class _SubmissionUncertain(MESSExecutionError):
    """A scheduler call may have accepted a job without returning its ID."""


def _scheduler_unavailable(error):
    return bool(_SCHEDULER_UNAVAILABLE.search(str(error)))


def _atomic_json(path, payload):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(
        payload, indent=2, sort_keys=True) + '\n')
    os.replace(temporary, path)


def mess_command(parameters):
    """Read an executable and optional arguments, without shell expansion."""
    value = (os.environ.get('KINBOT_MESS_COMMAND', '').strip()
             or parameters.get('mess_command', 'mess'))
    if not isinstance(value, str) or not value.strip():
        raise MESSExecutionError('mess_command must specify an executable and optional arguments.')
    try:
        command = shlex.split(value)
    except ValueError as error:
        raise MESSExecutionError(f'Invalid mess_command: {error}') from error
    if not command or not command[0]:
        raise MESSExecutionError('mess_command must specify a nonempty executable.')
    return command


def _run(command, *, timeout=None):
    return subprocess.run(command, capture_output=True, text=True, check=False,
                          timeout=timeout)


def _submit(queue, script):
    try:
        result = _run([constants.qsubmit[queue], str(script)],
                      timeout=_QUERY_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired as error:
        raise _SubmissionUncertain(
            'MESS scheduler submission timed out; inspect the queue before '
            'retrying because the submission outcome is unknown.') from error
    except OSError as error:
        raise MESSExecutionError(
            f'MESS scheduler submission could not start: {error}') from error
    if result.returncode:
        raise MESSExecutionError(f'MESS submission failed: {result.stderr.strip()}')
    if queue == 'slurm':
        match = re.search(r'\bSubmitted batch job (\d+)\b', result.stdout)
    else:
        match = re.fullmatch(r'\s*(\d+(?:\.[\w.-]+)?)\s*', result.stdout)
    if not match:
        raise _SubmissionUncertain(
            f'Cannot read MESS job ID from {result.stdout!r}; inspect the '
            'queue before retrying because the submission outcome is '
            'unknown.')
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


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


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
    return {
        'output': str(output), 'output_sha256': _sha256(output),
        'output_bytes': current_output[1],
    }


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
    solver = mess_command(writer.par)
    for index in indices:
        extension = '.sh' if queue == 'local' else constants.qext[queue]
        script = Path('me') / f'run_mess_{index}{extension}'
        if queue == 'local':
            command = shlex.join(solver + [f'mess_{index}.inp'])
            script.write_text(f'#!/bin/sh\ncd me || exit 1\nexec {command}\n')
        else:
            writer.write_submitscript(str(script), index)
        scripts.append(script)
    command = 'sh' if queue == 'local' else constants.qsubmit[queue]
    batch = Path('batch_me.sub')
    batch.write_text(''.join(f'{command} {shlex.quote(str(script))}\n' for script in scripts))
    batch.chmod(0o700)
    if not writer.par['run_me']:
        return 0

    automatic = bool(writer.par.get('anl_handoff_manifest'))
    destination = Path('me/mess_execution.json')
    expected_jobs = {f'mess_{index}' for index in indices}
    input_hashes = {
        index: _sha256(Path('me') / f'mess_{index}.inp')
        for index in indices
    }
    if destination.is_file():
        try:
            existing = json.loads(destination.read_text())
            jobs = existing['jobs']
            reusable = (
                existing.get('schema') == 1
                and existing.get('status') == 'complete'
                and existing.get('queue') == queue
                and existing.get('command') == solver
                and set(jobs) == expected_jobs)
            for name, details in jobs.items():
                output = Path(details['output'])
                input_path = Path(details['input'])
                reusable &= (
                    output.is_file()
                    and _sha256(output) == details['output_sha256']
                    and input_path.is_file()
                    and _sha256(input_path) == details['input_sha256'])
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            reusable = False
        if reusable:
            logger.info('Reusing verified completed MESS execution record.')
            return 0

    progress_path = Path('me/mess_execution.running.json')
    active = {}
    completed = {}
    submitting = None
    if progress_path.is_file():
        try:
            progress = json.loads(progress_path.read_text())
            identity_matches = (
                progress.get('schema') == 1
                and progress.get('queue') == queue
                and progress.get('command') == solver
                and progress.get('indices') == indices
                and progress.get('input_sha256') == input_hashes)
            if not identity_matches:
                if progress.get('status') != 'failed':
                    raise ValueError(
                        'saved running MESS execution identity changed')
                os.replace(
                    progress_path,
                    progress_path.with_name(
                        progress_path.name + '.previous-failure'))
                progress = None
            if progress is None:
                pass
            elif progress.get('status') == 'failed':
                uncertain = progress.get('submitting')
                if uncertain is not None:
                    raise _SubmissionUncertain(
                        'MESS submission outcome is uncertain for '
                        f'{uncertain.get("index", "unknown")}; inspect the '
                        'scheduler before clearing the saved progress record.')
                raise MESSExecutionError(
                    'Previous MESS execution failed: ' +
                    ' | '.join(progress.get('failures', ())))
            else:
                submitting = progress.get('submitting')
                if submitting is not None:
                    raise MESSExecutionError(
                        'MESS submission outcome is uncertain for '
                        f'{submitting.get("index", "unknown")}; inspect the '
                        'scheduler before clearing the saved progress record.')
                completed = dict(progress.get('completed', {}))
                if not set(completed).issubset(indices):
                    raise ValueError(
                        'saved MESS completed job set is invalid')
                restored_indices = []
                for pid, value in progress.get('active', {}).items():
                    previous = value.get('previous_output')
                    if value['index'] not in indices:
                        raise ValueError(
                            'saved active MESS job is not requested')
                    restored_indices.append(value['index'])
                    active[str(pid)] = (
                        value['index'], tuple(previous) if previous else None)
                if (len(restored_indices) != len(set(restored_indices))
                        or set(restored_indices) & set(completed)):
                    raise ValueError(
                        'saved MESS jobs contain duplicate execution state')
        except MESSExecutionError:
            raise
        except (OSError, AttributeError, KeyError, TypeError, ValueError,
                json.JSONDecodeError) as error:
            raise MESSExecutionError(
                f'Invalid saved MESS execution state at {progress_path}: '
                f'{error}') from error

    ended = {}
    query_failures = {}
    failures = []
    limit = max(1, int(writer.par['uq_max_runs']))

    def save_progress(status='running'):
        _atomic_json(progress_path, {
            'schema': 1, 'status': status, 'queue': queue,
            'command': solver, 'indices': indices,
            'input_sha256': input_hashes,
            'active': {
                str(pid): {
                    'index': index,
                    'previous_output': list(previous)
                    if previous is not None else None,
                }
                for pid, (index, previous) in active.items()
            },
            'completed': completed,
            'submitting': submitting,
            'failures': failures,
        })

    save_progress()

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
                    if automatic and _scheduler_unavailable(error):
                        if count == _MAX_QUERY_FAILURES:
                            logger.warning(
                                '%s Automatic MESS execution retains job %s '
                                'and will keep waiting for scheduler recovery.',
                                error, pid)
                        continue
                    raise MESSExecutionError(
                        f'{error} Monitoring stopped after {count} consecutive query failures; '
                        f'submitted jobs may still be active: {", ".join(active)}.') from error
                query_failures.pop(pid, None)
                if running:
                    continue
                ended[pid] = time.monotonic()
            try:
                completed[index] = _verify(index, previous)
            except _ResultPending as error:
                if time.monotonic() - ended[pid] < _RESULT_GRACE_SECONDS:
                    continue
                failures.append(f'{error} Result files did not become available within '
                                f'{_RESULT_GRACE_SECONDS:g} seconds after queue completion.')
            except MESSExecutionError as error:
                failures.append(str(error))
            del active[pid]
            ended.pop(pid, None)
            save_progress()

    scheduled = {index for index, _ in active.values()}
    for index, script in zip(indices, scripts):
        if index in completed or index in scheduled:
            continue
        output = Path('me') / f'mess_{index}.out'
        previous = _stamp(output)
        if queue == 'local':
            with open(f'me/mess_{index}.stdout', 'w') as stdout, \
                    open(f'me/mess_{index}.err', 'w') as stderr:
                try:
                    result = subprocess.run(solver + [f'mess_{index}.inp'], cwd='me',
                                            stdout=stdout, stderr=stderr, check=False)
                except OSError as error:
                    failures.append(
                        f'Cannot start MESS calculation {index} with {solver[0]!r}: {error}. '
                        'Check mess_command and the executable in the local environment.')
                    break
            try:
                completed[index] = _verify(
                    index, previous, result.returncode)
            except MESSExecutionError as error:
                failures.append(str(error))
                break
            save_progress()
        else:
            while len(active) >= limit:
                poll()
            if failures:
                break
            (Path('me') / f'mess_{index}.exitcode').unlink(missing_ok=True)
            while True:
                submitting = {'index': index, 'script': str(script)}
                save_progress()
                try:
                    pid = _submit(queue, script)
                    break
                except MESSExecutionError as error:
                    if isinstance(error, _SubmissionUncertain):
                        failures.append(str(error))
                        break
                    submitting = None
                    save_progress()
                    if automatic and _scheduler_unavailable(error):
                        logger.warning(
                            '%s No MESS job was accepted; retrying after '
                            'scheduler recovery.', error)
                        time.sleep(5)
                        continue
                    failures.append(str(error))
                    break
            if failures:
                break
            submitting = None
            active[pid] = (index, previous)
            scheduled.add(index)
            save_progress()
    while active:
        poll()
    if failures:
        save_progress(status='failed')
        raise MESSExecutionError('\n'.join(failures))
    if set(completed) != set(indices):
        raise MESSExecutionError(
            'MESS execution ended without every requested job result.')
    for index, details in completed.items():
        input_path = Path('me') / f'mess_{index}.inp'
        output_path = Path(details['output'])
        if (not output_path.is_file()
                or _sha256(output_path) != details['output_sha256']):
            raise MESSExecutionError(
                f'MESS calculation {index}: saved output changed before '
                'the completion record was published.')
        details['input'] = str(input_path)
        details['input_sha256'] = input_hashes[index]
    record = {
        'schema': 1, 'status': 'complete', 'queue': queue,
        'command': solver,
        'jobs': {
            f'mess_{index}': details
            for index, details in sorted(completed.items())
        },
    }
    networks = Path('me/mess_networks.json')
    if networks.is_file():
        record['network_manifest'] = str(networks)
        record['network_manifest_sha256'] = _sha256(networks)
    _atomic_json(destination, record)
    progress_path.unlink(missing_ok=True)
    return 0
