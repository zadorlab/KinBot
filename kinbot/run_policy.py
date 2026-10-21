"""Live, shared Slurm admission for cooperating KinBot calculation trees.

A policy in an ancestor directory applies to every descendant controller.
Only scheduling is changed: native task files and running jobs stay intact.
"""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import threading
import time
import uuid

from kinbot.execution import atomic_json, slurm_job_record, SLURM_TERMINAL_STATES

NAME = '.kinbot_resource_policy.json'
_LOCK = threading.RLock()


def find(directory=None):
    path = Path(directory or Path.cwd()).resolve()
    found = None
    for parent in (path, *path.parents):
        candidate = parent / NAME
        if candidate.is_file():
            if found is not None:
                raise ValueError(f'Nested shared policies would split admission: {found} and {candidate}')
            found = candidate
    return found


def _nodes(values):
    if not isinstance(values, (list, tuple, set)) or any(
            not isinstance(v, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', v)
            for v in values):
        raise ValueError('Node lists require explicit individual node names.')
    return sorted(set(values))


def _read(path):
    value = json.loads(path.read_text())
    if value.get('schema') != 1 or type(value.get('max_jobs')) is not int or value['max_jobs'] < 0:
        raise ValueError(f'{path}: invalid shared scheduling policy.')
    _nodes(value['excluded_nodes'])
    if 'minimum_walltime' in value:
        from kinbot.run_walltime import validate_minimum
        validate_minimum(value['minimum_walltime'])
    return value


def read(directory=None):
    path = find(directory)
    return _read(path) if path else None


def limit(default, directory=None):
    policy = read(directory)
    return policy['max_jobs'] if policy is not None else default


@contextmanager
def _locked(path):
    with _LOCK, (path.parent / '.kinbot_resource_policy.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def _initial_exclusions(root):
    """Preserve input and learned exclusions when a pool is first enabled."""
    nodes = set()
    inputs = [root/'input.json']
    prepared = root/'prepared.json'
    if prepared.is_file():
        record = json.loads(prepared.read_text())
        inputs.extend(Path(p).resolve() for p in record.get('inputs', []))
    for path in inputs:
        if path.is_file():
            nodes.update(_nodes(json.loads(path.read_text()).get('slurm_exclude_nodes', [])))
    scopes = {root, *(path.parent for path in inputs)}
    for scope in scopes:
        if not scope.is_relative_to(root):
            raise ValueError('Prepared input is outside the shared pool.')
        for prefix in ('', '[0-9]*/'):
            for pattern in ('.kinbot_slurm_nodes.json',
                            'anl_automatic/states/*/composite/state.json',
                            'vrctst/*/*/state.json'):
                for path in scope.glob(prefix + pattern):
                    state = json.loads(path.read_text())
                    nodes.update(_nodes(state.get('excluded_nodes', [])))
                    nodes.update(_nodes(state.get('scratch_excluded_nodes', [])))
    return sorted(nodes)


def update(root, *, max_jobs=None, exclude=(), allow=(), minimum_walltime=None):
    root = Path(root).resolve(strict=True)
    path = root/NAME
    ancestor = find(root.parent)
    if ancestor:
        raise ValueError(f'An ancestor pool already manages this tree: {ancestor}')
    if max_jobs is not None and (type(max_jobs) is not int or max_jobs < 0):
        raise ValueError('max_jobs must be a nonnegative integer; zero pauses new admissions.')
    exclude, allow = _nodes(exclude), _nodes(allow)
    if minimum_walltime not in (None, ''):
        from kinbot.run_walltime import validate_minimum
        validate_minimum(minimum_walltime)
    if set(exclude) & set(allow):
        raise ValueError('A node cannot be excluded and allowed in the same update.')
    with _locked(path):
        previous = _read(path) if path.exists() else None
        if previous is None:
            if max_jobs is None:
                raise ValueError('Specify --max-jobs when enabling a shared pool.')
            # Nested policies would partition the supposedly shared pool.
            if any(p != path for p in root.rglob(NAME)):
                raise ValueError('Child policies already exist; use their common pool instead.')
        value = dict(previous or {'schema': 1, 'max_jobs': max_jobs,
                                 'excluded_nodes': _initial_exclusions(root), 'revision': 0})
        if max_jobs is not None:
            value['max_jobs'] = max_jobs
        if minimum_walltime == '':
            value.pop('minimum_walltime', None)
        elif minimum_walltime is not None:
            value['minimum_walltime'] = minimum_walltime
        value['excluded_nodes'] = sorted((set(value['excluded_nodes']) | set(exclude)) - set(allow))
        value.update(revision=value['revision'] + 1,
                     updated_at=datetime.now(timezone.utc).isoformat())
        archive = root/'.kinbot_policy_history'/f'{value["revision"]:06d}.json'
        atomic_json(archive, {'before': previous, 'after': value})
        atomic_json(path, value)
    return value


def validate_options(options):
    if options == {}:
        return
    if not isinstance(options, dict) or set(options) - {'max_jobs', 'excluded_nodes', 'minimum_walltime'}:
        raise ValueError('queue_live_policy accepts max_jobs, excluded_nodes and minimum_walltime.')
    if type(options.get('max_jobs')) is not int or options['max_jobs'] < 0:
        raise ValueError('queue_live_policy.max_jobs must be a nonnegative integer.')
    _nodes(options.get('excluded_nodes', []))
    if 'minimum_walltime' in options:
        from kinbot.run_walltime import validate_minimum
        validate_minimum(options['minimum_walltime'])


def initialize(directory, options):
    """Optional single-input defaults; later live edits take precedence."""
    validate_options(options)
    if not options:
        return
    if not find(directory):
        update(directory, max_jobs=options['max_jobs'], exclude=options.get('excluded_nodes', []),
               minimum_walltime=options.get('minimum_walltime'))


def exclude_failed_node(directory, node):
    path = find(directory)
    if path and node not in _read(path)['excluded_nodes']:
        update(path.parent, exclude=[node])


def _queue(root):
    from kinbot.qc_scratch import account_identity
    from kinbot.anl.dispatch import SchedulerUnavailableError
    try:
        result = subprocess.run(['squeue', '-h', '-u', account_identity()['name'],
                                 '-o', '%i|%D|%Z'], capture_output=True, text=True,
                                check=False, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise SchedulerUnavailableError(f'Shared pool cannot query Slurm: {error}') from error
    if result.returncode:
        raise SchedulerUnavailableError(f'Shared pool cannot query Slurm: {result.stderr.strip()}')
    jobs = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split('|', 2)
        if len(fields) != 3:
            raise RuntimeError(f'Unrecognized shared-pool squeue row: {line!r}')
        ident, nodes, directory = [f.strip() for f in fields]
        if not directory or directory in ('(null)', 'N/A'):
            raise RuntimeError('Slurm omitted a work directory; pool occupancy is unknown.')
        if Path(directory).resolve().is_relative_to(root):
            if not ident.isdecimal() or not nodes.isdecimal() or int(nodes) < 1:
                raise RuntimeError('Shared pools currently require individual single-node jobs, not arrays.')
            jobs[ident] = int(nodes)
    return jobs


def _reconcile(root, jobs, claims):
    for token, claim in list(claims.items()):
        ident = claim.get('job_id')
        if ident in jobs:
            del claims[token]  # Visible jobs are now counted by Slurm itself.
        elif ident:
            record = slurm_job_record(ident)
            if record.get('JobId') == ident and record.get('JobState', '').split('+')[0] in SLURM_TERMINAL_STATES:
                del claims[token]
    return sum(jobs.values()) + len(claims)


def _submit_once(command, directory, path, kwargs):
    with _locked(path):
        return _admit_locked(command, directory, path, kwargs)


def submit(command, *, cwd=None, wait=False, **kwargs):
    """Serialize count + sbatch across threads/processes, retaining ambiguity.

    Return None when the pool is full and wait=False. No scheduler submission
    has occurred in that case. Submission timeouts reserve a slot until review;
    they are never automatically resubmitted.
    """
    directory = Path(cwd or Path.cwd()).resolve()
    path = find(directory)
    if path is None:
        return subprocess.run(command, cwd=cwd, **kwargs)
    if Path(command[0]).name != 'sbatch':
        raise ValueError('Shared resource policies currently support Slurm submissions only.')
    # All KinBot callers submit a script, optionally with --parsable. Refuse
    # competing command-line resource options rather than silently overriding
    # either a caller's requirements or the shared pool's policy.
    if (len(command) < 2 or str(command[-1]).startswith('-')
            or any(arg != '--parsable' for arg in command[1:-1])):
        raise ValueError('Shared admission requires sbatch [--parsable] SCRIPT.')
    from kinbot.anl.dispatch import SchedulerUnavailableError
    warned = False
    while True:
        try:
            result = _submit_once(command, directory, path, kwargs)
        except SchedulerUnavailableError as error:
            if not wait:
                raise
            if not warned:
                import sys
                print(f'Shared pool waiting for scheduler: {error}', file=sys.stderr, flush=True)
                warned = True
        else:
            warned = False
            if result is not None or not wait:
                return result
        time.sleep(20)


def _admit_locked(command, directory, path, kwargs):
    policy = _read(path)
    ledger = path.parent/'.kinbot_policy_admissions.json'
    claims = json.loads(ledger.read_text()) if ledger.exists() else {}
    occupied = _reconcile(path.parent, _queue(path.parent), claims)
    atomic_json(ledger, claims)
    if any(c['directory'] == str(directory) and c['command'] == list(command)
           for c in claims.values()):
        raise RuntimeError('A previous pool submission is unconfirmed; reconcile it before resubmitting.')
    if occupied < policy['max_jobs']:
        override = {}
        if policy.get('minimum_walltime'):
            from kinbot.run_walltime import submission_override
            from kinbot.anl.dispatch import SchedulerUnavailableError
            try:
                override = submission_override(directory/command[-1], policy['minimum_walltime'],
                                               env=kwargs.get('env'))
            except (subprocess.SubprocessError, OSError, RuntimeError) as error:
                # No sbatch has occurred. A failed partition lookup must not
                # become an uncertain submission or consume a pool slot.
                raise SchedulerUnavailableError(f'Cannot resolve walltime policy: {error}') from error
        token = uuid.uuid4().hex
        claim = {'directory': str(directory), 'command': list(command),
                 'created_at': datetime.now(timezone.utc).isoformat(),
                 'policy_revision': policy['revision'], 'excluded_nodes': policy['excluded_nodes'],
                 'status': 'submitting'}
        if override:
            claim['walltime_override'] = override
        claims[token] = claim
        atomic_json(ledger, claims)  # Durable before the first sbatch byte.
        # Command-line policy overrides stale script/environment exclusions.
        actual = [*command[:-1], '--nodes=1', '--exclusive',
                  '--exclude=' + ','.join(policy['excluded_nodes']),
                  *(['--time='+override['effective']] if override else []), command[-1]]
        try:
            response = subprocess.run(actual, cwd=directory, **kwargs)
        except OSError:
            del claims[token]  # Failed exec: sbatch never started.
            atomic_json(ledger, claims)
            raise
        # Timeout/BaseException intentionally leaves an uncertain claim.
        if response.returncode:
            del claims[token]
            atomic_json(ledger, claims)
            return response
        stdout = response.stdout.decode() if isinstance(response.stdout, bytes) else response.stdout
        match = re.fullmatch(r'(?:Submitted batch job\s+)?(\d+)(?:;[^\s]+)?\s*', stdout or '')
        if not match:
            raise RuntimeError('Shared pool received an unrecognized sbatch response; submission is uncertain.')
        claim.update(job_id=match.group(1), status='submitted')
        atomic_json(ledger, claims)
        atomic_json(path.parent/'.kinbot_policy_submissions'/f'{claim["job_id"]}.json', claim)
        return response
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    change = sub.add_parser('set')
    change.add_argument('root', type=Path)
    change.add_argument('--max-jobs', type=int)
    change.add_argument('--exclude', action='append', default=[], help='Comma-separated explicit node names to exclude.')
    change.add_argument('--allow', action='append', default=[], help='Remove explicit/learned exclusions for these nodes.')
    walltime = change.add_mutually_exclusive_group()
    walltime.add_argument('--minimum-walltime', help='Raise future requests up to this time, capped by partition MaxTime.')
    walltime.add_argument('--clear-minimum-walltime', action='store_true')
    show = sub.add_parser('show'); show.add_argument('root', type=Path)
    launch = sub.add_parser('submit'); launch.add_argument('script')
    args = parser.parse_args(argv)
    if args.action == 'set':
        result = update(args.root, max_jobs=args.max_jobs,
                        minimum_walltime='' if args.clear_minimum_walltime else args.minimum_walltime,
                        exclude=[n for v in args.exclude for n in v.split(',')],
                        allow=[n for v in args.allow for n in v.split(',')])
        print(json.dumps(result, indent=2))
        print('Policy saved. Running and pending jobs were not modified.')
    elif args.action == 'show':
        path = find(args.root)
        if not path:
            raise SystemExit('No shared policy found.')
        with _locked(path):
            ledger = path.parent/'.kinbot_policy_admissions.json'
            claims = json.loads(ledger.read_text()) if ledger.exists() else {}
            jobs = _queue(path.parent)
            occupied = _reconcile(path.parent, jobs, claims)  # Read-only in this command.
            print(json.dumps({'root': str(path.parent), 'policy': _read(path),
                              'live_nodes': sum(jobs.values()), 'reserved_slots': occupied,
                              'unreconciled_admissions': claims}, indent=2))
    else:
        result = submit(['sbatch', args.script], wait=True, capture_output=True,
                        text=True, check=False, timeout=30)
        print(result.stdout, end='')
        if result.stderr:
            import sys
            print(result.stderr, end='', file=sys.stderr)
        raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
