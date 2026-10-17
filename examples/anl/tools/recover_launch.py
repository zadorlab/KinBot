"""Finish an interrupted handoff into a fresh, distinct production directory.

Never resume scientific results from a retired run. A repeated invocation
leaves live new drivers alone and starts only cases that never started.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import socket

from kinbot import run_cleanup
from kinbot.execution import atomic_json
from kinbot.run_format import _existing_results
import run_examples as handoff


def progress(root, phase, **details):
    record = {'time_utc': datetime.now(timezone.utc).isoformat(),
              'host': socket.gethostname(), 'pid': os.getpid(),
              'phase': phase, **details}
    atomic_json(root/'restart_status.json', record)
    print(f'{record["time_utc"]} {phase}', flush=True)


def argv_for(record):
    """Read argument boundaries rather than matching a shell command substring."""
    try:
        return (Path('/proc')/str(record['pid'])/'cmdline').read_bytes().decode().rstrip('\0').split('\0')
    except (OSError, UnicodeError):
        return []


def old_handoffs(previous):
    """Only this package's exact launcher invocation for the retired root."""
    selected = []
    for process in run_cleanup._processes(Path('/')):
        argv = argv_for(process)
        if len(argv) < 3 or not Path(argv[0]).name.startswith('python'):
            continue
        # The Python invocation can contain -u before the script filename.
        script = next((i for i, arg in enumerate(argv)
                       if arg in (str(handoff.PACKAGE/'run_examples.py'),
                                  str(handoff.PACKAGE/'recover_launch.py'))), None)
        if script is None:
            continue
        args = argv[script + 1:]
        if argv[script].endswith('/run_examples.py'):
            if len(args) < 2 or args[0] not in ('prepare', 'preflight', 'launch', 'resume'):
                continue
            target = args[1]
        else:
            if not args:
                continue
            target = args[0]
        if Path(target).expanduser().resolve() == previous:
            selected.append(process)
    return selected


def case_states(root):
    states = []
    for path in handoff.verify(root)['inputs']:
        directory = Path(path).parent
        live = run_cleanup.inspect(directory)  # Query failure must propagate.
        launch_file, status_file = directory/'launch.json', directory/'pes_status.json'
        launch = json.loads(launch_file.read_text()) if launch_file.exists() else {}
        status = json.loads(status_file.read_text()) if status_file.exists() else {}
        # /proc only sees this host. Never decide that a remote driver died.
        host = status.get('host', launch.get('host'))
        if host and host.split('.')[0] != socket.gethostname().split('.')[0]:
            raise RuntimeError(f'{directory}: check/retry on driver host {host}; no jobs changed.')
        drivers = []
        for process in live['processes']:
            argv = argv_for(process)
            if any(argv[i:i+3] == ['-m', 'kinbot.pes', path] for i in range(len(argv)-2)):
                drivers.append(process['pid'])
        if len(drivers) > 1:
            raise RuntimeError(f'{directory}: multiple PES drivers {drivers}; no new driver started.')
        if drivers:
            state = 'running'
        elif status.get('status') == 'complete' and not (live['jobs'] or live['processes']):
            state = 'complete'
        elif (live['jobs'] or live['processes'] or launch or status
              or _existing_results(directory) is not None):
            raise RuntimeError(f'{directory}: previously started but no live PES driver; '
                               'inspect the parent monitor. No duplicate or blind scientific restart was launched.')
        else:
            state = 'unstarted'
        states.append({'input': path, 'state': state, 'pids': drivers})
    return states


def recover(root, previous):
    root = Path(root).expanduser().resolve()
    previous = [Path(path).expanduser().resolve() for path in previous]
    for path in previous:
        if path == root or path in root.parents or root in path.parents:
            raise ValueError('Fresh and retired roots must be separate directories.')
    root.mkdir(parents=True, exist_ok=True)
    with (root/'.handoff.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(f'Another handoff is still active for {root}; no duplicate started. '
                  'See restart_status.json and the existing restart log.', flush=True)
            return {'status': 'handoff_already_active'}
        try:
            progress(root, 'Preparing/verifying single-JSON inputs')
            handoff.prepare(root)
            states = case_states(root)
            progress(root, 'Checking production backends (may take a few minutes)', cases=states)
            # The same check is safe after an interrupted preflight or retirement.
            checks = handoff.preflight(root)
            report = {'root': str(root), 'preflight': checks, 'retired': [], 'cases_before': states}
            for path in previous:
                try:
                    path.stat()
                except FileNotFoundError:
                    print(f'Already absent; leaving existing archives intact: {path}', flush=True)
                    report['retired'].append({'previous': str(path), 'status': 'already_absent'})
                    continue
                progress(root, f'Stopping and archiving {path.name}')
                # Includes an old prepare/preflight/launch helper whose cwd is
                # the repository, plus its children, by PID/starttime identity.
                identities = old_handoffs(path)
                report['retired'].append(handoff.retire_production(
                    root, path, known=identities, allow_partial=True))
                atomic_json(root/'restart_report.json', report)
            report['archived_repository_layout'] = handoff.archive_legacy_example_layout(root)
            states = case_states(root)  # Recheck before starting anything.
            pending = [item['input'] for item in states if item['state'] == 'unstarted']
            progress(root, f'Starting {len(pending)} unstarted PES driver(s)', cases=states)
            report['launched'] = handoff.start_inputs(pending)
            report['retained'] = [item for item in states if item['state'] != 'unstarted']
            atomic_json(root/'restart_report.json', report)
            progress(root, 'Handoff finished; check both PES with the monitor',
                     launched=report['launched'], retained=report['retained'])
            print(json.dumps({'run': str(root), 'launched': report['launched'],
                              'retained': report['retained']}, indent=2), flush=True)
            return report
        except Exception as error:
            progress(root, 'Handoff stopped; rerun after the reported issue is resolved',
                     error=f'{type(error).__name__}: {error}')
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--retire', type=Path, action='append', default=[])
    args = parser.parse_args()
    recover(args.root, args.retire)
