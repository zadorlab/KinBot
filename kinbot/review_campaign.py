"""Run an explicit inventory of independent orbital and relaxation campaigns.

Scientific recipes and molecular identities belong in supplied manifests.
Workers retain their own locks, native dependencies and allocation histories.
This supervisor never accepts a model or replaces a production controller.
"""
import argparse
from contextlib import ExitStack
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from kinbot import requested_orbital_recon, run_policy
from kinbot.execution import atomic_json


MODULES = {'orbitals': 'kinbot.requested_orbital_recon',
           'relaxations': 'kinbot.vrc_relaxation_recon'}


def inventory(root, records):
    root = Path(root).resolve(strict=True)
    result, names, directories = [], set(), set()
    for row in records:
        name, kind = row['name'], row['kind']
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', name):
            raise ValueError('Campaign names must be simple, unique log labels.')
        directory = (root/row['directory']).resolve(strict=True)
        if (kind not in MODULES or name in names or directory in directories
                or not directory.is_relative_to(root)):
            raise ValueError('Invalid or duplicate review campaign.')
        names.add(name); directories.add(directory)
        result.append(dict(name=name, kind=kind, directory=str(directory)))
    if not result:
        raise ValueError('No review campaigns were specified.')
    return result


def supervise(root, records, *, code=None, allow_nodes=(), queue_all=False):
    """Caller holds the supervisor lock; each worker independently adopts jobs."""
    root = Path(root).resolve(strict=True)
    records = inventory(root, records)
    allowed = run_policy._nodes(allow_nodes)
    if allowed and not queue_all:
        raise ValueError('Node overrides require explicit --queue-all review scheduling.')
    children = []
    with ExitStack() as logs:
        for row in records:
            log = logs.enter_context((root/('queue_all_'+row['name']+'.log')).open('a'))
            command = [sys.executable, '-B', '-u', '-m', MODULES[row['kind']],
                       'run', row['directory']]
            if queue_all:
                command.append('--queue-all')
            if row['kind'] == 'orbitals':
                command.append('--report-failures')
            for node in allowed:
                command += ['--allow-node', node]
            env = dict(os.environ)
            if code is not None:
                env['PYTHONPATH'] = str(code)
            child = subprocess.Popen(command, cwd=root, stdin=subprocess.DEVNULL,
                stdout=log, stderr=subprocess.STDOUT, env=env)
            children.append((row['name'], child))
            print(f'{row["name"]}: owner {child.pid}; log {log.name}', flush=True)
        marker = root/'queue_all_refresh.json'
        journal = json.loads(marker.read_text()) if marker.exists() else {}
        journal.update(phase='running', campaigns=records,
                       owners={name:p.pid for name,p in children})
        atomic_json(marker, journal)
        while any(p.poll() is None for _,p in children):
            print('Review owners (null = active): '+json.dumps(
                {name:p.poll() for name,p in children}), flush=True)
            time.sleep(20)
        result = {name:p.returncode for name,p in children}
    for row in records:
        if row['kind'] == 'orbitals':
            folder = Path(row['directory'])
            (folder/'orbital_evidence.txt').write_text(requested_orbital_recon.evidence(folder)+'\n')
    atomic_json(root/'queue_all_finished.json', dict(
        owner_exit_codes=result, production_model_accepted=False))
    print('Diagnostic owners finished; native results still require review: '+json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path, help='JSON with schema=1 and a campaigns list.')
    parser.add_argument('--queue-all', action='store_true')
    parser.add_argument('--allow-node', action='append', default=[])
    args = parser.parse_args()
    manifest = args.manifest.resolve(strict=True)
    config = json.loads(manifest.read_text())
    if config.get('schema') != 1:
        raise ValueError('Expected review campaign schema 1.')
    root = manifest.parent
    records = inventory(root, config['campaigns'])
    # A repeat adopts completed/live native jobs after the previous owner exits.
    # Never remove a lock file or kill an owner merely because a PID looks old.
    with (root/'.queue_all.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = supervise(root, records, allow_nodes=args.allow_node, queue_all=args.queue_all)
    raise SystemExit(int(any(result.values())))


if __name__ == '__main__':
    main()
