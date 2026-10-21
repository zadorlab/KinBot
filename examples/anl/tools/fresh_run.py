"""Prepare/review a fresh two-PES run, then explicitly retire its predecessor.

Repeat the same invocation after an interruption. A completed launch is never
repeated; an uncertain process launch requires inspection instead of guessing.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import socket
import subprocess

from kinbot import run_cache_cleanup, run_cleanup, run_policy
from kinbot.execution import atomic_json, file_sha256
from kinbot.run_format import ensure_current_run
import run_examples as handoff


def node_plan(partition, reserve, preferred=()):
    """Use an expanded single-partition snapshot; unknown states are excluded."""
    if type(reserve) is not int or reserve < 0:
        raise ValueError('reserve_healthy_nodes must be a nonnegative integer.')
    q = subprocess.run(['sinfo', '-N', '-h', '-p', partition, '-o', '%N|%T'],
                       capture_output=True, text=True, check=True, timeout=20)
    states = {}
    for line in q.stdout.splitlines():
        name, state = (v.strip() for v in line.split('|'))
        run_policy._nodes([name])
        states.setdefault(name, set()).add(state.lower())
    healthy = {n for n, s in states.items() if s <= {'idle', 'allocated', 'mixed', 'completing'}}
    if len(healthy) <= reserve:
        raise RuntimeError('Too few healthy partition nodes after the requested reserve.')
    selected = sorted(healthy, key=lambda n: (n not in preferred, states[n] != {'idle'}, n))[:reserve]
    return dict(partition=partition, states={n: sorted(s) for n, s in states.items()},
                reserved_healthy_nodes=selected, max_jobs=len(healthy)-reserve,
                excluded_nodes=sorted((set(states)-healthy) | set(selected)))


def pin_code(root):
    import kinbot
    source = Path(kinbot.__file__).parent
    hashes = {'kinbot/'+str(p.relative_to(source)): file_sha256(p)
              for p in source.rglob('*') if p.is_file()
              and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    digest = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    pin = root/('.fresh_code.'+digest[:16])
    if not pin.exists():
        import tempfile
        temp = Path(tempfile.mkdtemp(prefix='.fresh_code.preparing.', dir=root))
        shutil.copytree(source, temp/'kinbot', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        atomic_json(temp/'manifest.json', hashes)
        # Check the copied bytes before publishing the snapshot.
        verify_pin(temp)
        temp.rename(pin)
    verify_pin(pin)
    return pin


def verify_pin(pin):
    for name, digest in json.loads((pin/'manifest.json').read_text()).items():
        if file_sha256(pin/name) != digest:
            raise RuntimeError('Pinned source changed; refusing a mixed-version launch.')


def retire(previous, root, journal):
    """A journaled parent move follows confirmed termination of both trees."""
    archive = Path(journal['archive'])
    if archive.exists():
        if previous.exists():
            raise RuntimeError('Both previous and archive exist; inspect before proceeding.')
        if not (archive/'fresh_retirement.json').exists() or json.loads(
                (archive/'fresh_retirement.json').read_text()).get('replacement') != str(root):
            raise RuntimeError('Archive does not match this fresh-run handoff.')
        return
    if not previous.is_dir():
        raise RuntimeError('Previous run is missing and its archive is unconfirmed.')
    directories = [previous/c for c in handoff.CASES]
    for d in directories:
        ensure_current_run(d)
        for name in ('launch.json', 'pes_status.json'):
            p = d/name
            if p.exists():
                host = json.loads(p.read_text()).get('host')
                if host and host.split('.')[0] != socket.gethostname().split('.')[0]:
                    raise RuntimeError(f'Retire on controller host {host}; not on this host.')
    # All inventory queries must succeed before cancelling either tree.
    inventories = [run_cleanup.inspect(d) for d in directories]
    atomic_json(previous/'fresh_retirement.json', dict(replacement=str(root),
                archive=str(archive), inventories=inventories))
    for d in directories:
        run_cleanup.stop(d, timeout=300)
    for d in directories:
        live = run_cleanup.inspect(d)
        if live['jobs'] or live['processes']:
            raise RuntimeError(f'{d}: work remains; not archived.')
    previous.rename(archive)


def run(root, previous, config, *, apply=False, skip_cleanup=False):
    run_cache_cleanup.validate(config.get('cleanup', {}))
    root = Path(root).expanduser().resolve()
    previous = Path(previous).expanduser().resolve()
    if root == previous or root in previous.parents or previous in root.parents:
        raise ValueError('Fresh and previous runs must be separate trees.')
    root.mkdir(parents=True, exist_ok=True)
    with (root/'.handoff.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = root/'fresh_run.json'
        if path.exists():
            journal = json.loads(path.read_text())
            if journal['previous'] != str(previous) or journal['config'] != config:
                raise ValueError('Existing fresh-run journal has different arguments/configuration.')
            if journal['phase'] == 'launched':
                return journal  # Never turn a repeated handoff into another launch.
        else:
            unexpected = [p for p in root.iterdir() if p.name != '.handoff.lock']
            if unexpected:
                raise ValueError('Fresh destination contains unrecognized files.')
            node = node_plan(config['partition'], config['reserve_healthy_nodes'],
                             config.get('preferred_reserve_nodes', []))
            policy = dict(config['policy'], max_jobs=node['max_jobs'],
                          excluded_nodes=node['excluded_nodes'])
            run_policy.validate_options(policy)
            journal = dict(schema=1, phase='preparing', previous=str(previous), config=config,
                nodes=node, policy=policy, launches={}, archive=str(previous.with_name(
                    previous.name+'.stopped.'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))))
            atomic_json(path, journal)
        if journal['phase'] == 'preparing':
            overrides = dict(config.get('input_overrides', {}), queue_name=config['partition'],
                             slurm_exclude_nodes=journal['nodes']['excluded_nodes'],
                             queue_live_policy=journal['policy'])
            handoff.prepare(root, overrides=overrides)
            run_policy.initialize(root, journal['policy'])
            from kinbot import live_control
            live_control.initialize(root, overrides.get('queue_live_control', {}))
            journal.update(phase='prepared', code=str(pin_code(root)))
            atomic_json(path, journal)
        verify_pin(Path(journal['code']))
        if journal['phase'] in ('prepared', 'reviewed'):
            handoff.assert_not_started(root)
            print('Checking replacement inputs and backend environment before retirement...', flush=True)
            handoff.preflight(root)
            journal.update(phase='reviewed', code=str(pin_code(root)))
            atomic_json(path, journal)
        if not apply:
            return journal
        if journal['phase'] == 'reviewed':
            print('Stopping only the previous two PES trees; retaining their files in the archive.', flush=True)
            retire(previous, root, journal)
            journal['phase'] = 'retired'
            atomic_json(path, journal)
        if journal['phase'] == 'retired':
            handoff.verify(root)  # Reviewed WFU seeds must survive cache removal.
            if skip_cleanup:
                # This is an invocation-level opt-out, not a change to the
                # immutable prepared chemistry/configuration. Keep any partial
                # cleanup journal and remaining cache contents as evidence.
                print('Skipping optional home-cache cleanup; remaining cache files are retained.', flush=True)
                journal['cleanup'] = dict(status='skipped', reason='Explicit --skip-cleanup',
                    previous_cleanup_record=(str(root/'cache_cleanup.json')
                                             if (root/'cache_cleanup.json').exists() else None))
            else:
                journal['cleanup'] = run_cache_cleanup.cleanup(
                    journal['archive'], root, config.get('cleanup', {}), protected=[handoff.REPO])
            handoff.verify(root)
            journal['phase'] = 'cleaned'
            atomic_json(path, journal)
        for case in handoff.CASES:
            d = root/case
            p = d/'launch.json'
            if p.exists():
                record = json.loads(p.read_text())
                if record.get('code') != journal['code']:
                    raise RuntimeError(f'{case}: launch snapshot differs; refusing duplicate launch.')
                journal['launches'][case] = record
            elif case in journal['launches']:
                raise RuntimeError(f'{case}: launch outcome is uncertain. Inspect processes before retrying; no duplicate started.')
            else:
                live = run_cleanup.inspect(d)
                if live['jobs'] or live['processes']:
                    raise RuntimeError(f'{case}: existing work has no launch receipt; no duplicate started.')
                journal['launches'][case] = {'status': 'launching'}
                atomic_json(path, journal)
                journal['launches'][case] = handoff.start_inputs(
                    [str(d/'input.json')], code=journal['code'])[0]
            atomic_json(path, journal)
        journal['phase'] = 'launched'
        atomic_json(path, journal)
        return journal


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--replace', required=True, type=Path)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('fresh_run.json'))
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--skip-cleanup', action='store_true',
                        help='Retain remaining home caches and resume launch; all run/seed/launch checks still apply.')
    args = parser.parse_args()
    result = run(args.root, args.replace, json.loads(args.config.read_text()),
                 apply=args.apply, skip_cleanup=args.skip_cleanup)
    print(json.dumps({k: result[k] for k in ('phase', 'previous', 'archive', 'nodes', 'code', 'launches')}, indent=2))
