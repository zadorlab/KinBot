"""Register an existing ordinary PES launch for the controller-refresh workflow.

This command writes management receipts only. It never signals processes,
changes calculation inputs, or submits/cancels jobs. The subsequent refresh
still performs its complete process, scheduler, input and lock preflight.
"""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import socket

from kinbot.execution import atomic_json, file_sha256


def register(root, *, apply=False):
    root = Path(root).resolve(strict=True)
    marker = root/'ordinary_adoption.json'
    if marker.exists():
        record = json.loads(marker.read_text())
        if record.get('schema') != 1 or record.get('root') != str(root):
            raise ValueError('Ordinary launch registration does not match this run.')
        if record['phase'] == 'registered':
            if not all((root/name).is_file() for name in ('prepared.json', 'salvage.json')):
                raise RuntimeError('Registered run management receipts are missing.')
            return record
    else:
        if any((root/name).exists() for name in ('prepared.json', 'salvage.json',
                                                 'fresh_run.json', 'controller_refresh.json')):
            raise ValueError('Run already has management metadata; use its existing refresh workflow.')
        inputs = sorted(root.glob('*/input.json'))
        if not inputs or (root/'input.json').exists():
            raise ValueError('Specify the parent containing ordinary PES case directories.')
        observed = {}
        for inp in inputs:
            if inp.is_symlink() or inp.parent.is_symlink() or not inp.resolve().is_relative_to(root):
                raise ValueError(f'Input escapes the controlled run: {inp}')
            if not isinstance(json.loads(inp.read_text()), dict):
                raise ValueError(f'Invalid PES input: {inp}')
            status = inp.parent/'pes_status.json'
            if status.is_symlink():
                raise ValueError(f'Symlinked controller receipt: {status}')
            native = status.read_bytes()
            saved = json.loads(native)
            if (saved.get('schema') != 1 or type(saved.get('pid')) is not int or saved['pid'] <= 0
                    or not saved.get('started_utc') or saved.get('status') not in ('running', 'complete', 'failed')
                    or str(saved.get('host', '')).split('.')[0] != socket.gethostname().split('.')[0]):
                raise ValueError(f'Cannot establish an ordinary launch on this host: {status}')
            pid = inp.parent/'launch.pid'
            if pid.exists() and (pid.is_symlink() or int(pid.read_text().strip()) != saved['pid']):
                raise ValueError(f'Launch PID and PES receipt disagree: {inp.parent}')
            observed[str(status.relative_to(root))] = dict(sha256=hashlib.sha256(native).hexdigest(), record=saved)
        record = dict(schema=1, root=str(root), phase='registering', observed_launches=observed,
                      prepared=dict(inputs=[str(p) for p in inputs],
                                    sha256={str(p.relative_to(root)): file_sha256(p) for p in inputs}))
    # Exact input hashes bound during registration; original native data stays
    # untouched. A partial metadata write can be completed idempotently.
    for name, digest in record['prepared']['sha256'].items():
        if file_sha256(root/name) != digest:
            raise ValueError(f'Input changed during registration: {name}')
    if not apply:
        return {**record, 'phase': 'reviewed'}
    atomic_json(marker, record)
    history = dict(status='launched', inputs=record['prepared']['inputs'], code=None,
                   original_code_revision='unknown: ordinary launch without an immutable source snapshot',
                   ordinary_adoption=str(marker))
    for name, content in (('prepared.json', record['prepared']), ('salvage.json', history)):
        path = root/name
        if path.exists() and json.loads(path.read_text()) != content:
            raise RuntimeError(f'Management receipt collision: {path}')
        if not path.exists():
            atomic_json(path, content)
    record['phase'] = 'registered'
    atomic_json(marker, record)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--apply', action='store_true', help='Write management receipts; start/stop nothing.')
    args = parser.parse_args(argv)
    root = args.root.resolve(strict=True)
    with root.with_suffix('.salvage.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = register(root, apply=args.apply)
    print(json.dumps(dict(phase=result['phase'], root=str(root),
                         cases=result['prepared']['inputs'], processes_changed=False,
                         jobs_changed=False), indent=2))


if __name__ == '__main__':
    main()
