"""Run, resume or read the fixed-chain confirmation using KinBot's executor.

No SSH, automatic chemical acceptance or generation of production selections.
The first point of each chain naturally gates its dependent points, avoiding
a second execution of the same point as a separate pilot.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import fcntl
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import tarfile

from kinbot.execution import atomic_json, file_sha256
from kinbot import molpro_reconnaissance as recon
from kinbot.orbital_campaign import _settings_environment
from kinbot.orbital_monitor import snapshot, render


def validated_record(root):
    record = json.loads((root / 'confirmation.json').read_text())
    if file_sha256(root / 'catalog.json') != record['catalog_sha256']:
        raise ValueError('Confirmation catalogue changed; prepare a new run.')
    for system in record['systems']:
        if file_sha256(system['configuration']) != system['configuration_sha256']:
            raise ValueError('Immutable confirmation configuration changed.')
    return record


def run(root):
    root = Path(root).resolve()
    path = root / 'confirmation.json'
    with (root / '.confirmation.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        record = validated_record(root)
        settings = record['execution']
        manifests = [system['manifest'] for system in record['systems']]
        original_directory = Path.cwd()
        previous_handler = signal.getsignal(signal.SIGTERM)
        def interrupted(signum, frame):
            raise KeyboardInterrupt('Controller interrupted; native jobs preserved.')
        signal.signal(signal.SIGTERM, interrupted)
        os.chdir(root)
        record.update(status='running', controller_pid=os.getpid(), max_jobs=settings['max_jobs'])
        record.pop('error', None)
        atomic_json(path, record)
        try:
            with _settings_environment(settings):
                results = recon.run_many(manifests, max_jobs=settings['max_jobs'],
                    keep_going=settings['keep_going'],
                    allocation_retries=settings['allocation_retries'])
            if not all(j.get('status') == 'complete' for m in results for j in m['jobs']):
                raise RuntimeError('Confirmation has unresolved jobs.')
            record['status'] = 'complete_review_required'
        except BaseException as error:
            record.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'failed',
                          error=str(error))
            raise
        finally:
            record.pop('controller_pid', None)
            atomic_json(path, record)
            os.chdir(original_directory)
            signal.signal(signal.SIGTERM, previous_handler)
    return record


def collect(root, destination):
    """Export a stopped native snapshot, including WFU and original seed bytes.

    A partial snapshot is useful diagnostic evidence, never production approval.
    Scheduler failure or live allocations prevent collection. No native result
    is removed, and an existing archive is never overwritten.
    """
    root = Path(root).resolve()
    destination = Path(destination).resolve()
    if destination == root or root in destination.parents:
        raise ValueError('Write the archive outside the confirmation directory.')
    if destination.exists():
        raise FileExistsError(destination)
    with (root / '.confirmation.lock').open('a+') as lock, ExitStack() as stack:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        record = validated_record(root)
        manifests = sorted({Path(s['manifest']) for s in record['systems']})
        for path in manifests:
            stack.enter_context(recon._manifest_lock(path.parent))
        queued = subprocess.run(['squeue', '-h', '-o', '%A'], check=True,
                                capture_output=True, text=True, timeout=30)
        live = set(queued.stdout.split())
        for path in manifests:
            for job in json.loads(path.read_text())['jobs']:
                if job.get('submission_status') == 'submitting':
                    raise RuntimeError('Reconcile uncertain submissions before collection.')
                if str(job.get('job_id')) in live:
                    raise RuntimeError('Native allocations remain live; collect after they finish.')
        audits = [recon._audit(path, live_job_ids=live) for path in manifests]
        # Export ordinary files only; links could point outside this evidence tree.
        paths = {}
        for path in sorted(root.rglob('*')):
            if path.is_symlink():
                raise ValueError('Evidence export refuses symbolic links: ' + str(path))
            if path.is_file() and not path.name.endswith('.lock'):
                paths[path.relative_to(root).as_posix()] = path
        catalog = json.loads((root / 'catalog.json').read_text())
        for case in catalog['cases']:
            seed = case['seed']
            path = Path(seed['path'])
            if file_sha256(path) != seed['sha256']:
                raise ValueError('Original reviewed seed is missing or changed: ' + str(path))
            name = 'reviewed_seeds/' + case['id'] + '.wfu'
            if name in paths:
                raise ValueError('Reserved seed archive name already exists: ' + name)
            paths[name] = path
        inventory = {'schema': 1, 'captured_utc': datetime.now(timezone.utc).isoformat(),
            'source_directory': str(root), 'production_approved': False,
            'complete_jobs': sum(j['status'] == 'complete' for m in audits for j in m['jobs']),
            'total_jobs': sum(len(m['jobs']) for m in audits),
            'files': [{'archive_path': name, 'source_path': str(path),
                       'sha256': file_sha256(path), 'bytes': path.stat().st_size}
                      for name, path in sorted(paths.items())]}
        destination.parent.mkdir(parents=True, exist_ok=True)
        # The inventory is inside the archive; the archive digest is a sidecar.
        with tarfile.open(destination, 'x:gz') as archive:
            for artifact in inventory['files']:
                path = paths[artifact['archive_path']]
                if file_sha256(path) != artifact['sha256']:
                    raise RuntimeError('Evidence changed during collection: ' + str(path))
                archive.add(path, arcname=artifact['archive_path'], recursive=False)
            payload = json.dumps(inventory, indent=2).encode()
            header = tarfile.TarInfo('export_inventory.json')
            header.size = len(payload)
            archive.addfile(header, io.BytesIO(payload))
        digest = file_sha256(destination)
        destination.with_name(destination.name + '.sha256').write_text(
            digest + '  ' + destination.name + '\n')
    return {'archive': str(destination), 'sha256': digest,
            'complete_jobs': inventory['complete_jobs'], 'total_jobs': inventory['total_jobs'],
            'production_approved': False}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['run', 'monitor', 'collect'])
    p.add_argument('root', type=Path)
    p.add_argument('--output', type=Path, help='Fresh archive outside the run tree, for collect.')
    args = p.parse_args()
    if args.command == 'monitor':
        print(render(snapshot(args.root, registry=args.root / 'confirmation.json'), details=True))
    elif args.command == 'run':
        print(json.dumps(run(args.root), indent=2))
    else:
        if args.output is None:
            p.error('collect requires --output')
        print(json.dumps(collect(args.root, args.output), indent=2))
