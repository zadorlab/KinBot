"""Snapshot, pilot, run and collect arbitrary Molpro recon configurations.

System controllers run sequentially and wait for their Slurm allocations to
leave, so max_jobs is a campaign-wide one-node job cap. No SSH is performed.
"""
import argparse
from contextlib import contextmanager, ExitStack
from copy import deepcopy
import fcntl
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tarfile

from kinbot.execution import atomic_json, file_sha256
from kinbot import molpro_reconnaissance as recon


@contextmanager
def _lock(root):
    with (root / '.campaign.lock').open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('Another controller owns this orbital campaign.') from error
        yield


def prepare(root, configurations):
    root = Path(root).expanduser().resolve()
    if root.exists():
        raise FileExistsError(f'{root}: choose a fresh campaign path, or resume it.')
    root.mkdir(parents=True)
    from kinbot.run_format import ensure_current_run
    ensure_current_run(root, create=True)
    records = []
    for index, source in enumerate(configurations):
        source = Path(source).expanduser().resolve()
        data = json.loads(source.read_text())
        folder = root / 'inputs' / f'{index + 1:02d}'
        folder.mkdir(parents=True)
        if data.get('kinbot_input'):
            linked = (source.parent / data['kinbot_input']).resolve()
            snapshot = folder / 'electronic_context.json'
            snapshot.write_bytes(linked.read_bytes())
            data['kinbot_input'] = str(snapshot)
        snapshot = folder / 'reconnaissance_input.json'
        atomic_json(snapshot, data)
        loaded = recon._load(snapshot)
        pilot = deepcopy(data)
        pilot.pop('coordinates', None)
        pilot.pop('scan_policy', None)
        # A native pilot checks every basis/channel's requested spin/root
        # block before launching its full radial/orientation campaign.
        first = {}
        for anchor in loaded['anchors']:
            if anchor.get('_channel'):
                first.setdefault(anchor['_channel'], anchor)
        pilot['anchors'] = list(first.values()) or loaded['anchors'][:1]
        pilot['title'] += '_pilot'
        pilot_path = folder / 'pilot_input.json'
        atomic_json(pilot_path, pilot)
        pilot_manifest = recon.prepare(pilot_path, root / f'{index + 1:02d}_pilot')
        full_manifest = recon.prepare(snapshot, root / f'{index + 1:02d}_recon')
        records.append({'configuration': str(snapshot), 'sha256': file_sha256(snapshot),
            'context_sha256': (file_sha256(data['kinbot_input']) if data.get('kinbot_input') else None),
            'pilot_manifest': str(pilot_manifest), 'manifest': str(full_manifest),
            'title': data['title']})
    atomic_json(root / 'campaign.json', {'schema': 1, 'status': 'prepared', 'systems': records,
                                        'python': sys.executable})
    return root / 'campaign.json'


def _validate_inputs(record):
    source = Path(record['configuration'])
    if file_sha256(source) != record['sha256']:
        raise ValueError('Campaign configuration changed; use a fresh campaign.')
    data = json.loads(source.read_text())
    if data.get('kinbot_input') and file_sha256(data['kinbot_input']) != record['context_sha256']:
        raise ValueError('Campaign electronic context changed; use a fresh campaign.')


def preflight(root):
    """Check the prepared module setup and scheduler without submitting work."""
    root = Path(root).expanduser().resolve()
    campaign = json.loads((root / 'campaign.json').read_text())
    for record in campaign['systems']:
        _validate_inputs(record)
    missing = [program for program in ('bash', 'sbatch', 'squeue')
               if not shutil.which(program)]
    if missing:
        raise RuntimeError('Required launch commands unavailable: ' + ', '.join(missing))
    setups = []
    for record in campaign['systems']:
        for key in ('pilot_manifest', 'manifest'):
            setup = Path(record[key]).parent / 'site_setup.sh'
            script = ('set -eo pipefail\nexport KINBOT_BACKEND=molpro\n'
                      f'source {shlex.quote(str(setup))}\n'
                      'command -v molpro\n')
            result = subprocess.run(['bash', '-c', script], capture_output=True,
                                    text=True, check=False, timeout=60)
            if result.returncode:
                raise RuntimeError(f'{setup}: Molpro site setup failed: '
                                   + (result.stderr.strip() or result.stdout.strip()))
            setups.append({'setup': str(setup), 'molpro': result.stdout.strip()})
    subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                   text=True, check=True, timeout=30)
    return {'status': 'ready_for_native_pilots', 'setups': setups}


def run(root, *, max_jobs=10, max_refinements=3, pilot_only=False):
    root = Path(root).expanduser().resolve()
    if isinstance(max_jobs, bool) or not isinstance(max_jobs, int) or max_jobs < 1:
        raise ValueError('max_jobs must be positive.')
    if isinstance(max_refinements, bool) or not isinstance(max_refinements, int) or max_refinements < 0:
        raise ValueError('max_refinements must be nonnegative.')
    with _lock(root):
        previous_directory = Path.cwd()
        os.chdir(root)  # run-scoped cleanup can identify the live controller
        path = root / 'campaign.json'
        campaign = json.loads(path.read_text())
        campaign.update(status='running', controller_pid=os.getpid(), max_jobs=max_jobs)
        campaign.pop('error', None)
        atomic_json(path, campaign)
        try:
            preflight(root)
            # Complete every system's native pilot first. No full recon jobs
            # start if any requested basis/root block fails its pilot.
            for record in campaign['systems']:
                _validate_inputs(record)
                print(f'Native pilot: {record["title"]}', flush=True)
                report = recon.run(record['pilot_manifest'], max_jobs=max_jobs)
                if not all(job['status'] == 'complete' for job in report['jobs']):
                    raise RuntimeError('Native reconnaissance pilot did not complete.')
            if not pilot_only:
                for record in campaign['systems']:
                    print(f'Orbital reconnaissance: {record["title"]}', flush=True)
                    for attempt in range(max_refinements + 1):
                        report = recon.run(record['manifest'], max_jobs=max_jobs)
                        if report['status'] != 'coverage_refinement_required' or attempt == max_refinements:
                            break
                        recon.refine(record['manifest'])
                    record['status'] = report['status']
                    print(f'{record["title"]}: {report["status"]}', flush=True)
                    atomic_json(path, campaign)
            campaign['status'] = 'pilot_complete' if pilot_only else 'complete_review_required'
            campaign.pop('error', None)
        except BaseException as error:
            campaign.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'failed',
                            error=str(error))
            raise
        finally:
            campaign.pop('controller_pid', None)
            atomic_json(path, campaign)
            os.chdir(previous_directory)
        return campaign


def progress(root):
    """Read progress without fighting a live controller's manifest lock."""
    root = Path(root).expanduser().resolve()
    data = json.loads((root / 'campaign.json').read_text())
    systems = []
    for record in data['systems']:
        item = {'title': record['title']}
        for key in ('pilot_manifest', 'manifest'):
            manifest = json.loads(Path(record[key]).read_text())
            item[key] = {'status': manifest['status'], 'jobs': len(manifest['jobs']),
                         'complete': sum(job.get('status') == 'complete' for job in manifest['jobs'])}
            if manifest.get('excluded_jobs'):
                item[key]['excluded_jobs'] = len(manifest['excluded_jobs'])
                item[key]['excluded_bases'] = manifest['excluded_bases']
        systems.append(item)
    return {'status': data['status'], 'max_jobs': data.get('max_jobs'),
            'error': data.get('error'), 'systems': systems}


def retry_failed(root):
    """Reconcile failed native jobs without discarding successful probes."""
    root = Path(root).expanduser().resolve()
    with _lock(root):
        campaign = json.loads((root / 'campaign.json').read_text())
        results = []
        for record in campaign['systems']:
            _validate_inputs(record)
            for key in ('pilot_manifest', 'manifest'):
                results.append(recon.retry_failed(record[key]))
        campaign['status'] = 'prepared_retry'
        campaign.pop('error', None)
        atomic_json(root / 'campaign.json', campaign)
        return results


def collect(root, destination):
    """Archive review evidence; large binary WFU records remain on the cluster."""
    root = Path(root).expanduser().resolve()
    destination = Path(destination).expanduser().resolve()
    if root == destination or root in destination.parents:
        raise ValueError('Collection archive must be outside the campaign tree.')
    with _lock(root):
        data = json.loads((root / 'campaign.json').read_text())
        paths = {root / 'campaign.json'}
        for record in data['systems']:
            paths.update(Path(record['configuration']).parent.glob('*.json'))
            for key in ('pilot_manifest', 'manifest'):
                manifest = recon.audit(record[key])
                paths.add(Path(record[key]))
                for job in manifest['jobs'] + manifest.get('excluded_jobs', []):
                    directory, name = Path(job['directory']), job['name']
                    paths.update(path for path in directory.glob(name + '.*')
                                 if path.suffix in ('.inp', '.out', '.stdout', '.stderr', '.molden', '.xml', '.bin'))
                    paths.add(directory / 'execution.json')
                    attempts = directory / 'attempts'
                    if attempts.is_dir():
                        paths.update(path for path in attempts.rglob('*') if path.is_file()
                            and (path.suffix in ('.inp', '.out', '.stdout', '.stderr', '.molden', '.xml', '.bin')
                                 or path.name == 'execution.json'))
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(destination, 'w:gz') as archive:
            for path in sorted(paths):
                if path.is_file():
                    archive.add(path, arcname=str(path.relative_to(root)), recursive=False)
    return str(destination)


def exclude_bases(root, bases):
    """Reduce a stopped campaign without modifying/deleting its evidence.

    Exact native basis identifiers are explicit user choices. Every manifest
    is checked before any mutation; the selection persists during refinement.
    Repeating the command is idempotent, including after an interrupted write.
    """
    root = Path(root).expanduser().resolve()
    requested = set(bases)
    if not requested or any(not isinstance(basis, str) or not basis.strip() for basis in requested):
        raise ValueError('Specify at least one exact basis name.')
    with _lock(root), ExitStack() as locks:
        data = json.loads((root / 'campaign.json').read_text())
        manifests = []
        for record in data['systems']:
            _validate_inputs(record)
            for key in ('pilot_manifest', 'manifest'):
                path = Path(record[key])
                locks.enter_context(recon._manifest_lock(path.parent))
                manifests.append((path, json.loads(path.read_text())))
        queue = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                               text=True, check=True, timeout=30)
        live = set(queue.stdout.split())
        known = {job['basis'] for _, manifest in manifests
                 for job in manifest['jobs'] + manifest.get('excluded_jobs', [])}
        if requested - known:
            raise ValueError('Unknown reconnaissance bases: ' + ', '.join(sorted(requested - known)))
        for _, manifest in manifests:
            if any(str(job.get('job_id')) in live for job in manifest['jobs']):
                raise RuntimeError('Wait for all campaign allocations to leave before changing bases.')
            if any(job.get('submission_status') == 'submitting' for job in manifest['jobs']):
                raise RuntimeError('Reconcile uncertain Slurm submissions before changing bases.')
            if not any(job['basis'] not in requested for job in manifest['jobs']):
                raise ValueError('Basis selection would leave a manifest empty.')
        reports = []
        for path, manifest in manifests:
            removed = [job for job in manifest['jobs'] if job['basis'] in requested]
            manifest['jobs'] = [job for job in manifest['jobs'] if job['basis'] not in requested]
            manifest.setdefault('excluded_jobs', []).extend(removed)
            manifest['excluded_bases'] = sorted(set(manifest.get('excluded_bases', [])) | requested)
            manifest['status'] = 'prepared_basis_selection'
            atomic_json(path, manifest)
            audited = recon._audit(path)
            reports.append({'manifest': str(path), 'retained': len(manifest['jobs']),
                            'newly_excluded': len(removed),
                            'completed_preserved': sum(job['status'] == 'complete' for job in audited['jobs'])})
        data['status'] = 'prepared_basis_selection'
        data.pop('error', None)
        atomic_json(root / 'campaign.json', data)
        return reports


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    create = commands.add_parser('prepare')
    create.add_argument('root')
    create.add_argument('configurations', nargs='+')
    drive = commands.add_parser('run')
    drive.add_argument('root')
    drive.add_argument('--max-jobs', type=int, default=10)
    drive.add_argument('--max-refinements', type=int, default=3)
    drive.add_argument('--pilot-only', action='store_true')
    report = commands.add_parser('progress')
    report.add_argument('root')
    check = commands.add_parser('preflight')
    check.add_argument('root')
    retry = commands.add_parser('retry-failed')
    retry.add_argument('root')
    selection = commands.add_parser('exclude-bases')
    selection.add_argument('root')
    selection.add_argument('--basis', action='append', required=True,
                           help='Exact native basis identifier to exclude; repeat as needed.')
    bundle = commands.add_parser('collect')
    bundle.add_argument('root')
    bundle.add_argument('destination')
    args = parser.parse_args(argv)
    if args.action == 'prepare':
        result = str(prepare(args.root, args.configurations))
    elif args.action == 'run':
        result = run(args.root, max_jobs=args.max_jobs,
                     max_refinements=args.max_refinements, pilot_only=args.pilot_only)
    elif args.action == 'collect':
        result = collect(args.root, args.destination)
    elif args.action == 'preflight':
        result = preflight(args.root)
    elif args.action == 'retry-failed':
        result = retry_failed(args.root)
    elif args.action == 'exclude-bases':
        result = exclude_bases(args.root, args.basis)
    else:
        result = progress(args.root)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
