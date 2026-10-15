"""Snapshot, pilot, run and collect arbitrary Molpro recon configurations.

Independent probes share one campaign-wide Slurm cap. Native pilots gate
full panels; checkpoint dependencies are retained. No SSH is performed.
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
import signal
import subprocess
import sys
import tarfile
import time

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


def stop_controller(root, *, timeout=30., proc=Path('/proc')):
    """Stop only this Linux campaign controller, preserving Slurm allocations."""
    root = Path(root).expanduser().resolve(strict=True)
    if not (root / 'campaign.json').is_file() or not (root / '.kinbot_run.json').is_file():
        raise ValueError('Expected a prepared orbital campaign.')
    with (root / '.campaign.lock').open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            saved = json.loads((root / 'campaign.json').read_text())
            pid = saved.get('controller_pid')
            if not isinstance(pid, int) or pid <= 1:
                raise RuntimeError('Busy campaign has no identified controller; inspect before stopping.')
            process = proc / str(pid)
            command = (process / 'cmdline').read_bytes().split(b'\0')
            signature = [b'-m', b'kinbot.orbital_campaign', b'run']
            if (process.stat().st_uid != os.getuid() or (process / 'cwd').resolve() != root
                    or not any(command[i:i + 3] == signature for i in range(len(command) - 2))):
                raise RuntimeError('Controller PID identity does not match this campaign; no signal sent.')
            os.kill(pid, signal.SIGTERM)
            deadline = time.monotonic() + timeout
            while True:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    if time.monotonic() >= deadline:
                        raise RuntimeError('Controller has not released the campaign lock; do not start another.')
                    time.sleep(.2)
        saved = json.loads((root / 'campaign.json').read_text())
        if saved.get('status') == 'running':
            saved.update(status='controller_stopped', error=None)
            saved.pop('controller_pid', None)
            atomic_json(root / 'campaign.json', saved)
        return {'status': 'controller_stopped', 'slurm_jobs_preserved': True}


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


def run(root, *, max_jobs=10, max_refinements=0, pilot_only=False, probe_method=None):
    root = Path(root).expanduser().resolve()
    if isinstance(max_jobs, bool) or not isinstance(max_jobs, int) or max_jobs < 1:
        raise ValueError('max_jobs must be positive.')
    if isinstance(max_refinements, bool) or not isinstance(max_refinements, int) or max_refinements < 0:
        raise ValueError('max_refinements must be nonnegative.')
    if probe_method not in (None, 'casci', 'casscf'):
        raise ValueError('probe_method must be casci or casscf.')
    with _lock(root):
        previous_directory = Path.cwd()
        os.chdir(root)  # run-scoped cleanup can identify the live controller
        path = root / 'campaign.json'
        campaign = json.loads(path.read_text())
        campaign.update(status='running', controller_pid=os.getpid(), max_jobs=max_jobs,
                        selected_probe_method=probe_method)
        campaign.pop('error', None)
        atomic_json(path, campaign)
        try:
            preflight(root)
            selected = []
            for record in campaign['systems']:
                _validate_inputs(record)
                method = json.loads(Path(record['configuration']).read_text()).get('probe_method', 'casci')
                if probe_method is None or method == probe_method:
                    selected.append(record)
            if not selected:
                raise ValueError(f'No systems request probe_method={probe_method!r}.')
            if probe_method is not None:
                queue = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                       text=True, check=True, timeout=30)
                live = set(queue.stdout.split())
                for record in campaign['systems']:
                    if record in selected:
                        continue
                    for key in ('pilot_manifest', 'manifest'):
                        manifest = json.loads(Path(record[key]).read_text())
                        jobs = manifest['jobs'] + manifest.get('excluded_jobs', [])
                        if any(str(job.get('job_id')) in live or job.get('submission_status') == 'submitting'
                               for job in jobs):
                            raise RuntimeError('Stop/reconcile unselected allocations before a scoped run.')
                        # Refresh old completion labels using current validation,
                        # but neither submit nor retry an unselected control.
                        recon.audit(record[key])
            all_manifests = [record[key] for record in selected for key in ('pilot_manifest', 'manifest')]
            # Independent pilots share the cap. All must pass before opening
            # full panels, which then need no second one-job startup gate.
            for record in selected:
                print(f'Native pilot: {record["title"]}', flush=True)
            pilots = recon.run_many([record['pilot_manifest'] for record in selected],
                                    max_jobs=max_jobs, capacity_manifests=all_manifests)
            for report in pilots:
                if not all(job['status'] == 'complete' for job in report['jobs']):
                    raise RuntimeError('Native reconnaissance pilot did not complete.')
            if not pilot_only:
                pending = selected
                for attempt in range(max_refinements + 1):
                    for record in pending:
                        print(f'Orbital reconnaissance: {record["title"]}', flush=True)
                    reports = recon.run_many([record['manifest'] for record in pending],
                        max_jobs=max_jobs, startup_checked=True, capacity_manifests=all_manifests)
                    refine = []
                    for record, report in zip(pending, reports):
                        record['status'] = report['status']
                        print(f'{record["title"]}: {report["status"]}', flush=True)
                        if report['status'] == 'coverage_refinement_required':
                            refine.append(record)
                    atomic_json(path, campaign)
                    if not refine or attempt == max_refinements:
                        break
                    for record in refine:
                        recon.refine(record['manifest'])
                    pending = refine
            campaign['status'] = 'pilot_complete' if pilot_only else 'complete_review_required'
            if probe_method is not None:
                campaign['status'] = 'selected_' + campaign['status']
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
    draining = False
    for record in data['systems']:
        item = {'title': record['title']}
        if data.get('selected_probe_method'):
            method = json.loads(Path(record['configuration']).read_text()).get('probe_method', 'casci')
            item['selected'] = method == data['selected_probe_method']
        for key in ('pilot_manifest', 'manifest'):
            manifest = json.loads(Path(record[key]).read_text())
            item[key] = {'status': manifest['status'], 'jobs': len(manifest['jobs']),
                         'complete': sum(job.get('status') == 'complete' for job in manifest['jobs'])}
            failures = [{'name': job['name'], 'directory': job['directory'],
                         'job_id': job.get('job_id'), 'reason': job['terminal_failure']}
                        for job in manifest['jobs'] if job.get('terminal_failure')]
            if failures:
                item[key].update(failed=len(failures), failures=failures)
                draining |= item.get('selected', True)
            if manifest.get('excluded_jobs'):
                item[key]['excluded_jobs'] = len(manifest['excluded_jobs'])
                item[key]['excluded_bases'] = manifest['excluded_bases']
        systems.append(item)
    status = 'draining_after_failure' if draining and data['status'] == 'running' else data['status']
    return {'status': status, 'max_jobs': data.get('max_jobs'),
            'selected_probe_method': data.get('selected_probe_method'),
            'error': data.get('error'), 'systems': systems}


def retry_failed(root, *, include_interrupted=False):
    """Reconcile failed native jobs without discarding successful probes."""
    root = Path(root).expanduser().resolve()
    with _lock(root):
        campaign = json.loads((root / 'campaign.json').read_text())
        results = []
        for record in campaign['systems']:
            _validate_inputs(record)
            for key in ('pilot_manifest', 'manifest'):
                results.append(recon.retry_failed(record[key], include_interrupted=include_interrupted))
        campaign['status'] = 'prepared_retry'
        campaign.pop('error', None)
        atomic_json(root / 'campaign.json', campaign)
        return results


@contextmanager
def _stopped_snapshot(root):
    """Hold campaign/manifests stable and reject still-owned native outputs."""
    with _lock(root), ExitStack() as locks:
        data = json.loads((root / 'campaign.json').read_text())
        manifests, jobs = {}, []
        for record in data['systems']:
            for key in ('pilot_manifest', 'manifest'):
                path = Path(record[key]).resolve()
                if path in manifests:
                    continue
                locks.enter_context(recon._manifest_lock(path.parent))
                manifest = json.loads(path.read_text())
                manifests[path] = manifest
                jobs.extend(manifest['jobs'] + manifest.get('excluded_jobs', []))
        if any(job.get('submission_status') == 'submitting' for job in jobs):
            raise RuntimeError('Reconcile uncertain submissions before exporting or reviewing this campaign.')
        recorded = {str(job['job_id']) for job in jobs if job.get('job_id')}
        if recorded:
            queue = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                   text=True, check=True, timeout=30)
            if recorded & set(queue.stdout.split()):
                raise RuntimeError('Wait for campaign allocations to leave Slurm before exporting or reviewing.')
        yield data, manifests


def collect(root, destination):
    """Archive review evidence; large binary WFU records remain on the cluster."""
    root = Path(root).expanduser().resolve()
    destination = Path(destination).expanduser().resolve()
    if root == destination or root in destination.parents:
        raise ValueError('Collection archive must be outside the campaign tree.')
    with _stopped_snapshot(root) as (data, _):
        paths = {root / 'campaign.json', root / 'review_plan.json', root / 'review_index.json'}
        paths.update(root.glob('review_report*.json'))
        for record in data['systems']:
            paths.update(Path(record['configuration']).parent.glob('*.json'))
            for key in ('pilot_manifest', 'manifest'):
                # _stopped_snapshot already holds every manifest lock.
                manifest = recon._audit(record[key])
                paths.add(Path(record[key]))
                for job in manifest['jobs'] + manifest.get('excluded_jobs', []):
                    directory = Path(job.get('archived_directory', job['directory']))
                    name = job['name']
                    paths.update(path for path in directory.glob(name + '.*')
                                 if path.suffix in ('.inp', '.out', '.stdout', '.stderr', '.molden', '.xml', '.bin'))
                    paths.add(directory / 'execution.json')
                    for attempts in (directory / 'attempts', directory / '.kinbot_slurm_io_retries'):
                        if attempts.is_dir():
                            paths.update(path for path in attempts.rglob('*') if path.is_file()
                                and (path.suffix in ('.inp', '.out', '.stdout', '.stderr', '.molden', '.xml', '.bin')
                                     or path.name in ('execution.json', 'stdout.txt', 'stderr.txt')))
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


def archive_excluded(root):
    """Move excluded basis directories out of the active tree, keeping evidence.

    Original paths remain provenance; archived_directory locates the retained
    files for collection. Repeating after an interrupted rename is safe.
    """
    root = Path(root).expanduser().resolve()
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
        moves = []
        for path, manifest in manifests:
            all_jobs = manifest['jobs'] + manifest.get('excluded_jobs', [])
            if any(str(job.get('job_id')) in live for job in all_jobs):
                raise RuntimeError('Wait for all campaign allocations to leave before archiving.')
            if any(job.get('submission_status') == 'submitting' for job in all_jobs):
                raise RuntimeError('Reconcile uncertain Slurm submissions before archiving.')
            for job in manifest.get('excluded_jobs', []):
                source = Path(job['directory'])
                destination = root / 'archived_excluded' / path.parent.name / source.name
                if (source.is_symlink() or source.parent.resolve() != path.parent.resolve()
                        or root not in source.resolve().parents):
                    raise ValueError('Excluded job directory is outside its manifest tree.')
                if job.get('archived_directory') not in (None, str(destination)):
                    raise ValueError('Excluded job has a different archive location.')
                if root not in destination.resolve().parents:
                    raise ValueError('Excluded archive destination is outside the campaign tree.')
                if destination.is_symlink() or source.exists() and destination.exists():
                    raise ValueError('Excluded archive destination already exists; preserve both copies.')
                if not source.is_dir() and not destination.is_dir():
                    raise FileNotFoundError(f'Excluded job artifacts are missing: {source}')
                moves.append((path, manifest, job, source, destination))
        for path, manifest, job, source, destination in moves:
            if source.is_dir():
                destination.parent.mkdir(parents=True, exist_ok=True)
                source.rename(destination)
            if job.get('archived_directory') != str(destination):
                job['archived_directory'] = str(destination)
                atomic_json(path, manifest)
        return {'archive': str(root / 'archived_excluded'),
                'excluded_directories': len(moves), 'files_deleted': 0}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    create = commands.add_parser('prepare')
    create.add_argument('root')
    create.add_argument('configurations', nargs='+')
    drive = commands.add_parser('run')
    drive.add_argument('root')
    drive.add_argument('--max-jobs', type=int, default=10)
    drive.add_argument('--max-refinements', type=int, default=0)
    drive.add_argument('--pilot-only', action='store_true')
    drive.add_argument('--probe-method', choices=['casci', 'casscf'],
                       help='Run only this probe method; retain other systems as unresolved controls.')
    report = commands.add_parser('progress')
    report.add_argument('root')
    check = commands.add_parser('preflight')
    check.add_argument('root')
    stop = commands.add_parser('stop-controller', help='Stop only the controller, keeping Slurm jobs running.')
    stop.add_argument('root')
    retry = commands.add_parser('retry-failed')
    retry.add_argument('root')
    retry.add_argument('--include-interrupted', action='store_true',
                       help='Also archive/retry stopped known allocations without a final execution record.')
    archive = commands.add_parser('archive-excluded')
    archive.add_argument('root')
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
                     max_refinements=args.max_refinements, pilot_only=args.pilot_only,
                     probe_method=args.probe_method)
    elif args.action == 'collect':
        result = collect(args.root, args.destination)
    elif args.action == 'preflight':
        result = preflight(args.root)
    elif args.action == 'stop-controller':
        result = stop_controller(args.root)
    elif args.action == 'retry-failed':
        result = retry_failed(args.root, include_interrupted=args.include_interrupted)
    elif args.action == 'archive-excluded':
        result = archive_excluded(args.root)
    elif args.action == 'exclude-bases':
        result = exclude_bases(args.root, args.basis)
    else:
        result = progress(args.root)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
