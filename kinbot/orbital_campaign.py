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


def load_settings(root, settings=None):
    """Load optional, editable execution settings without changing scientific inputs."""
    path = Path(settings).expanduser().resolve() if settings else Path(root) / 'campaign_settings.json'
    if settings is None and not path.is_file():
        return {}
    data = json.loads(path.read_text())
    allowed = {'schema', 'max_jobs', 'max_refinements', 'pilot_only', 'probe_method',
               'keep_going', 'allocation_retries', 'startup_retries', 'io_retries',
               'memory_mb', 'excluded_nodes', 'reference_solver_options', 'scratch_root'}
    if not isinstance(data, dict) or type(data.get('schema')) is not int or data['schema'] != 1 or set(data) - allowed:
        raise ValueError('Unknown campaign settings or schema; expected schema 1.')
    for key in ('max_jobs', 'max_refinements', 'allocation_retries', 'startup_retries',
                'io_retries', 'memory_mb'):
        if key not in data:
            continue
        minimum = 1 if key in ('max_jobs', 'memory_mb') else 0
        if type(data[key]) is not int or data[key] < minimum:
            raise ValueError(f'campaign settings {key} must be an integer >= {minimum}.')
        if key.endswith('_retries') and data[key] > 10:
            raise ValueError(f'campaign settings {key} must be from 0 to 10.')
    for key in ('pilot_only', 'keep_going'):
        if key in data and type(data[key]) is not bool:
            raise ValueError(f'campaign settings {key} must be boolean.')
    if data.get('probe_method') not in (None, 'casci', 'casscf'):
        raise ValueError('campaign settings probe_method must be casci, casscf or null.')
    nodes = data.get('excluded_nodes', [])
    if not isinstance(nodes, list) or any(not isinstance(n, str) or not recon._SAFE.fullmatch(n) for n in nodes):
        raise ValueError('excluded_nodes must list individual Slurm node names.')
    if 'scratch_root' in data and (not isinstance(data['scratch_root'], str)
            or not data['scratch_root'].strip() or not Path(data['scratch_root']).is_absolute()):
        raise ValueError('scratch_root must be an absolute directory path.')
    controls = data.get('reference_solver_options', {})
    recon.validate_mr_correlation_options(controls)
    if set(controls) - {'max_reference_iterations', 'max_reference_vectors'}:
        raise ValueError('reference_solver_options supports reference iterations and vectors only; '
                         'change scientific controls in a new reconnaissance input.')
    return data


@contextmanager
def _settings_environment(settings):
    """Explicit JSON resource/retry settings override environment for this controller."""
    changes = {}
    if 'memory_mb' in settings:
        changes.update(SBATCH_MEM_PER_NODE=str(settings['memory_mb']),
                       SBATCH_MEM_PER_CPU=None, SBATCH_MEM_PER_GPU=None)
    if 'io_retries' in settings:
        changes['KINBOT_SLURM_IO_RETRIES'] = str(settings['io_retries'])
    if 'scratch_root' in settings:
        changes['KINBOT_MOLPRO_SCRATCH'] = settings['scratch_root']
    previous = {key: os.environ.get(key) for key in changes}
    try:
        for key, value in changes.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


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


def append(root, addition):
    """Adopt an unstarted nested campaign under the parent's single job cap.

    The parent controller must be stopped; existing allocations keep running.
    Child configurations, seeds and review provenance remain immutable.
    """
    root, addition = (Path(p).expanduser().resolve(strict=True) for p in (root, addition))
    if root not in addition.parents:
        raise ValueError('The addition must be prepared inside the parent campaign tree.')
    with _lock(root), _lock(addition), ExitStack() as locks:
        path = root / 'campaign.json'
        parent = json.loads(path.read_text())
        child = json.loads((addition / 'campaign.json').read_text())
        if not child.get('systems') or len({s['title'] for s in child['systems']}) != len(child['systems']):
            raise ValueError('Addition must contain distinct systems.')
        existing = {s['configuration']: s for s in parent['systems']}
        if all(all(existing.get(s['configuration'], {}).get(k) == s[k]
                   for k in ('configuration', 'sha256', 'context_sha256', 'pilot_manifest', 'manifest', 'title'))
               for s in child['systems']):
            return {'added': 0, 'systems': len(parent['systems'])}
        if child.get('adopted_by') not in (None, str(root)):
            raise ValueError('This addition is already owned by another campaign.')
        titles = {s['title'] for s in parent['systems']}
        if any(s['title'] in titles or s['configuration'] in existing for s in child['systems']):
            raise ValueError('Addition conflicts with existing campaign systems.')
        index_path = root / 'review_index.json'
        child_index = addition / 'review_index.json'
        if index_path.is_file() != child_index.is_file():
            raise ValueError('Both campaigns need review metadata, or neither: incomplete metadata breaks review exports.')
        if index_path.is_file():
            for owner, index_file in ((parent, index_path), (child, child_index)):
                indexed = [s['title'] for s in json.loads(index_file.read_text())['systems']]
                if len(indexed) != len(set(indexed)) or not {s['title'] for s in owner['systems']} <= set(indexed):
                    raise ValueError('Review metadata must identify every campaign system uniquely.')
        for record in child['systems']:
            _validate_inputs(record)
            for key in ('configuration', 'pilot_manifest', 'manifest'):
                if addition not in Path(record[key]).resolve().parents:
                    raise ValueError('Addition contains paths outside its own tree.')
            for key in ('pilot_manifest', 'manifest'):
                manifest_path = Path(record[key])
                locks.enter_context(recon._manifest_lock(manifest_path.parent))
                manifest = json.loads(manifest_path.read_text())
                if any(j.get('job_id') or j.get('submission_status') or
                       (Path(j['directory']) / 'execution.json').exists()
                       for j in manifest['jobs'] + manifest.get('excluded_jobs', [])):
                    raise ValueError('Only unstarted additions can be adopted.')
        # Metadata is published first. Repeating after an interrupted write
        # merges by title, so neither duplicate review entries nor jobs appear.
        if child_index.is_file():
            index = json.loads(index_path.read_text())
            entries = {s['title']: s for s in index['systems']}
            entries.update({s['title']: s for s in json.loads(child_index.read_text())['systems']})
            index['systems'] = list(entries.values())
            atomic_json(index_path, index)
        child['adopted_by'] = str(root)
        atomic_json(addition / 'campaign.json', child)
        parent['systems'].extend(child['systems'])
        parent.setdefault('additions', []).append(str(addition))
        atomic_json(path, parent)
        return {'added': len(child['systems']), 'systems': len(parent['systems'])}


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


def run(root, *, max_jobs=None, max_refinements=None, pilot_only=None, probe_method=None,
        keep_going=None, allocation_retries=None, settings=None):
    root = Path(root).expanduser().resolve()
    configured = load_settings(root, settings)
    max_jobs = configured.get('max_jobs', 10) if max_jobs is None else max_jobs
    max_refinements = configured.get('max_refinements', 0) if max_refinements is None else max_refinements
    pilot_only = configured.get('pilot_only', False) if pilot_only is None else pilot_only
    probe_method = configured.get('probe_method') if probe_method is None else probe_method
    keep_going = configured.get('keep_going', False) if keep_going is None else keep_going
    allocation_retries = configured.get('allocation_retries', 0) if allocation_retries is None else allocation_retries
    if isinstance(max_jobs, bool) or not isinstance(max_jobs, int) or max_jobs < 1:
        raise ValueError('max_jobs must be positive.')
    if isinstance(max_refinements, bool) or not isinstance(max_refinements, int) or max_refinements < 0:
        raise ValueError('max_refinements must be nonnegative.')
    if probe_method not in (None, 'casci', 'casscf'):
        raise ValueError('probe_method must be casci or casscf.')
    if type(allocation_retries) is not int or not 0 <= allocation_retries <= 10:
        raise ValueError('allocation_retries must be an integer from 0 to 10.')
    if configured.get('excluded_nodes'):
        refresh_exclusions(root, nodes=configured['excluded_nodes'], settings=settings)
    with _lock(root), _settings_environment(configured):
        previous_directory = Path.cwd()
        path = root / 'campaign.json'
        campaign = json.loads(path.read_text())
        if campaign.get('adopted_by'):
            raise ValueError('Run the owning campaign instead: ' + campaign['adopted_by'])
        # Validate resources for every panel before any submission or cwd change.
        for record in campaign['systems']:
            _validate_inputs(record)
            data = json.loads(Path(record['configuration']).read_text())
            if ('memory_mb' in configured and configured['memory_mb'] <
                    data.get('ppn', 8) * data.get('memory_mw', 1000) * 8):
                raise ValueError('campaign memory_mb does not cover all Molpro rank stacks.')
        if 'startup_retries' in configured or configured.get('reference_solver_options'):
            for manifest_path in sorted({Path(s[k]) for s in campaign['systems'] for k in ('pilot_manifest', 'manifest')}):
                with recon._manifest_lock(manifest_path.parent):
                    manifest = json.loads(manifest_path.read_text())
                    if 'startup_retries' in configured:
                        manifest['startup_retries'] = configured['startup_retries']
                    if 'reference_solver_options' in configured:
                        manifest['reference_solver_options'] = configured['reference_solver_options']
                    atomic_json(manifest_path, manifest)
        campaign['execution_settings'] = {**configured, 'max_jobs': max_jobs,
            'max_refinements': max_refinements, 'pilot_only': pilot_only, 'probe_method': probe_method,
            'keep_going': keep_going, 'allocation_retries': allocation_retries}
        os.chdir(root)  # run-scoped cleanup can identify the live controller
        campaign.update(status='running', controller_pid=os.getpid(), max_jobs=max_jobs,
                        selected_probe_method=probe_method, keep_going=keep_going,
                        allocation_retries=allocation_retries)
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
                                    max_jobs=max_jobs, capacity_manifests=all_manifests,
                                    **({'keep_going': True} if keep_going else {}),
                                    **({'allocation_retries': allocation_retries} if allocation_retries else {}))
            for report in pilots:
                if not all(job['status'] == 'complete' for job in report['jobs']):
                    raise RuntimeError('Native reconnaissance pilot did not complete.')
            if not pilot_only:
                pending = selected
                for attempt in range(max_refinements + 1):
                    for record in pending:
                        print(f'Orbital reconnaissance: {record["title"]}', flush=True)
                    reports = recon.run_many([record['manifest'] for record in pending],
                        max_jobs=max_jobs, startup_checked=True, capacity_manifests=all_manifests,
                        **({'keep_going': True} if keep_going else {}),
                        **({'allocation_retries': allocation_retries} if allocation_retries else {}))
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


def progress(root, *, failures_only=False):
    """Read progress without fighting a live controller's manifest lock."""
    root = Path(root).expanduser().resolve()
    data = json.loads((root / 'campaign.json').read_text())
    systems = []
    active_failures = []
    excluded_nodes, startup_retries, allocation_retries = set(), 0, 0
    for record in data['systems']:
        item = {'title': record['title']}
        if data.get('selected_probe_method'):
            method = json.loads(Path(record['configuration']).read_text()).get('probe_method', 'casci')
            item['selected'] = method == data['selected_probe_method']
        for key in ('pilot_manifest', 'manifest'):
            manifest = json.loads(Path(record[key]).read_text())
            excluded_nodes.update(manifest.get('slurm_excluded_nodes', []))
            startup_retries += sum(len(job.get('startup_retry_history', []))
                                   for job in manifest['jobs'])
            allocation_retries += sum(len(job.get('allocation_retry_history', []))
                                      for job in manifest['jobs'])
            item[key] = {'status': manifest['status'], 'jobs': len(manifest['jobs']),
                         'complete': sum(job.get('status') == 'complete' for job in manifest['jobs'])}
            failures = [{'name': job['name'], 'directory': job['directory'],
                         'job_id': job.get('job_id'), 'reason': job['terminal_failure']}
                        for job in manifest['jobs'] if job.get('terminal_failure')]
            if failures:
                item[key].update(failed=len(failures), failures=failures)
                if item.get('selected', True):
                    active_failures.extend({'system': record['title'], 'phase': key, **failure}
                                           for failure in failures)
            if manifest.get('excluded_jobs'):
                item[key]['excluded_jobs'] = len(manifest['excluded_jobs'])
                item[key]['excluded_bases'] = manifest['excluded_bases']
        systems.append(item)
    status = data['status']
    if active_failures and status == 'running':
        status = 'running_with_failures' if data.get('keep_going') else 'draining_after_failure'
    result = {'status': status, 'max_jobs': data.get('max_jobs'),
            'selected_probe_method': data.get('selected_probe_method'),
            'error': data.get('error') or (active_failures[0]['reason'] if active_failures else None),
            'failed_jobs': len(active_failures), 'startup_retries': startup_retries,
            'allocation_retries': allocation_retries,
            'slurm_excluded_nodes': sorted(excluded_nodes)}
    result.update({'failures': active_failures} if failures_only else {'systems': systems})
    return result


def refresh_exclusions(root, *, nodes=(), settings=None):
    """Apply recorded exclusions to this campaign's pending jobs only.

    Requires a stopped controller. Running jobs keep their allocations; this
    never cancels or requeues them. Preserve existing scheduler exclusions.
    """
    root = Path(root).expanduser().resolve()
    nodes = [*nodes, *load_settings(root, settings).get('excluded_nodes', [])]
    if any(not isinstance(node, str) or not recon._SAFE.fullmatch(node) for node in nodes):
        raise ValueError('Exclusions must be individual Slurm node names.')
    with _lock(root), ExitStack() as locks:
        campaign = json.loads((root / 'campaign.json').read_text())
        paths = sorted({Path(record[key]).resolve() for record in campaign['systems']
                        for key in ('pilot_manifest', 'manifest')})
        for path in paths:
            locks.enter_context(recon._manifest_lock(path.parent))
        manifests = [json.loads(path.read_text()) for path in paths]
        excluded = sorted(set(nodes) | {node for m in manifests for node in m.get('slurm_excluded_nodes', [])})
        if any(not isinstance(node, str) or not recon._SAFE.fullmatch(node) for node in excluded):
            raise ValueError('Invalid campaign exclusion node.')
        if nodes:
            for path, manifest in zip(paths, manifests):
                manifest['slurm_excluded_nodes'] = excluded
                atomic_json(path, manifest)
        reports, seen = [], set()
        for manifest in manifests:
            for job in manifest['jobs'] + manifest.get('excluded_jobs', []):
                identifier = str(job.get('job_id', ''))
                if not identifier or identifier in seen:
                    continue
                seen.add(identifier)
                record = recon.slurm_job_record(identifier)
                entry = {'job_id': identifier, 'state': record.get('JobState', 'unavailable')}
                reports.append(entry)
                if (not excluded or record.get('JobId') != identifier
                        or record.get('WorkDir') != job['directory'] or record.get('JobState') != 'PENDING'):
                    entry['action'] = 'preserved'
                    continue
                old = record.get('ExcNodeList', '')
                values = ([old] if old not in ('', '(null)', 'None') else []) + excluded
                value = ','.join(dict.fromkeys(values))
                result = subprocess.run(['scontrol', 'update', 'JobId=' + identifier, 'ExcNodeList=' + value],
                    capture_output=True, text=True, check=False, timeout=15)
                after = recon.slurm_job_record(identifier)
                if result.returncode:
                    if (after.get('JobId') == identifier and after.get('WorkDir') == job['directory']
                            and after.get('JobState') in {'RUNNING', *recon.SLURM_TERMINAL_STATES}):
                        entry['action'] = 'allocation_changed_during_update; preserved'
                        continue
                    raise RuntimeError(f'{identifier}: exclusion update failed: {result.stderr.strip()}')
                if after.get('JobState') in ('PENDING', 'RUNNING'):
                    if after.get('JobId') != identifier or after.get('WorkDir') != job['directory']:
                        raise RuntimeError(f'{identifier}: scheduler identity changed after exclusion update.')
                    actual = after.get('ExcNodeList', '')
                    if '[' in actual:
                        expanded = subprocess.run(['scontrol', 'show', 'hostnames', actual],
                            capture_output=True, text=True, check=True, timeout=15)
                        actual_nodes = set(expanded.stdout.split())
                    else:
                        actual_nodes = set(actual.split(','))
                    if not set(excluded) <= actual_nodes:
                        raise RuntimeError(f'{identifier}: Slurm did not retain exclusions: {actual!r}')
                    entry.update(action='updated', requested=value, scheduler_exclusions=actual)
                else:
                    entry.update(action='allocation_left_queue; update_unverified', requested=value)
        return {'excluded_nodes': excluded, 'jobs': reports}


def retry_startup_failures(root, *, settings=None):
    """Explicit recovery after infrastructure repair; preserve retry evidence.

    Reuse the controller's strict pre-shell classifier (scheduler identity,
    signal 53, no worker/guard artifacts, publication grace). Grant the configured
    startup allowance (default two) only to explicitly recovered stopped jobs.
    """
    root = Path(root).expanduser().resolve()
    allowance = load_settings(root, settings).get('startup_retries', 2)
    with _lock(root):
        campaign = json.loads((root / 'campaign.json').read_text())
        queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                 text=True, check=True, timeout=30)
        live = set(queued.stdout.split())
        reports, seen = [], set()
        for record in campaign['systems']:
            _validate_inputs(record)
            for key in ('pilot_manifest', 'manifest'):
                path = Path(record[key]).resolve()
                if path in seen:
                    continue
                seen.add(path)
                with recon._manifest_lock(path.parent):
                    manifest = json.loads(path.read_text())
                    manifest['startup_retries'] = allowance
                    recon._recover_slurm_startup_failures(path, manifest, live)
                    selected = {str(j['job_id']): j['name'] for j in manifest['jobs']
                        if j.get('startup_failure') and j.get('job_id')
                        and str(j['job_id']) not in live
                        and time.time() - j.get('missing_since', time.time()) > 60.}
                    if not selected:
                        continue
                    result = recon._retry_failed(path, include_interrupted=True,
                                                  job_ids=set(selected), live_job_ids=live)
                    manifest = json.loads(path.read_text())
                    for identifier, name in selected.items():
                        if name not in result['retry_ready']:
                            continue
                        job = next(j for j in manifest['jobs'] if j['name'] == name)
                        limit = len(job.get('startup_retry_history', [])) + allowance
                        job.setdefault('manual_startup_retry_history', []).append({
                            'job_id': identifier, 'time': time.time(),
                            'previous_limit': job.get('startup_retry_limit', allowance), 'new_limit': limit})
                        job['startup_retry_limit'] = limit
                        job.pop('startup_failure', None)
                    atomic_json(path, manifest)
                    reports.append(result)
        campaign['status'] = 'prepared_retry'
        campaign.pop('error', None)
        atomic_json(root / 'campaign.json', campaign)
        return {'manual_retries': reports, 'running_jobs_preserved': True}


def retry_failed(root, *, include_interrupted=False, job_ids=None, max_reference_vectors=None,
                 max_reference_iterations=None, settings=None):
    """Reconcile failed native jobs without discarding successful probes."""
    root = Path(root).expanduser().resolve()
    controls = load_settings(root, settings).get('reference_solver_options', {}).copy()
    controls.update({key: value for key, value in (
        ('max_reference_vectors', max_reference_vectors),
        ('max_reference_iterations', max_reference_iterations)) if value is not None})
    if controls:
        recon.validate_mr_correlation_options(controls)
        if not job_ids:
            raise ValueError('A numerical retry override requires explicit --job-id selection.')
    with _lock(root):
        campaign = json.loads((root / 'campaign.json').read_text())
        if job_ids is not None:
            job_ids = {str(identifier) for identifier in job_ids}
            if not job_ids or any(not identifier.isdigit() for identifier in job_ids):
                raise ValueError('Retry job IDs must be nonempty numeric Slurm IDs.')
        results = []
        for record in campaign['systems']:
            _validate_inputs(record)
            for key in ('pilot_manifest', 'manifest'):
                if job_ids is not None:
                    manifest = json.loads(Path(record[key]).read_text())
                    if not any(str(j.get('job_id')) in job_ids for j in manifest['jobs']):
                        continue
                results.append(recon.retry_failed(record[key], include_interrupted=include_interrupted,
                    **({'job_ids': job_ids} if job_ids is not None else {}),
                    **controls))
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
    with _stopped_snapshot(root) as (data, snapshots):
        paths = {root / 'campaign.json', root / 'campaign_settings.json', root / 'review_plan.json', root / 'review_index.json',
                 root / 'controller.log'}
        for addition in data.get('additions', []):
            paths.update(Path(addition) / name for name in ('campaign.json', 'review_plan.json', 'review_index.json'))
        # Native output cannot establish why Slurm killed/requeued a process.
        # Keep allocation and step accounting with the same stopped snapshot;
        # sites without sacct still get all native evidence and an explicit
        # accounting-unavailable record. Never query unrelated allocations.
        identifiers = set()
        for manifest in snapshots.values():
            for job in manifest['jobs'] + manifest.get('excluded_jobs', []):
                for attempt in [job, *job.get('attempt_history', []),
                                *job.get('startup_retry_history', [])]:
                    identifier = str(attempt.get('job_id', ''))
                    if identifier.isdigit():
                        identifiers.add(identifier)
        accounting = {'schema': 1, 'captured_unix_time': time.time(), 'queries': []}
        ordered = sorted(identifiers, key=int)
        for start in range(0, len(ordered), 200):
            command = ['sacct', '--jobs=' + ','.join(ordered[start:start + 200]),
                       '--parsable2', '--duplicates',
                       '--format=JobIDRaw,State,ExitCode,Reason,Elapsed,MaxRSS,NodeList,Restarts']
            query = {'command': command}
            try:
                result = subprocess.run(command, capture_output=True, text=True, timeout=30)
                query.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
            except (OSError, subprocess.TimeoutExpired) as error:
                query['error'] = str(error)
            accounting['queries'].append(query)
        accounting_path = root / 'scheduler_accounting.json'
        atomic_json(accounting_path, accounting)
        paths.add(accounting_path)
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
                                 if path.suffix in ('.inp', '.out', '.stdout', '.stderr', '.molden', '.xml', '.bin', '.sbatch'))
                    paths.add(directory / 'execution.json')
                    for attempts in (directory / 'attempts', directory / '.kinbot_slurm_io_retries'):
                        if attempts.is_dir():
                            paths.update(path for path in attempts.rglob('*') if path.is_file()
                                and (path.suffix in ('.inp', '.out', '.stdout', '.stderr', '.molden', '.xml', '.bin', '.sbatch')
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
    add = commands.add_parser('append', help='Adopt an unstarted nested campaign under the same job cap.')
    add.add_argument('root')
    add.add_argument('addition')
    drive = commands.add_parser('run')
    drive.add_argument('root')
    drive.add_argument('--max-jobs', type=int)
    drive.add_argument('--max-refinements', type=int)
    drive.add_argument('--pilot-only', action=argparse.BooleanOptionalAction, default=None)
    drive.add_argument('--keep-going', action=argparse.BooleanOptionalAction, default=None,
                       help='Finish independent manifests after a probe failure; still exit with failure.')
    drive.add_argument('--allocation-retries', type=int,
                       help='Bounded retries of finalized I/O failures or interrupted known workers (0–10).')
    drive.add_argument('--probe-method', choices=['casci', 'casscf'],
                       help='Run only this probe method; retain other systems as unresolved controls.')
    report = commands.add_parser('progress')
    report.add_argument('root')
    report.add_argument('--failures-only', action='store_true',
                        help='Show failed selected jobs without the full system inventory.')
    monitor = commands.add_parser('monitor', help='Read scheduler facts and live native log tails; never modify jobs.')
    monitor.add_argument('root')
    monitor.add_argument('--details', action='store_true', help='Include native stage and latest iteration/output.')
    monitor.add_argument('--json', action='store_true', help='Emit the read-only snapshot as JSON.')
    refresh = commands.add_parser('refresh-exclusions', help='Apply campaign exclusions to pending allocations only.')
    refresh.add_argument('root')
    refresh.add_argument('--node', action='append', default=[],
                         help='Persist an additional individual node exclusion across the campaign; repeat as needed.')
    startup = commands.add_parser('retry-startup-failures',
        help='Explicitly recover confirmed pre-shell I/O failures after infrastructure repair; preserve live jobs.')
    startup.add_argument('root')
    check = commands.add_parser('preflight')
    check.add_argument('root')
    stop = commands.add_parser('stop-controller', help='Stop only the controller, keeping Slurm jobs running.')
    stop.add_argument('root')
    retry = commands.add_parser('retry-failed')
    retry.add_argument('root')
    retry.add_argument('--include-interrupted', action='store_true',
                       help='Also archive/retry stopped known allocations without a final execution record.')
    retry.add_argument('--job-id', action='append',
                       help='Retry only these allocation IDs; preserve unrelated running jobs. Repeat as needed.')
    retry.add_argument('--max-reference-vectors', type=int,
                       help='Explicit numerical retry: MR internal Davidson capacity (MAXVI); requires --job-id.')
    retry.add_argument('--max-reference-iterations', type=int,
                       help='Explicit numerical retry: MR internal CI iteration ceiling (MAXITI); requires --job-id.')
    for command in (drive, refresh, startup, retry):
        command.add_argument('--settings', help='Execution settings JSON; defaults to ROOT/campaign_settings.json.')
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
    if args.action == 'monitor':
        from kinbot.orbital_monitor import snapshot, render
        result = snapshot(args.root)
        print(json.dumps(result, indent=2) if args.json else render(result, details=args.details))
        return 0
    if args.action == 'prepare':
        result = str(prepare(args.root, args.configurations))
    elif args.action == 'append':
        result = append(args.root, args.addition)
    elif args.action == 'run':
        result = run(args.root, max_jobs=args.max_jobs,
                     max_refinements=args.max_refinements, pilot_only=args.pilot_only,
                     probe_method=args.probe_method, keep_going=args.keep_going,
                     allocation_retries=args.allocation_retries, settings=args.settings)
    elif args.action == 'collect':
        result = collect(args.root, args.destination)
    elif args.action == 'preflight':
        result = preflight(args.root)
    elif args.action == 'stop-controller':
        result = stop_controller(args.root)
    elif args.action == 'refresh-exclusions':
        result = refresh_exclusions(args.root, nodes=args.node, settings=args.settings)
    elif args.action == 'retry-startup-failures':
        result = retry_startup_failures(args.root, settings=args.settings)
    elif args.action == 'retry-failed':
        result = retry_failed(args.root, include_interrupted=args.include_interrupted, job_ids=args.job_id,
                              max_reference_vectors=args.max_reference_vectors,
                              max_reference_iterations=args.max_reference_iterations, settings=args.settings)
    elif args.action == 'archive-excluded':
        result = archive_excluded(args.root)
    elif args.action == 'exclude-bases':
        result = exclude_bases(args.root, args.basis)
    else:
        result = progress(args.root, failures_only=args.failures_only)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
