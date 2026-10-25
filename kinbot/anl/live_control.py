"""Owner-only dispatcher control with journalled archival and generation fences."""
from copy import deepcopy
import fcntl
import json
from pathlib import Path
import uuid

from kinbot import live_control as control
from kinbot.execution import atomic_json, file_sha256, slurm_job_record, SLURM_TERMINAL_STATES


def worker_load(graph, ident):
    """Unchanged allocations may start even if a controller died mid-edit."""
    from kinbot.anl import dispatch
    graph = Path(graph).resolve()
    # NFS implements flock through POSIX byte-range locks. A shared lock
    # requires a readable descriptor there (write-only append raises EBADF).
    with (graph/'.manifest.lock').open('a+') as handle:
        fcntl.flock(handle, fcntl.LOCK_SH)
        journal = graph/'.live_control_transaction.json'
        txn = json.loads(journal.read_text()) if journal.exists() else {}
        if txn.get('status') != 'applying':
            return dispatch._load_unlocked(graph)
        plan = txn['plan']
        if ident in plan['reset']:
            raise RuntimeError('Retired task generation cannot start during a live edit.')
        before = {}
        for name, digest in plan['before'].items():
            saved = Path(txn['archive'])/'graphs'/name
            source = saved if saved.exists() else graph/name
            if file_sha256(source) != digest:
                raise RuntimeError('Live-edit previous manifest verification failed.')
            before[name] = json.loads(source.read_text())
        return graph, before['workflow.json'], before['state.json']


def resume(graph):
    marker = Path(graph)/'.live_control_transaction.json'
    if not marker.exists():
        return
    transaction = json.loads(marker.read_text())
    if transaction['status'] == 'applying':
        from kinbot.run_recover import _apply_graph
        with (Path(graph)/'.manifest.lock').open('a') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX)
            _apply_graph(transaction['plan'], Path(graph).resolve(), Path(transaction['archive']))
            transaction['status'] = 'applied'
            atomic_json(marker, transaction)
    if transaction['status'] == 'applied' and transaction.get('request'):
        path = Path(transaction['request'])
        r = json.loads(path.read_text())
        if r['status'] == 'applying':
            control.acknowledge(path, r, 'applied', archive=transaction['archive'])


def check_commit(graph, spec, state, reset):
    """Check manifest compatibility before cancelling any affected work."""
    from kinbot.anl import dispatch
    dispatch.validate_spec(spec)
    original = json.loads((graph/'workflow.json').read_text())
    if original != spec:
        # Old pinned workers do not know the manifest read lock. They may
        # start at any moment, so defer manifest edits until they finish.
        for name,entry in state['tasks'].items():
            if entry['status'] in ('submitted','submitting') and name not in reset:
                record = json.loads((graph/'tasks'/name/'task.json').read_text())
                if record.get('live_control_schema') != 1:
                    raise RuntimeError('Resource update deferred until legacy workers in this graph finish.')


def commit(graph, spec, state, reset, reasons, request_path=None):
    check_commit(graph, spec, state, reset)
    archive = graph/'attempts'/('live_'+uuid.uuid4().hex)
    plan = dict(graph=str(graph), spec=spec, state=state, reset=sorted(reset), reasons=reasons,
                before={n: file_sha256(graph/n) for n in ('workflow.json', 'state.json')},
                complete=[n for n,e in state['tasks'].items() if e['status']=='complete'],
                active={n:e['job_id'] for n,e in state['tasks'].items()
                        if e['status']=='submitted' and n not in reset})
    atomic_json(graph/'.live_control_transaction.json', dict(
        status='applying', plan=plan, archive=str(archive),
        request=str(request_path) if request_path else None))
    resume(graph)


def _ownership(directory, job, *, record=None):
    r = slurm_job_record(str(job)) if record is None else record
    if (r.get('JobId') != str(job) or not r.get('JobState') or not r.get('WorkDir')
            or Path(r['WorkDir']).resolve() != directory.resolve()):
        raise RuntimeError('Scheduler ownership/state unavailable; no retry will be submitted.')
    if r.get('UserId'):
        from kinbot.qc_scratch import account_identity
        if r['UserId'].split('(')[0] != account_identity()['name']:
            raise RuntimeError('Allocation belongs to a different account.')
    return r


def _finished_graph_attempt(directory, job, *, allow_completed):
    """Prove a saved graph attempt ended even after slurmctld purged it.

    This never cancels a job, infers success, or substitutes for a failed
    queue query. Ordinary generated workers retain their existing rules.
    """
    from kinbot.anl import dispatch
    directory = Path(directory).resolve()
    if directory.parent.name != 'tasks':
        return False
    graph, spec, state = dispatch._load(directory.parent.parent)
    entry = state['tasks'].get(directory.name, {})
    if (str(entry.get('job_id')) != str(job)
            or entry.get('status') not in ('failed', 'unavailable', 'complete')):
        return False
    task = next(t for t in spec['tasks'] if t['id'] == directory.name)
    if not (directory/'execution.json').is_file():
        # SIGKILL at a walltime limit prevents the worker receipt. Reuse the
        # general job-specific timeout proof (including its saved ownership
        # and timestamp checks), never infer termination from missing files.
        from kinbot import timeout_retry
        dispatch._verify_stage_files(graph, task, entry)
        return (not dispatch._job_active(str(job)) and timeout_retry.evidence(
            directory, job, [directory/'slurm.stderr']) is not None)
    outcome = json.loads((directory/'execution.json').read_text())
    if (str(outcome.get('slurm_job_id')) != str(job)
            or outcome.get('status') not in ({'failed', 'executed'} if allow_completed else {'failed'})):
        return False
    dispatch._verify_stage_files(graph, task, entry)
    dispatch._verify_execution(graph, task, entry, outcome)
    return not dispatch._job_active(str(job))


def stop_attempt(directory, job, request_path, record, *, allow_completed=False):
    """An explicit stop is reconciled before any file or ledger can be reset."""
    from kinbot.run_refresh import _cancel_held_pending, _pending_request
    retained = slurm_job_record(str(job))
    if not retained and _finished_graph_attempt(directory, job, allow_completed=allow_completed):
        return True
    if not retained and not allow_completed and record.get('kind') == 'task':
        # This is restricted to an explicit request; never automatically
        # replace arbitrary vanished allocations or reinterpret native errors.
        from kinbot.anl import dispatch, lost_job
        if Path(directory).parent.name == 'tasks':
            _, spec, state = dispatch._load(Path(directory).parent.parent)
            task = next(t for t in spec['tasks'] if t['id'] == Path(directory).name)
            proof = lost_job.evidence(directory, job, task,
                                      state['tasks'].get(task['id'], {}))
            if proof is not None:
                record['lost_job_evidence'] = proof
                atomic_json(request_path, record)
                return True
    current = _ownership(directory, job, record=retained)
    phase = current['JobState'].split('+')[0]
    if phase == 'COMPLETED' and not allow_completed:
        control.acknowledge(request_path, record, 'retained', detail='Allocation completed before the request was applied; reconcile its result first.')
        return False
    if phase in SLURM_TERMINAL_STATES:
        return True
    if record['pending_only']:
        _cancel_held_pending([dict(job_id=str(job), directory=str(directory))], record, request_path)
        hold = record.get('pending_holds', {}).get(str(job), {})
        if hold.get('phase') == 'retained':
            control.acknowledge(request_path, record, 'retained', detail='Allocation started; pending-only request left it running.')
            return False
        return _ownership(directory, job)['JobState'] == 'CANCELLED'
    if not record['cancel_running']:
        raise ValueError('Allocation is live. Use pending-only or explicitly request cancel-running.')
    if not record.get('cancel_requested'):
        record['cancel_requested'] = True
        control.acknowledge(request_path, record, 'applying')
        record['cancel_diagnostic'] = _pending_request(['scancel', str(job)])
        atomic_json(request_path, record)
    current = _ownership(directory, job)
    if current['JobState'].split('+')[0] == 'COMPLETED' and not allow_completed:
        control.acknowledge(request_path,record,'retained',detail='Allocation completed during cancellation; reconcile its result first.')
        return False
    return current['JobState'].split('+')[0] in SLURM_TERMINAL_STATES


def process(graph):
    """Called with the graph drive lock, before ordinary reconciliation."""
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import descendants
    graph = Path(graph).resolve()
    from kinbot.pool_drain import pulse
    pulse(graph)
    resume(graph)
    from kinbot.anl import protocol
    # Pending previews bind the previous manifest. Apply those first, then
    # append the default QZ anharmonic branch without retiring other work.
    protocol.process(graph)
    if (graph/'workflow.json').exists() and not protocol.requests(graph):
        protocol.migrate(graph)
    if control.find(graph) is None:
        return
    for path, r in control.requests(graph/'tasks', kinds={'task', 'refinement'}):
        target = Path(r['target'])
        if target.parent != graph/'tasks':
            continue
        try:
            _, spec, state = dispatch._load(graph)
            ident = target.name
            current = control.generation(target, 'task')
            if any(current.get(k) != v for k,v in r['expected'].items() if k != 'status'):
                raise ValueError('Attempt changed after the request; refusing to affect a replacement generation.')
            if r['action'] == 'release':
                control.acknowledge(path, r, 'applied')
                continue
            entry = state['tasks'].get(ident, {})
            if entry.get('manual_rejection'):
                raise ValueError('Reaction is excluded; include its channel before requesting a task retry.')
            if r['kind'] == 'refinement':
                from kinbot.anl.refinement import process as refine
                refine(graph, path, r, spec, state)
                continue
            from kinbot.vrc_scan_recovery import rejected_point
            rejected = (spec.get('intent', {}).get('vrc_scan_distance_angstrom') is not None
                        and rejected_point(graph, ident))
            if entry.get('status') == 'submitting' or (entry.get('status') == 'complete' and not rejected):
                raise ValueError('Completed or ambiguous submissions cannot be retried through live controls.')
            unsent = (entry.get('status') == 'staged' and not entry.get('job_id')
                      and bool(r.get('fallback')))
            if entry.get('status') not in ('failed', 'unavailable', 'submitted') and not unsent and not rejected:
                raise ValueError('Only failed/unavailable or explicitly stopped allocations can be retried.')
            outcome = target/'execution.json'
            if outcome.exists() and json.loads(outcome.read_text()).get('status') == 'executed' and not rejected:
                raise ValueError('A native result exists; reconcile/reparse it before requesting a retry.')
            children = descendants(spec, ident)
            if any(state['tasks'].get(n, {}).get('status') in
                   ('complete', 'submitted', 'submitting', 'staged') for n in children):
                raise ValueError('Downstream work exists; changing its source requires a dependency-aware recovery, not a live retry.')
            # Validate everything before requesting any scheduler mutation.
            task = next(t for t in spec['tasks'] if t['id'] == ident)
            replacement = resource_spec(graph, spec, state, only={ident})
            if r.get('fallback'):
                from kinbot.manual_fallback import task_options
                replacement['tasks']=[task_options(t,r['fallback']) if t['id']==ident else t
                                      for t in replacement['tasks']]
            dispatch.validate_spec(replacement)
            if entry.get('job_id'):
                if not stop_attempt(target, entry['job_id'], path, r, allow_completed=rejected):
                    continue
            elif not unsent:
                from kinbot.anl.timeout_recovery import _job_id
                job = _job_id(target, ident, entry)
                if not job:
                    raise ValueError('Failed task has no verifiable allocation ID; reconcile it before retrying.')
                if not stop_attempt(target, job, path, r, allow_completed=rejected):
                    continue
            if outcome.exists() and json.loads(outcome.read_text()).get('status') == 'executed' and not rejected:
                control.acknowledge(path,r,'retained',detail='Native completion arrived during retry preparation; reconcile it first.')
                continue
            control.acknowledge(path, r, 'applying')
            reset = {ident} | {n for n in children if state['tasks'].get(n, {}).get('status') in ('blocked', 'skipped')}
            state.setdefault('manual_retries', {}).setdefault(ident, []).append(
                dict(request=r['id'], reason=r['reason'], previous=deepcopy(entry),
                     fallback=r.get('fallback'), requested_at=r['created_at']))
            commit(graph, replacement, state, reset, {n:r['reason'] for n in reset}, path)
        except (ValueError, RuntimeError, OSError) as error:
            # Once a transactional reset began, its journal must be resumed;
            # a scheduler uncertainty likewise remains pending for review.
            marker = graph/'.live_control_transaction.json'
            in_transaction = marker.exists() and json.loads(marker.read_text()).get('status') == 'applying'
            if in_transaction:
                raise
            control.acknowledge(path, r, 'blocked' if isinstance(error, ValueError) else r['status'],
                                error=str(error))
    _, spec, state = dispatch._load(graph)
    try:
        revised = resource_spec(graph, spec, state)
    except ValueError as error:
        atomic_json(graph/'.live_control_status.json', dict(
            status='blocked', error=str(error), updated_at=control.now()))
        return  # A bad resource edit must not terminate otherwise valid work.
    changed = {a['id'] for a,b in zip(spec['tasks'], revised['tasks']) if a != b}
    if changed:
        reset = {n for n in changed if n in state['tasks']}
        try:
            commit(graph, revised, state, reset, {n:'Live resource rule changed before submission.' for n in reset})
        except RuntimeError as error:
            if 'legacy workers' not in str(error):
                raise
            atomic_json(graph/'.live_control_status.json',dict(
                status='deferred',error=str(error),updated_at=control.now()))
            return
    status_file = graph/'.live_control_status.json'
    if status_file.exists() and json.loads(status_file.read_text()).get('status') != 'ready':
        atomic_json(status_file,dict(status='ready',updated_at=control.now()))


def resource_spec(graph, spec, state, only=None):
    """Preserve base resources so removing a rule restores unsubmitted work."""
    result = deepcopy(spec)
    baseline = graph/'.live_resource_baseline.json'
    saved = json.loads(baseline.read_text()) if baseline.exists() else {}
    changed_baseline = False
    for task in result['tasks']:
        ident = task['id']
        if state['tasks'].get(ident, {}).get('manual_rejection'):
            continue
        if only is not None and ident not in only:
            continue
        if only is None and state['tasks'].get(ident, {}).get('status') not in (None, 'staged', 'blocked', 'skipped'):
            continue
        options = control.overrides(graph/'tasks'/ident)
        if not options and ident not in saved:
            continue
        if ident not in saved:
            saved[ident] = deepcopy(task['resources'])
            changed_baseline = True
        task['resources'] = {**saved[ident], **options}
        manual_resources = task.get('manual_fallback', {}).get('options', {}).get('resources', {})
        task['resources'].update(manual_resources)
        # A live rule must not undo a longer allocation already justified by
        # timeout recovery on the next polling iteration.
        from kinbot.run_walltime import seconds
        for retry in state.get('timeout_retries',{}).get(ident,[]):
            if isinstance(retry,dict) and retry.get('walltime'):
                if seconds(task['resources']['walltime']) < seconds(retry['walltime']):
                    task['resources']['walltime'] = retry['walltime']
                if retry.get('partition'):
                    task['resources']['partition'] = retry['partition']
        if options or manual_resources:
            validate_native_resources(task,{**options, **manual_resources})
            task['resources']['live_override'] = True
            if 'memory_mb' in options or 'memory_mb' in manual_resources:
                task['resources']['use_all_node_memory'] = False
    from kinbot.anl.dispatch import validate_spec
    validate_spec(result)
    if changed_baseline:
        atomic_json(baseline, saved)
    return result


def validate_native_resources(task,options):
    """Custom fixed native decks cannot safely inherit allocation-only edits."""
    import re
    if task['kind'] != 'external':
        return
    backend = task['backend'].lower()
    text = task['input_template']
    command = ' '.join(task['command'])
    if backend == 'molpro':
        if 'cores' in options and '{cores}' not in command:
            raise ValueError('Molpro CPU edits require a {cores} command placeholder.')
        if set(options)-{'walltime'} and ('{molpro_stack_mw}' not in command
                or re.search(r'(?im)^\s*memory\s*,',text)):
            raise ValueError('Molpro resource edits require dynamic command-line stack memory and no fixed MEMORY card.')
    elif backend == 'gaussian':
        if 'cores' in options and '{{CORES}}' not in text:
            raise ValueError('Gaussian CPU edits require {{CORES}} in the native input.')
        if 'memory_mb' in options and not any(t in text for t in ('{{MEMORY_MB}}','{{WORK_MEMORY_MB}}')):
            raise ValueError('Gaussian memory edits require a native memory placeholder.')
    elif backend in ('mrcc','cfour'):
        if 'memory_mb' in options and not any(t in text for t in ('{{MEMORY_MB}}','{{WORK_MEMORY_MB}}')):
            raise ValueError(f'{backend} memory edits require a native memory placeholder.')
