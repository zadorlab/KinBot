"""Recover a stopped prepared PES run after an outage, preserving owned jobs.

Default: verify artifacts and print a recovery plan. --apply journals and
archives obsolete attempts. --launch additionally starts replacement drivers.
This tool never cancels jobs or changes electronic-structure parameters.
"""
import argparse
from collections import Counter
from contextlib import ExitStack
from copy import deepcopy
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import socket

from kinbot import run_cleanup, run_refresh, run_salvage
from kinbot.anl import dispatch
from kinbot.anl.salvage import descendants, geometry_invalidations
from kinbot.anl.validation import require_primary_geometries
from kinbot.execution import atomic_json, file_sha256


def _stopped(root, prepared=None):
    # This also works while input migration is interrupted between writing
    # an input and its new prepared.json hash. No PID signalling is needed.
    prepared = prepared or json.loads((root/'prepared.json').read_text())
    for filename in prepared['inputs']:
        case = Path(filename).resolve().parent
        if not case.is_relative_to(root):
            raise RuntimeError('Prepared input is outside the recovery root.')
        for directory in (case, *sorted(case.glob('[0-9]*'))):
            if not directory.is_dir():
                continue
            for name in ('launch.json', 'pes_status.json', 'kinbot_status.json'):
                path = directory/name
                if path.exists():
                    host = json.loads(path.read_text()).get('host')
                    if host and host.split('.')[0] != socket.gethostname().split('.')[0]:
                        raise RuntimeError(f'Run recovery on the original controller host {host}.')
    if run_cleanup._processes(root):
        raise RuntimeError('Controllers are still alive; outage recovery leaves them untouched.')


def _verify_result(graph, task, entry, outcome):
    dispatch._verify_stage_files(graph, task, entry)
    final_hash = dispatch._verify_execution(graph, task, entry, outcome)
    if task['kind'] == 'external':
        dispatch._validate_external_outputs(graph/'tasks'/task['id'], task)
        if task.get('result_parser'):
            from kinbot.anl.results import parse_result
            parser = task['result_parser']
            parse_result((graph/'tasks'/task['id']/parser['file']).read_text(errors='replace'), parser)
    return final_hash


def graph_plan(graph, jobs):
    """Reconcile in memory. Failed/missing results never become successes."""
    graph, spec, state = dispatch._load(graph)
    before = {name: file_sha256(graph/name) for name in ('workflow.json', 'state.json')}
    revised, updated = deepcopy(spec), deepcopy(state)
    originals = state['tasks']
    tasks = {task['id']: task for task in spec['tasks']}
    active, reset, reasons, verified = {}, set(), {}, set()
    for directory in (graph/'tasks').iterdir():
        if directory.is_symlink():
            raise RuntimeError(f'{directory}: task directory is a symlink; review before recovery.')
        if directory.is_dir() and directory.name not in tasks:
            raise RuntimeError(f'{directory}: unrecognized task directory; review before recovery.')
    for ident, task in tasks.items():
        entry = updated['tasks'].get(ident)
        directory = graph/'tasks'/ident
        owners = [j for j in jobs if Path(j['directory']).resolve() == directory]
        if len(owners) > 1:
            raise RuntimeError(f'{directory}: multiple live allocations; nothing adopted.')
        if owners:
            if entry is None:
                raise RuntimeError(f'{directory}: live allocation lacks committed input provenance.')
            dispatch._verify_stage_files(graph, task, entry)
            owner = owners[0]
            if entry.get('job_id') and str(entry['job_id']) != owner['job_id']:
                raise RuntimeError(f'{directory}: live job differs from the recorded owner.')
            entry.update(status='submitted', job_id=owner['job_id'])
            active[ident] = owner['job_id']
            continue
        if entry and entry.get('job_id'):
            other = [j for j in jobs if j['job_id'] == str(entry['job_id'])]
            if other:
                raise RuntimeError(f'{directory}: recorded job belongs to another task directory.')
        outcome_file = directory/'execution.json'
        try:
            outcome = json.loads(outcome_file.read_text()) if outcome_file.exists() else None
        except (OSError, ValueError) as error:
            raise RuntimeError(f'{outcome_file}: cannot read outcome; restore storage before recovery.') from error
        if outcome and outcome.get('status') == 'executed':
            if entry is None:
                raise RuntimeError(f'{directory}: result lacks a committed task entry.')
            # A successful sbatch RPC whose response was lost can leave
            # "submitting" alongside a fully valid completed execution.
            try:
                final_hash = _verify_result(graph, task, entry, outcome)
            except (ValueError, RuntimeError) as error:
                reset.add(ident)
                reasons[ident] = f'Artifact/native validation failed: {error}'
                entry.update(status='failed')
            else:
                entry.update(status='complete', execution='executed')
                for key in ('error', 'unavailable_reason', 'missing_since', 'failure_phase'):
                    entry.pop(key, None)
                if final_hash is not None:
                    entry['final_geometry_sha256'] = final_hash
                verified.add(ident)
            continue
        status = (entry or {}).get('status')
        if outcome and outcome.get('status') == 'failed':
            if entry is None:
                raise RuntimeError(f'{directory}: failed output lacks task provenance.')
            dispatch._verify_execution(graph, task, entry, outcome)
            phase = outcome.get('failure_phase') or entry.get('failure_phase')
            if phase not in ('execution', 'scheduler'):
                raise RuntimeError(f'{directory}: non-execution failure needs explicit review: {outcome.get("error")}')
            reset.add(ident)
            reasons[ident] = outcome.get('error', 'Failed native execution')
            entry.update(status='failed')
        elif status in ('submitted', 'submitting', 'complete'):
            reset.add(ident)
            reasons[ident] = 'No live allocation and no completed execution record; retry interrupted attempt.'
            entry.update(status='failed')
        elif status in ('failed', 'unavailable'):
            if entry.get('failure_phase') == 'preflight':
                raise RuntimeError(f'{directory}: unresolved environment/preflight failure.')
            reset.add(ident)
            reasons[ident] = entry.get('error') or entry.get('unavailable_reason') or 'Failed attempt'
        elif status == 'staged' or (entry is None and directory.exists()):
            reset.add(ident)
            reasons[ident] = 'Restage unsubmitted work with current code and local-storage policy.'

    invalid = geometry_invalidations(spec, updated)
    # Corrupt completed results also invalidate their actual consumers.
    for ident in list(reset):
        if ident in verified or (originals.get(ident, {}).get('status') == 'complete'):
            invalid.update(descendants(spec, ident))
    if invalid & active.keys():
        raise RuntimeError(f'{graph}: incompatible live geometry tasks {sorted(invalid & active.keys())}; leave running for explicit review.')
    reset.update(invalid)
    for ident in invalid:
        reasons.setdefault(ident, 'Upstream geometry must be recomputed; old dependent result cannot be retained.')
    keep = verified - reset
    require_primary_geometries(revised, preserve_tasks=keep | active.keys())
    dispatch.validate_spec(revised, allow_obsolete=True)
    # Stage all unstarted tasks afresh, including previously blocked/skipped
    # descendants. Preserve deliberate fallback skips with no failed source.
    while True:
        released = {ident for ident, entry in updated['tasks'].items()
                    if entry.get('status') == 'blocked' or
                    (entry.get('status') == 'skipped' and
                     set(entry.get('skipped_by', ())) & reset)}
        if released <= reset:
            break
        reset.update(released)
        for ident in released:
            reasons.setdefault(ident, 'Release the dependency after its source retry.')
    return {'graph': str(graph), 'before': before, 'spec': revised, 'state': updated,
            'reset': sorted(reset), 'reasons': reasons, 'complete': sorted(keep),
            'active': active, 'geometry_invalidations': sorted(invalid)}


def plan(root):
    root = Path(root).resolve(strict=True)
    _stopped(root)
    prepared, _ = run_salvage._cases(root)
    jobs = run_cleanup._jobs(root)  # Failed query raises, never treated as empty.
    graphs = []
    for path in run_refresh._graphs(root):
        if not path.resolve().is_relative_to(root):
            raise RuntimeError(f'{path}: graph points outside the recovery root.')
        print(f'Verifying {path.parent.relative_to(root)}', flush=True)
        graphs.append(graph_plan(path.parent, jobs))
    claimed = {j for graph in graphs for j in graph['active'].values()}
    if claimed != {j['job_id'] for j in jobs}:
        raise RuntimeError('Live job outside a verified ANL/VRC graph; no restart prepared.')
    # No controller can submit between these queries. A completed live job
    # may disappear; newly appearing or relocated jobs require a new plan.
    _stopped(root)
    current = run_cleanup._jobs(root)
    previous = {(j['job_id'], str(Path(j['directory']).resolve())) for j in jobs}
    if not {(j['job_id'], str(Path(j['directory']).resolve())) for j in current} <= previous:
        raise RuntimeError('Queue ownership changed during verification; repeat recovery.')
    return {'root': str(root), 'prepared': prepared, 'graphs': graphs, 'jobs': jobs,
            'time': datetime.now(timezone.utc).isoformat()}


def _apply_graph(item, root, archive):
    """Idempotent graph transaction. Old files always precede new ledgers."""
    graph = Path(item['graph'])
    saved = archive/'graphs'/graph.relative_to(root)
    done = saved/'committed.json'
    if done.exists():
        return
    saved.mkdir(parents=True, exist_ok=True)
    for name, expected in item['before'].items():
        backup = saved/name
        if not backup.exists():
            if file_sha256(graph/name) != expected:
                raise RuntimeError(f'{graph/name}: changed since recovery plan.')
            temporary = backup.with_suffix(backup.suffix+'.preparing')
            temporary.write_bytes((graph/name).read_bytes())
            if file_sha256(temporary) != expected:
                raise RuntimeError(f'{temporary}: recovery backup verification failed.')
            temporary.replace(backup)
        elif file_sha256(backup) != expected:
            raise RuntimeError(f'{backup}: recovery backup changed.')
    state = deepcopy(item['state'])
    pending = state.setdefault('pending_retries', {})
    for ident in item['reset']:
        source, target = graph/'tasks'/ident, saved/'tasks'/ident
        if source.exists():
            if target.exists():
                raise RuntimeError(f'{source}: task reappeared during recovery; refusing overwrite.')
            if source.is_symlink():
                raise RuntimeError(f'{source}: task directory is a symlink.')
            target.parent.mkdir(parents=True, exist_ok=True)
            source.rename(target)
        old = state['tasks'].pop(ident, {})
        was_started = old.get('status') not in (None, 'waiting', 'staged', 'blocked', 'skipped')
        lineage = {'attempt': old.get('attempt', 1) + int(was_started),
                   'outage_recovery_reason': item['reasons'].get(ident, 'Upstream recovery')}
        if target.exists():
            lineage['previous_attempt'] = str(target)
        pending[ident] = lineage
    atomic_json(graph/'workflow.json', item['spec'])
    state['spec_sha256'] = file_sha256(graph/'workflow.json')
    atomic_json(graph/'state.json', state)
    atomic_json(done, {'complete': item['complete'], 'active': item['active'], 'reset': item['reset']})


def _resume_inventory(root, record):
    if record.get('controller_host', socket.gethostname()) != socket.gethostname():
        raise RuntimeError('Resume recovery on its original controller host.')
    _stopped(root, record['prepared'])
    allowed = {(job['job_id'], str(Path(job['directory']).resolve())) for job in record['jobs']}
    live = run_cleanup._jobs(root)
    if not {(job['job_id'], str(Path(job['directory']).resolve())) for job in live} <= allowed:
        raise RuntimeError('New or relocated live allocation during recovery; refusing task changes.')


def recover(root, *, apply=False, launch=False):
    root = Path(root).resolve(strict=True)
    journal = root/'outage_recovery.json'
    record = json.loads(journal.read_text()) if journal.exists() else None
    if record and record['phase'] == 'launched':
        return record  # Never restart failed replacements in an unbounded loop.
    if not apply:
        if record:
            return record
        return {**plan(root), 'phase': 'reviewed'}
    if record is None:
        reviewed = plan(root)
        # Full environment checks before moving anything. Native scratch is
        # verified on each compute node by the native worker, not on login.
        run_refresh._preflight(root)
        _resume_inventory(root, reviewed)
        pin = run_salvage._pin_code(root, refresh=True)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive = root.with_name(root.name+'.outage_archive.'+stamp)
        archive.mkdir()
        record = {**reviewed, 'phase': 'preparing', 'archive': str(archive), 'code': str(pin),
                  'controller_host': socket.gethostname()}
        atomic_json(journal, record)
    if record['phase'] == 'preparing':
        _resume_inventory(root, record)
        archive = Path(record['archive'])
        for item in record['graphs']:
            print(f'Recovering {Path(item["graph"]).relative_to(root)}', flush=True)
            _apply_graph(item, root, archive)
        # Refresh the resource policy in the single user input; preserve all
        # chemistry, thresholds, indices, spin levels, and sampling settings.
        prepared = deepcopy(record['prepared'])
        for filename in prepared['inputs']:
            path = Path(filename).resolve()
            run_salvage._save(path, root, archive)
            values = json.loads((archive/path.relative_to(root)).read_text())
            values.update(queue_exclusive=1, rotdpy_exclusive=1)
            values.setdefault('qc_scratch_fallbacks', [])
            values.setdefault('qc_scratch_min_mb', 4096)
            atomic_json(path, values)
        run_salvage._save(root/'prepared.json', root, archive)
        for name in prepared['sha256']:
            prepared['sha256'][name] = file_sha256(root/name)
        atomic_json(root/'prepared.json', prepared)
        run_refresh._archive_controllers(root, archive)
        for filename in prepared['inputs']:
            for worker in Path(filename).parent.glob('[0-9]*'):
                if worker.is_dir():
                    for relative in ('.kinbot_post_reaction_branches.json',
                                     'anl_automatic/final_report.json', 'anl_automatic/execution.json',
                                     'anl_automatic/handoff.json'):
                        # The branch record was already copied above; keep
                        # its original evidence and remove only the stale copy.
                        path = worker/relative
                        if path.exists():
                            target = archive/path.relative_to(root)
                            if target.exists() and path.is_file() and file_sha256(target) == file_sha256(path):
                                path.unlink()
                            else:
                                run_salvage._save(path, root, archive, move=True)
        run_salvage._save(root/'salvage.json', root, archive)
        history = json.loads((root/'salvage.json').read_text())
        history.update(status='prepared', code=record['code'], inputs=prepared['inputs'])
        atomic_json(root/'salvage.json', history)
        record['phase'] = 'prepared'
        atomic_json(journal, record)
    if launch and record['phase'] == 'prepared':
        # launch() records each PID separately and never duplicates an
        # already launched case if the second launch is interrupted.
        run_salvage.launch(root)
        record['phase'] = 'launched'
        atomic_json(journal, record)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--launch', action='store_true')
    args = parser.parse_args(argv)
    if args.launch and not args.apply:
        parser.error('--launch requires --apply')
    root = args.root.resolve(strict=True)
    with ExitStack() as locks:
        if args.apply:
            handle = locks.enter_context(root.with_suffix('.salvage.lock').open('a'))
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with ExitStack() as worker_locks:
            journal = root/'outage_recovery.json'
            phase = json.loads(journal.read_text())['phase'] if journal.exists() else None
            # A prior invocation can have launched one case before interruption.
            # Its worker owns this lock legitimately; launch() adopts it below.
            if args.apply and phase not in ('prepared', 'launched'):
                for path in root.glob('*/[0-9]*/.kinbot_run.json'):
                    handle = worker_locks.enter_context(path.open())
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            result = recover(root, apply=args.apply)
        # New workers need these inode locks immediately on startup.
        if args.launch:
            result = recover(root, apply=True, launch=True)
    print(json.dumps({key: result[key] for key in ('phase', 'archive', 'code') if key in result}, indent=2))
    cases, retained = {}, []
    for graph in result['graphs']:
        case = Path(graph['graph']).relative_to(root).parts[0]
        cases.setdefault(case, Counter()).update(
            verified=len(graph['complete']), retained=len(graph['active']),
            restage=len(graph['reset']))
        retained.extend(graph['active'].values())
    for case, counts in sorted(cases.items()):
        print(case + ': ' + ', '.join(f'{key}={value}' for key, value in counts.items()))
    print('Retained allocation IDs: ' + (', '.join(sorted(retained)) or 'none'))
    if args.apply:
        print(f'Full recovery record: {root / "outage_recovery.json"}')


if __name__ == '__main__':
    main()
