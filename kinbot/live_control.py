"""Opt-in live controls; commands are data consumed by the owning controller.

No client edits a dispatcher ledger or an ASE database. Resource rules apply
only before submission; one-shot requests are fenced to an observed attempt.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import fnmatch
import json
import os
from pathlib import Path
import re
import threading
import uuid
import hashlib

from kinbot.execution import atomic_json, file_sha256

NAME = '.kinbot_live_control.json'
_LOCK = threading.RLock()


def now():
    return datetime.now(timezone.utc).isoformat()


def find(directory=None):
    directory = Path(directory or Path.cwd()).resolve()
    matches = [p/NAME for p in (directory, *directory.parents) if (p/NAME).is_file()]
    if len(matches) > 1:
        raise ValueError('Nested live-control roots are not supported.')
    return matches[0] if matches else None


def read(directory=None):
    path = find(directory)
    if path is None:
        return None
    value = json.loads(path.read_text())
    if (value.get('schema') != 1 or type(value.get('hold_failures')) is not bool
            or type(value.get('hold_task_failures', False)) is not bool):
        raise ValueError(f'Invalid live control: {path}')
    return value


@contextmanager
def locked(root):
    with _LOCK, (Path(root)/'.kinbot_live_control.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def configure(root, *, hold_failures=None, hold_task_failures=None):
    root = Path(root).resolve(strict=True)
    existing = find(root)
    if existing and existing.parent != root:
        raise ValueError(f'Use the existing live-control root: {existing.parent}')
    if not existing and any(root.rglob(NAME)):
        raise ValueError('Child live controls already exist; do not split ownership.')
    if hold_failures is not None and type(hold_failures) is not bool:
        raise ValueError('hold_failures must be boolean.')
    if hold_task_failures is not None and type(hold_task_failures) is not bool:
        raise ValueError('hold_task_failures must be boolean.')
    with locked(root):
        value = read(root) or dict(schema=1, revision=0, hold_failures=False, rules=[])
        if hold_failures is not None:
            value['hold_failures'] = hold_failures
        if hold_task_failures is not None:
            value['hold_task_failures'] = hold_task_failures
        return _save(root, value)


def _save(root, value):
    value.update(revision=value['revision']+1, updated_at=now())
    atomic_json(root/'.kinbot_controls/history'/f'{value["revision"]:06d}.json', value)
    atomic_json(root/NAME, value)
    return value


def validate_resources(values):
    if not isinstance(values, dict) or set(values)-{'cores', 'memory_mb', 'walltime'}:
        raise ValueError('Resources accept cores, memory_mb, walltime only.')
    for key in ('cores', 'memory_mb'):
        if key in values and (type(values[key]) is not int or values[key] <= 0):
            raise ValueError(f'{key} must be a positive integer.')
    if 'walltime' in values:
        from kinbot.run_walltime import validate_minimum
        validate_minimum(values['walltime'])


def validate_selector(selector):
    if (not isinstance(selector, str) or not selector or selector.startswith('/')
            or '..' in Path(selector).parts or '\x00' in selector):
        raise ValueError('Selectors must be paths/globs relative to the control root, without ..')


def resources(root, name, selectors=(), values=None, *, remove=False):
    root = Path(root).resolve(strict=True)
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', name):
        raise ValueError('Rule names accept letters, digits, dots, underscores and hyphens.')
    for selector in selectors:
        validate_selector(selector)
    values = values or {}
    validate_resources(values)
    if not remove and (not selectors or not values):
        raise ValueError('Specify selectors and at least one resource setting.')
    if find(root) != root/NAME:
        raise ValueError('Enable live controls at this root first.')
    with locked(root):
        config = read(root)
        config['rules'] = [r for r in config['rules'] if r['name'] != name]
        if not remove:
            config['rules'].append(dict(name=name, selectors=list(selectors), resources=values))
        return _save(root, config)


def overrides(target):
    path = find(Path(target).parent)
    if path is None:
        return {}
    relative = Path(target).resolve().relative_to(path.parent).as_posix()
    result = {}
    for rule in read(path.parent)['rules']:
        if any(fnmatch.fnmatchcase(relative, s) for s in rule['selectors']):
            result.update(rule['resources'])
    validate_resources(result)
    return result


@contextmanager
def _request_lock(root, target, kind):
    from contextlib import nullcontext
    from kinbot import pes_phase
    directory = target if kind=='network' else target.parent
    parent = directory.parent
    with locked(root):
        with pes_phase.locked(parent) if (parent/pes_phase.NAME).exists() else nullcontext():
            yield


def request(root, target, *, kind='task', action='retry', cancel_running=False,
            pending_only=False, reason='', expected=None, options=None, fallback=None):
    """Publish one explicit request, never submit/cancel from this function."""
    root = Path(root).resolve(strict=True)
    if find(root) != root/NAME:
        raise ValueError('Enable live controls at this root first.')
    target = Path(target)
    target = target.resolve() if target.is_absolute() else (root/target).resolve()
    if not target.is_relative_to(root) or target == root:
        raise ValueError('Target must be inside the live-control root.')
    if kind not in ('task', 'refinement', 'qc', 'reaction', 'conformer', 'network') or action not in ('retry', 'release'):
        raise ValueError('Unknown live-control target/action.')
    if pending_only and cancel_running:
        raise ValueError('Pending-only and cancel-running are mutually exclusive.')
    if not reason.strip():
        raise ValueError('Give a reason for the operator request.')
    if fallback is not None:
        from kinbot.manual_fallback import validate
        if action != 'retry' or kind not in ('task','qc','conformer','reaction'):
            raise ValueError('Per-attempt fallbacks require task, qc, conformer or reaction retry targets.')
        fallback=validate(fallback,kind)
    if kind == 'task' and options is not None:
        raise ValueError('Use a refinement request for geometry/frequency options.')
    if kind == 'refinement':
        from kinbot.anl.refinement import plan
        if (not isinstance(options, dict) or set(options) != {'refinement_plan'} or action != 'retry'
                or fallback is not None or pending_only):
            raise ValueError('Task refinement cannot be combined with other retry options.')
        reviewed = options['refinement_plan']
        current = plan(target.parent.parent, target.name, reviewed['options'])
        current.pop('replacement')
        if current != reviewed:
            raise ValueError('Refinement scope changed before publication; preview it again.')
    if kind == 'network':
        from kinbot.network_revisit import settings
        if not isinstance(options, dict) or set(options)-{'vdW_detection','rejected','vrc_tst_surface_models'} or action!='retry':
            raise ValueError('Network revisit accepts vdW_detection, rejected and reviewed VRC models.')
        settings(options.get('vdW_detection'),options.get('rejected',False))
        if 'vrc_tst_surface_models' in options:
            from kinbot.vrc_multistate import normalize_inline_model_paths, validate_surface_model
            normalize_inline_model_paths(options,target)
            for model in options['vrc_tst_surface_models'].values():
                validate_surface_model(model)
    if expected is None:
        expected = generation(target, kind)
    with _request_lock(root, target, kind):
        if kind == 'refinement':
            from kinbot.anl.refinement import require_open
            require_open(target.parent.parent)
        if kind in ('reaction','conformer','network'):
            generation(target,kind)  # Search closure and request publication share this lock.
        directory = root/'.kinbot_controls/requests'
        for p in directory.glob('*.json'):
            old = json.loads(p.read_text())
            if (old['target'] == str(target) and old['kind'] == kind
                    and old['status'] in ('queued', 'applying')):
                raise ValueError(f'A request for this target is already pending: {old["id"]}')
        ident = uuid.uuid4().hex
        record = dict(schema=1, id=ident, kind=kind, target=str(target), action=action,
                      expected=expected, cancel_running=cancel_running, pending_only=pending_only,
                      reason=reason, status='queued', created_at=now())
        if options is not None:
            record['options'] = options
        if fallback is not None:
            record['fallback'] = fallback
        atomic_json(directory/(ident+'.json'), record)
    return record


def generation(target, kind):
    target = Path(target).resolve()
    if kind == 'network':
        from kinbot.pes_phase import active_worker
        if not active_worker(target):
            raise ValueError('Network revisit requires a live revision-aware PES worker.')
        return {'input_sha256':file_sha256(target/(target.name+'.json'))}
    if kind in ('task', 'refinement'):
        from kinbot.anl.dispatch import _load
        graph = target.parent.parent
        if target.parent.name != 'tasks':
            raise ValueError('Dispatcher target must be GRAPH/tasks/TASK_ID.')
        _, spec, state = _load(graph)
        if target.name not in {t['id'] for t in spec['tasks']}:
            raise ValueError('Task does not exist in the workflow.')
        entry = state['tasks'].get(target.name, {})
        return {k: entry.get(k) for k in ('attempt', 'job_id', 'task_sha256', 'status')}
    if kind in ('reaction','conformer'):
        marker = scope_file(target,kind)
        if marker is None or not marker.exists():
            raise ValueError('This search has not registered a live retry point. Enable controls before the search starts.')
        scope = json.loads(marker.read_text())
        if scope['phase'] != 'active':
            from kinbot.pes_phase import active_worker
            if kind != 'reaction' or not active_worker(target.parent):
                raise ValueError('This search already handed off its results and has no live PES owner that can reopen exploration.')
    worker = Path(str(target)+'.py')
    if kind == 'qc' and not worker.is_file():
        raise ValueError(f'Original worker is unavailable: {worker}')
    result = {'worker_sha256': file_sha256(worker)} if worker.is_file() else {}
    attempt = qc_record(target)
    if attempt and attempt.is_file():
        result['job_id'] = json.loads(attempt.read_text())['job_id']
    marker = Path(str(target)+'.exitcode')
    if marker.is_file():
        result['exitcode'] = marker.read_text().strip()
    return result


def qc_record(target):
    marker = find(Path(target).parent)
    if marker is None:
        return None
    token = hashlib.sha256(str(Path(target).resolve()).encode()).hexdigest()[:24]
    return marker.parent/'.kinbot_controls/qc'/(token+'.json')


def observe_qc(target, job_id, status):
    path = qc_record(target)
    if path is None:
        return
    record = dict(target=str(Path(target).resolve()),kind='qc',job_id=str(job_id),status=status)
    previous = json.loads(path.read_text()) if path.exists() else None
    if previous != record:
        atomic_json(path,record)


def scope_file(target, kind):
    marker = find(Path(target).parent)
    if marker is None:
        return None
    token = hashlib.sha256((kind+str(Path(target).resolve())).encode()).hexdigest()[:24]
    return marker.parent/'.kinbot_controls/searches'/(token+'.json')


def search_scope(target,kind,status,*,close=False):
    """Close the retry window before handing an accepted network downstream."""
    path = scope_file(target,kind)
    if path is None:
        return True
    root = find(Path(target).parent).parent
    with locked(root):
        if close and any(r['target']==str(Path(target).resolve())
                         for _,r in requests(Path(target).parent,kinds={kind})):
            return False
        record = dict(target=str(Path(target).resolve()),kind=kind,status=status,
                      phase='closed' if close else 'active',pid=os.getpid())
        previous = json.loads(path.read_text()) if path.exists() else None
        if record != previous:
            atomic_json(path,record)
    return True


def requests(directory, *, kinds=None):
    marker = find(directory)
    if marker is None:
        return []
    directory = Path(directory).resolve()
    results = []
    for path in sorted((marker.parent/'.kinbot_controls/requests').glob('*.json')):
        value = json.loads(path.read_text())
        batch=value.get('input_revision')
        if batch:
            publication=json.loads(Path(batch['path']).read_text())
            if publication.get('id')!=batch['id'] or publication.get('status')!='published':
                continue  # A partially published input batch is never visible.
        if (value['status'] in ('queued', 'applying')
                and Path(value['target']).is_relative_to(directory)
                and (kinds is None or value['kind'] in kinds)):
            results.append((path, value))
    return results


def acknowledge(path, record, status, **details):
    record.update(status=status, updated_at=now(), **details)
    atomic_json(path, record)


def holding(target, *, kind='task'):
    """An opt-in failed-stage hold consumes no allocation and is releasable."""
    config = read(Path(target).parent)
    if not config or not (config['hold_failures'] or (kind=='task' and config.get('hold_task_failures',False))):
        return False
    marker = find(Path(target).parent)
    for p in (marker.parent/'.kinbot_controls/requests').glob('*.json'):
        r = json.loads(p.read_text())
        if (r['target'] == str(Path(target).resolve()) and r['kind'] == kind
                and r['action'] == 'release' and r['status'] == 'applied'):
            current = generation(target, kind)
            if all(current.get(k) == v for k,v in r['expected'].items() if k not in ('status','exitcode')):
                return False
    return True


def status(root):
    marker = find(root)
    if marker is None:
        return {'enabled': False}
    return {'enabled': True, 'root': str(marker.parent), 'config': read(root),
            'requests': [json.loads(p.read_text()) for p in
                         sorted((marker.parent/'.kinbot_controls/requests').glob('*.json'))],
            'scheduler_updates': [json.loads(p.read_text()) for p in
                         sorted((marker.parent/'.kinbot_controls/scheduler').glob('*.json'))]}


def summary(root):
    from collections import Counter
    marker = find(root)
    if marker is None:
        return None
    result = dict(hold_failures=read(root)['hold_failures'],
                  hold_task_failures=read(root).get('hold_task_failures',False), revision=read(root)['revision'])
    rows = [json.loads(p.read_text()) for p in
            sorted((marker.parent/'.kinbot_controls/requests').glob('*.json'))]
    result['requests'] = dict(Counter(r['status'] for r in rows))
    drain = marker.parent/'.kinbot_policy_drain_poll.json'
    if drain.exists():
        d = json.loads(drain.read_text())
        result['drain'] = {k:d.get(k) for k in ('status','error')}
    return result


def validate_options(options):
    if (not isinstance(options,dict) or set(options)-{'hold_failures','hold_task_failures','resource_rules'}
            or any(k in options and type(options[k]) is not bool for k in ('hold_failures','hold_task_failures'))):
        raise ValueError('queue_live_control accepts hold_failures, hold_task_failures and resource_rules.')
    rules=options.get('resource_rules',[])
    if not isinstance(rules,list):
        raise ValueError('resource_rules must be a list.')
    names=set()
    for rule in rules:
        if (not isinstance(rule,dict) or set(rule)!={'name','selectors','resources'}
                or not isinstance(rule['name'],str) or not re.fullmatch(r'[A-Za-z0-9_.-]+',rule['name'])
                or rule['name'] in names or not isinstance(rule['selectors'],list)
                or not rule['selectors'] or not rule['resources']):
            raise ValueError('Each resource rule needs a unique name, selectors, and resources.')
        names.add(rule['name'])
        for selector in rule['selectors']:validate_selector(selector)
        validate_resources(rule['resources'])


def initialize(directory, options):
    validate_options(options)
    if options and find(directory) is None:
        configure(directory,hold_failures=options.get('hold_failures',False),
                  hold_task_failures=options.get('hold_task_failures',False))
        for rule in options.get('resource_rules',[]):
            resources(Path(directory).resolve(),rule['name'],rule['selectors'],rule['resources'])
