"""Live operator controls for KinBot; no source edits or driver refresh required.

Use show to inspect receipts. Resource rules affect future attempts. Requests
are applied by their owning drivers; a stopped driver cannot consume requests.
"""
import argparse
import json
from pathlib import Path

from kinbot import live_control as control
from kinbot.execution import atomic_json


def inventory(root):
    """Use the monitor's bounded traversal, excluding archived attempts."""
    from kinbot.run_monitor import inventory as scan
    root = Path(root).resolve(strict=True)
    warnings = []
    rows = []
    for p in scan(root, warnings):
        if p.name != 'state.json' or not (p.parent/'workflow.json').is_file():
            continue
        state = json.loads(p.read_text())
        spec = json.loads((p.parent/'workflow.json').read_text())
        for task in spec['tasks']:
            entry = state['tasks'].get(task['id'], {})
            target = p.parent/'tasks'/task['id']
            rows.append(dict(target=target.relative_to(root).as_posix(), kind='task',
                             status=entry.get('status', 'waiting'), job_id=entry.get('job_id'),
                             resources=task['resources'], overrides=control.overrides(target)))
    for folder in ('qc','searches'):
        for p in (root/'.kinbot_controls'/folder).glob('*.json'):
            record = json.loads(p.read_text())
            target = Path(record['target'])
            if target.is_relative_to(root):
                rows.append({**record,'target':target.relative_to(root).as_posix(),
                             'overrides':control.overrides(target)})
    return {'targets': rows, 'warnings': warnings}


def pending_time(root, jobs, value, *, delta=False, apply=False):
    """Hold before modifying a pending limit, preserve pre-existing user holds."""
    from kinbot.run_walltime import jobs as owned_jobs, seconds, duration, partition_limits
    from kinbot.run_refresh import _pending_request
    root = Path(root).resolve(strict=True)
    live = {j['id']:j for j in owned_jobs(root)}
    selected = list(live) if jobs == ['all'] else jobs
    if not selected:
        raise ValueError('No jobs selected.')
    limits = partition_limits()
    plans = []
    for job in selected:
        if job not in live:
            raise ValueError(f'{job}: not an owned allocation under this root.')
        entry = live[job]
        if entry['state'] != 'PENDING':
            if jobs == ['all']:
                continue
            raise ValueError(f'{job}: only pending jobs can be changed by pending-time.')
        old = seconds(entry['limit'])
        requested = (old + (-1 if value.startswith('-') else 1)*seconds(value.lstrip('+-'))
                     if delta else seconds(value))
        if requested <= 0 or requested % 60 or requested > min(limits[x] for x in entry['partition'].split(',')):
            raise ValueError('Requested limit must be positive whole minutes and within partition MaxTime.')
        plans.append(dict(job_id=job, directory=entry['directory'], before=entry['limit'],
                          requested=duration(requested), status='reviewed'))
    if not apply:
        return plans
    directory = root/'.kinbot_controls/scheduler'
    import uuid
    for plan in plans:
        marker = directory/(uuid.uuid4().hex+'.json')
        plan['receipt'] = str(marker)
        atomic_json(marker, plan)
        from kinbot import run_policy
        policy = run_policy.find(root)
        from contextlib import nullcontext
        with run_policy._locked(policy) if policy else nullcontext():
            try:
                _pending_time_one(plan, marker, _pending_request)
            except (OSError,RuntimeError) as error:
                plan.update(status='unconfirmed',error=str(error))
                atomic_json(marker,plan)
    return plans


def _pending_time_one(plan, marker, command):
    from kinbot.run_walltime import seconds
    from kinbot.anl.live_control import _ownership
    job, directory = plan['job_id'], Path(plan['directory'])
    if plan['status'] in ('applied','retained'):
        return
    record = _ownership(directory, job)
    if record['JobState'] != 'PENDING' and 'original' not in plan:
        plan.update(status='retained', detail='Job started or finished before update.')
        atomic_json(marker, plan)
        return
    original = plan.setdefault('original', record)
    if not str(record.get('Priority', '')).isdecimal():
        raise RuntimeError('Cannot verify pending hold priority.')
    if seconds(record['TimeLimit']) != seconds(plan['before']) and plan['status']=='reviewed':
        raise RuntimeError('Time limit changed since review; no update performed.')
    try:
        if record['JobState'] == 'PENDING' and record['Priority'] != '0':
            plan['status'] = 'holding'; atomic_json(marker, plan)
            command(['scontrol', 'hold', job])
            record = _ownership(directory, job)
        if record['JobState'] != 'PENDING' or record.get('Priority') != '0':
            plan.update(status='retained', detail='Hold did not establish a pending job; no time update sent.')
            return
        plan['status'] = 'updating'; atomic_json(marker, plan)
        diagnostic = command(['scontrol', 'update', 'JobId='+job, 'TimeLimit='+plan['requested']])
        record = _ownership(directory, job)
        if seconds(record['TimeLimit']) != seconds(plan['requested']):
            plan.update(status='unconfirmed', error=diagnostic or 'Scheduler did not confirm the requested limit.')
        else:
            plan.update(status='applied', verified=record)
    except (RuntimeError,OSError) as error:
        plan.update(status='unconfirmed',error=str(error))
    finally:
        if original.get('Priority') != '0':
            # Only undo our hold. A running job can also retain requeue priority zero.
            try:
                current = _ownership(directory, job)
                if current.get('Priority') == '0' and current['JobState'] in ('PENDING','RUNNING','CONFIGURING','SUSPENDED'):
                    command(['scontrol', 'release', job])
                    current = _ownership(directory, job)
                    if current.get('Priority') == '0':
                        raise RuntimeError('Our hold still needs release; inspect scheduler receipt.')
            except (RuntimeError,OSError) as error:
                plan.update(status='release_unconfirmed',error=str(error))
        atomic_json(marker, plan)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    enable = sub.add_parser('enable'); enable.add_argument('root', type=Path)
    enable.add_argument('--hold-failures', action=argparse.BooleanOptionalAction, default=None)
    enable.add_argument('--hold-task-failures', action=argparse.BooleanOptionalAction, default=None)
    resubmit = sub.add_parser('resubmit', help='Validate an edited PES input; ask its existing owner to apply it.')
    resubmit.add_argument('input', type=Path)
    resubmit.add_argument('--run-directory', type=Path,
                          help='PES directory when its input is stored elsewhere.')
    resubmit.add_argument('--apply', action='store_true')
    resubmit.add_argument('--status', action='store_true')
    cap = sub.add_parser('capacity'); cap.add_argument('root', type=Path)
    cap.add_argument('--max-jobs', type=int)
    cap.add_argument('--exclude', action='append', default=[])
    cap.add_argument('--allow', action='append', default=[])
    cap.add_argument('--apply', action='store_true')
    resume = sub.add_parser('resume-scheduler'); resume.add_argument('receipt', type=Path)
    for name in ('show', 'list'):
        p = sub.add_parser(name); p.add_argument('root', type=Path)
    resource = sub.add_parser('resources'); resource.add_argument('root', type=Path)
    resource.add_argument('--name', required=True)
    resource.add_argument('--select', action='append', default=[])
    resource.add_argument('--cores', type=int); resource.add_argument('--memory-mb', type=int)
    resource.add_argument('--walltime'); resource.add_argument('--remove', action='store_true')
    for name in ('retry', 'release'):
        p = sub.add_parser(name); p.add_argument('root', type=Path)
        p.add_argument('--target', action='append', required=True,
                       help='Exact target path relative to ROOT; repeat for a batch.')
        p.add_argument('--kind', choices=['task','qc','reaction','conformer'], default='task')
        p.add_argument('--reason', required=True)
        p.add_argument('--fallback', type=Path,
                       help='JSON next-attempt native/optimizer/environment overrides; embedded in receipt.')
        g = p.add_mutually_exclusive_group()
        g.add_argument('--cancel-running', action='store_true')
        g.add_argument('--pending-only', action='store_true')
        p.add_argument('--apply', action='store_true', help='Publish the reviewed request; default only previews.')
    revisit = sub.add_parser('revisit')
    revisit.add_argument('root', type=Path)
    revisit.add_argument('--worker', required=True, help='Worker directory relative to ROOT.')
    revisit.add_argument('--vdw-detection', type=float)
    revisit.add_argument('--rejected', action='store_true')
    revisit.add_argument('--reload-vrc-models', action='store_true',
                         help='Read additional reviewed models from this worker input JSON.')
    revisit.add_argument('--reason', required=True)
    revisit.add_argument('--apply', action='store_true')
    wall = sub.add_parser('pending-time'); wall.add_argument('root', type=Path)
    group = wall.add_mutually_exclusive_group(required=True)
    group.add_argument('--job', action='append'); group.add_argument('--all-pending', action='store_true')
    group = wall.add_mutually_exclusive_group(required=True)
    group.add_argument('--walltime'); group.add_argument('--delta', help='Signed HH:MM:SS, e.g. --delta=-02:00:00')
    wall.add_argument('--apply', action='store_true')
    args = parser.parse_args(argv)
    if args.action == 'resubmit':
        from kinbot.input_refresh import submit, receipt, owner_directory
        if args.status:
            owner = owner_directory(args.input,directory=args.run_directory,require_live=False)
            if owner is None:
                raise ValueError('No PES input owner record was found.')
            result = receipt(owner)
        else:
            result = submit(args.input, apply=args.apply,directory=args.run_directory)
    elif args.action == 'resume-scheduler':
        from kinbot import run_policy
        from kinbot.run_refresh import _pending_request
        path = args.receipt.resolve(strict=True)
        plan = json.loads(path.read_text())
        policy = run_policy.find(Path(plan['directory']))
        from contextlib import nullcontext
        with run_policy._locked(policy) if policy else nullcontext():
            _pending_time_one(plan,path,_pending_request)
        result = plan
    elif args.action == 'capacity':
        from kinbot import run_policy, pool_drain
        root = args.root.resolve(strict=True)
        result = dict(current=run_policy.read(root), requested_max_jobs=args.max_jobs,
                      exclude=args.exclude,allow=args.allow)
        if args.apply:
            result = dict(policy=run_policy.update(root,max_jobs=args.max_jobs,
                exclude=[n for v in args.exclude for n in v.split(',')],
                allow=[n for v in args.allow for n in v.split(',')],manage_pending=True),
                drain=pool_drain.pulse(root,force=True))
    elif args.action == 'enable':
        result = control.configure(args.root, hold_failures=args.hold_failures,
                                   hold_task_failures=args.hold_task_failures)
    elif args.action == 'show':
        result = control.status(args.root)
    elif args.action == 'list':
        result = inventory(args.root)
    elif args.action == 'resources':
        values = {k:getattr(args,k) for k in ('cores','memory_mb','walltime') if getattr(args,k) is not None}
        result = control.resources(args.root, args.name, args.select, values, remove=args.remove)
    elif args.action == 'revisit':
        root = args.root.resolve(strict=True)
        target = (root/args.worker).resolve(strict=True)
        if not target.is_relative_to(root):
            raise ValueError('Worker must be within ROOT.')
        from kinbot.network_revisit import settings
        settings(args.vdw_detection,args.rejected)
        options = dict(vdW_detection=args.vdw_detection,rejected=args.rejected)
        if args.reload_vrc_models:
            saved=json.loads((target/(target.name+'.json')).read_text())
            if not saved.get('vrc_tst_surface_models'):
                raise ValueError('Worker input must contain reviewed vrc_tst_surface_models.')
            options['vrc_tst_surface_models']=saved['vrc_tst_surface_models']
        result = dict(target=str(target), options=options, expected=control.generation(target,'network'))
        if args.apply:
            result = control.request(root,target,kind='network',reason=args.reason,
                                     expected=result['expected'], options=options)
    elif args.action == 'pending-time':
        result = pending_time(args.root, ['all'] if args.all_pending else args.job,
                              args.delta or args.walltime, delta=args.delta is not None, apply=args.apply)
    else:
        root = args.root.resolve(strict=True)
        targets = [(root/target).resolve() for target in args.target]
        if any(not target.is_relative_to(root) for target in targets) or len(targets)!=len(set(targets)):
            raise ValueError('Targets must be distinct and inside ROOT.')
        expected = [control.generation(target, args.kind) for target in targets]
        fallback=None
        if args.fallback:
            from kinbot.manual_fallback import validate, task_options, worker_options
            fallback=validate(json.loads(args.fallback.read_text()),args.kind)
            if args.action!='retry' or args.kind not in ('task','qc','conformer','reaction'):
                raise ValueError('Fallbacks apply to task/qc/conformer/reaction retry targets.')
            # Validate every member before publishing any request in this batch.
            for target in targets:
                if args.kind=='task':
                    from kinbot.anl.dispatch import _load, validate_spec
                    _,spec,_=_load(target.parent.parent)
                    spec['tasks']=[task_options(t,fallback) if t['id']==target.name else t for t in spec['tasks']]
                    validate_spec(spec)
                elif args.kind=='reaction':
                    from kinbot.manual_fallback import reaction_workers
                    for worker in reaction_workers(target): worker_options(worker.read_text(),fallback)
                else:
                    from kinbot.qc_live import retry_sources
                    worker,_,_=retry_sources(target,expected[targets.index(target)].get('job_id'))
                    worker_options(worker.read_text(),fallback)
        result = [dict(target=str(t), kind=args.kind, action=args.action, expected=e) for t,e in zip(targets,expected)]
        if fallback is not None:
            for row in result: row['fallback']=fallback
        if args.apply:
            result = [control.request(root,t,kind=args.kind,action=args.action,
                                      expected=e,reason=args.reason,cancel_running=args.cancel_running,
                                      pending_only=args.pending_only,fallback=fallback) for t,e in zip(targets,expected)]
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
