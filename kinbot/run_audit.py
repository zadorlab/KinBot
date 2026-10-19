"""Read-only salvage, dependency, and resource-policy report for a prepared run.

No dispatcher advance/preflight, submission, cancellation, native executable,
or scratch probe is called. Optional reports must be outside the run tree.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path, PurePosixPath
import pwd
import os
import re
import shlex

from kinbot.execution import file_sha256, slurm_job_record
from kinbot.run_recap import read_json, scan


GRAPH_PATTERNS = (
    '[0-9]*/anl_automatic/states/*/composite',
    '[0-9]*/vrctst/molpro_scan/*', '[0-9]*/vrctst/molpro_pivots/*',
    '[0-9]*/vrctst/molpro/dispatch*',
)
ORDER_NOTE = (
    'Dependency levels are precedence, not a global queue or start-time prediction. '
    'ANL rotates among sorted species/TS graph keys and selects the first staged '
    'task in each persisted task map. Its cursor is not recorded. VRC runs '
    'concurrently; later scan/correction/sampling tasks are generated dynamically. '
    'Branch slots, completion/failure/retry decisions and Slurm determine actual starts.'
)


def hostnames(value):
    """Expand the numeric Slurm hostlist form without querying the scheduler."""
    if not value or value in ('(null)', 'N/A', 'None'):
        return set()
    result = set()
    for item in re.split(r',(?![^\[]*\])', value):
        match = re.fullmatch(r'([^\[\]]*)\[([0-9,\-]+)\](.*)', item)
        if match:
            for part in match[2].split(','):
                ends = part.split('-')
                if len(ends) > 2:
                    raise ValueError('Unsupported hostlist')
                lo, hi = int(ends[0]), int(ends[-1])
                if hi < lo or hi - lo > 10000:
                    raise ValueError('Unsupported hostlist range')
                width = len(ends[0]) if ends[0].startswith('0') else 0
                for number in range(lo, hi + 1):
                    result.update(hostnames(match[1] + str(number).zfill(width) + match[3]))
                    if len(result) > 10000:
                        raise ValueError('Hostlist exceeds audit limit')
        elif re.fullmatch(r'[\w.\-]+', item):
            result.add(item)
        else:
            raise ValueError('Unsupported hostlist: ' + value)
        if len(result) > 10000:
            raise ValueError('Hostlist exceeds audit limit')
    return result


def script_policy(path, expected, entry):
    if not path.is_file():
        return {'status': 'not_staged', 'detail': 'No submission script yet'}
    options = []
    for line in path.read_text().splitlines():
        if line.lstrip().startswith('#SBATCH '):
            options.extend(shlex.split(line.strip()[len('#SBATCH'):]))
        elif line.strip() and not line.lstrip().startswith('#'):
            break  # Slurm ignores directives after the first command.
    excluded = set()
    for index, option in enumerate(options):
        if option in ('--exclude', '-x') and index + 1 < len(options):
            excluded.update(hostnames(options[index + 1]))
        elif option.startswith('--exclude='):
            excluded.update(hostnames(option.split('=', 1)[1]))
    exclusive = '--exclusive' in options and not any(
        x in ('--share', '--oversubscribe', '-s') or x.startswith('--exclusive=')
        for x in options)
    problems = []
    if not exclusive:
        problems.append('missing/conflicting full-node exclusivity')
    if set(expected) - excluded:
        problems.append('missing exclusions: ' + ','.join(sorted(set(expected) - excluded)))
    if entry.get('job_sha256') and file_sha256(path) != entry['job_sha256']:
        problems.append('script differs from committed hash')
    return {'status': 'fail' if problems else 'pass', 'exclusive': exclusive,
            'excluded_nodes': sorted(excluded), 'problems': problems}


def scheduler_policy(job, detail, expected):
    """Use scheduler evidence independently of the on-disk submission script."""
    if not detail:
        return {'status': 'unknown', 'detail': 'scontrol record unavailable; job may have finished'}
    if str(detail.get('JobId')) != str(job['job_id']):
        return {'status': 'unknown', 'detail': 'scontrol job identity mismatch'}
    exclusive = detail.get('OverSubscribe', detail.get('Shared'))
    result = {'exclusive_field': exclusive, 'excluded_nodes': detail.get('ExcNodeList'),
              'nodes': detail.get('NodeList'), 'cpus': detail.get('NumCPUs'),
              'memory_per_node': detail.get('MinMemoryNode'), 'status': 'pass', 'problems': []}
    if exclusive not in ('NO', '0'):
        result['status'] = 'unknown' if exclusive is None else 'fail'
        result['problems'].append('exclusive allocation not confirmed')
    if detail.get('ExcNodeList') is None:
        if result['status'] != 'fail':
            result['status'] = 'unknown'
        result['problems'].append('scheduler exclusions unavailable')
    else:
        missing = set(expected) - hostnames(detail['ExcNodeList'])
        if missing:
            result['status'] = 'fail'
            result['problems'].append('missing exclusions: ' + ','.join(sorted(missing)))
    allocated = hostnames(detail.get('NodeList', ''))
    if allocated & set(expected):
        result['status'] = 'fail'
        result['problems'].append('allocation is on an excluded node')
    return result


def scratch_policy(directory, *, home, account, nodes=()):
    """Inspect compute-node provenance; never probe /scratch on the login host."""
    from kinbot.qc_scratch import LOCAL_FILESYSTEMS
    files = sorted(directory.glob('*.runtime.json'), key=lambda p: p.stat().st_mtime, reverse=True)
    if not files:
        return {'status': 'unobserved', 'detail': 'No native runtime record yet'}
    path = files[0]
    data = read_json(path)
    root = PurePosixPath(data.get('scratch_root', '/'))
    work = PurePosixPath(data.get('native_work_directory', '/'))
    scratch = PurePosixPath(data.get('scratch_directory', '/'))
    preferred = PurePosixPath('/scratch') / account
    problems, unknown = [], []
    if root == PurePosixPath(home) or PurePosixPath(home) in root.parents:
        problems.append('scratch is in HOME')
    if root == PurePosixPath(directory) or PurePosixPath(directory) in root.parents:
        problems.append('scratch is in shared results')
    if not data.get('scratch_io_checked'):
        unknown.append('local I/O probe not recorded')
    filesystem = data.get('scratch_mount', {}).get('filesystem')
    if filesystem is None:
        unknown.append('filesystem not recorded')
    elif filesystem not in LOCAL_FILESYSTEMS:
        problems.append('verified local filesystem not recorded')
    if any(x == PurePosixPath('/') for x in (root, work, scratch)):
        unknown.append('native work/scratch paths missing')
    elif root not in scratch.parents or scratch not in work.parents:
        problems.append('native work/scratch paths not within the selected root')
    host = str(data.get('scratch_host', '')).split('.')[0]
    if nodes and host not in {n.split('.')[0] for n in nodes}:
        unknown.append('runtime host differs from current allocation; possibly an earlier attempt')
    # A verified alternate local path is allowed, but must remain visible.
    return {'status': 'fail' if problems else ('unknown' if unknown else
                ('preferred_local' if root == preferred else 'alternate_local')),
            'file': str(path), 'host': data.get('scratch_host'),
            'scratch_root': str(root), 'native_work_directory': str(work),
            'source': data.get('scratch_source'), 'mount': data.get('scratch_mount'),
            'available_mb_at_start': data.get('scratch_available_mb'),
            'minimum_mb': data.get('scratch_minimum_mb'),
            'rejected_candidates': data.get('scratch_rejected_candidates', []),
            'problems': problems, 'unknown': unknown}


def pinned_policy(root, case):
    """Check policy source in the recorded driver snapshot, without importing it."""
    launch = case/'launch.json'
    if not launch.exists():
        return {'status': 'unknown', 'detail': 'No controller launch record'}
    record = read_json(launch)
    if not record.get('code'):
        return {'status': 'unknown', 'detail': 'Launch does not record pinned code'}
    pin = Path(record['code']).resolve()
    if not pin.is_relative_to(root):
        return {'status': 'unknown', 'detail': 'Pinned code is outside the run', 'code': str(pin)}
    manifest = read_json(pin/'manifest.json')
    files = ['qc_scratch.py', 'qc.py', 'anl/runtime.py', 'anl/dispatch.py',
             'ase_modules/calculators/molpro.py', 'ase_modules/calculators/gaussian.py']
    checked = {}
    for name in files:
        digest = file_sha256(pin/'kinbot'/name)
        checked[name] = (digest == manifest.get('kinbot/' + name) and
                         digest == file_sha256(Path(__file__).parent/name))
    return {'status': 'matches_audit_code' if all(checked.values()) else 'differs',
            'code': str(pin), 'files': checked,
            'scope': 'Recorded driver snapshot and these source files; not proof of live process environment'}


def task_rows(graph, historical, live, queue_error, expected, details, account, home):
    from kinbot.anl import dispatch
    graph, spec, state = dispatch._load(graph)
    entries = state['tasks']
    by_id = {t['id']: t for t in spec['tasks']}
    depths, visiting = {}, set()
    def depth(ident):
        if ident not in depths:
            if ident in visiting:
                raise ValueError('Conditional dependency cycle; no ordering inferred')
            visiting.add(ident)
            task = by_id[ident]
            deps = dispatch._task_dependencies(task, state) | set(task.get('fallback_on', []))
            depths[ident] = 1 + max((depth(dep) for dep in deps), default=-1)
            visiting.remove(ident)
        return depths[ident]
    first_staged = next((k for k, v in entries.items() if v['status'] == 'staged'), None)
    rows = []
    for index, task in enumerate(spec['tasks']):
        ident = task['id']; entry = entries.get(ident, {})
        directory = graph/'tasks'/ident
        owners = [j for j in live if Path(j['directory']).resolve() == directory]
        deps = dispatch._task_dependencies(task, state)
        needs = {dep: entries.get(dep, {}).get('status', 'unstarted') for dep in sorted(deps)
                 if entries.get(dep, {}).get('status') != 'complete'}
        triggers = {dep: entries.get(dep, {}).get('status', 'unstarted')
                    for dep in task.get('fallback_on', [])}
        status = entry.get('status', 'unstarted')
        if owners:
            conflict = (len(owners) != 1 or (entry.get('job_id') and
                        str(entry['job_id']) != owners[0]['job_id']))
            phase = 'ownership_conflict' if conflict else 'allocated'
        elif status in ('complete', 'skipped', 'failed', 'blocked', 'unavailable'):
            phase = status
        elif status in ('submitted', 'submitting'):
            phase = 'scheduler_unknown' if queue_error else 'awaiting_reconciliation'
        elif any(v == 'complete' for v in triggers.values()):
            phase = 'fallback_not_needed'
        elif triggers and not all(v in ('unavailable', 'skipped') for v in triggers.values()):
            phase = 'conditional_fallback'
        elif needs:
            phase = 'dependency_wait'
        else:
            phase = 'ready' if status == 'staged' else 'ready_to_stage'
        row = {'graph': str(graph), 'task': ident, 'molecule': spec['molecule'],
               'status': phase, 'recorded_status': status, 'spec_position': index + 1,
               'dependency_level': depth(ident), 'prerequisites': sorted(deps),
               'waiting_on': needs, 'fallback_on': triggers,
               'first_staged_in_graph': ident == first_staged,
               'salvaged_at_recovery': ident in historical.get('complete', []),
               'retained_job_at_recovery': historical.get('active', {}).get(ident),
               'reset_at_recovery': ident in historical.get('reset', []),
               'geometry_source': dispatch._geometry_source(task, state),
               'recorded_geometry_source': entry.get('geometry_from'),
               'resources': task.get('resources'), 'jobs': owners,
               'script_policy': script_policy(directory/'job.slurm', expected, entry)}
        if owners or phase in ('ready', 'ready_to_stage'):
            row['scratch_policy'] = scratch_policy(directory, home=home, account=account,
                nodes=hostnames(owners[0].get('nodes', '')) if len(owners) == 1 else ())
        if len(owners) == 1:
            row['scheduler_policy'] = scheduler_policy(owners[0], details.get(owners[0]['job_id'], {}), expected)
        rows.append(row)
    return rows


def audit(root, *, scheduler=None, details=None):
    from kinbot.run_monitor import queue
    from kinbot.run_cleanup import _processes
    root = Path(root).expanduser().resolve(strict=True)
    prepared = read_json(root/'prepared.json')
    journal = read_json(root/'outage_recovery.json') if (root/'outage_recovery.json').exists() else {}
    historical = {str(Path(g['graph']).resolve()): g for g in journal.get('graphs', [])}
    all_jobs, error = queue() if scheduler is None else scheduler
    jobs = [j for j in all_jobs if Path(j['directory']).resolve().is_relative_to(root)]
    if details is None:
        with ThreadPoolExecutor(max_workers=4) as pool:
            records = list(pool.map(lambda j: slurm_job_record(j['job_id'], timeout=8), jobs))
        details = {j['job_id']: record for j, record in zip(jobs, records)}
    report = {'time': datetime.now(timezone.utc).isoformat(), 'root': str(root),
              'recovery_phase': journal.get('phase'), 'recovery_time': journal.get('time'),
              'archive': journal.get('archive'), 'code': journal.get('code'),
              'scheduler_error': error, 'order_note': ORDER_NOTE, 'cases': [], 'errors': []}
    try:
        processes = _processes(root)
    except (OSError, RuntimeError) as exc:
        processes = []
        report['errors'].append(str(exc))
    account = pwd.getpwuid(os.getuid())
    for filename in prepared['inputs']:
        path = Path(filename).resolve(); case = path.parent
        if not path.is_relative_to(root):
            raise ValueError('Prepared input is outside this run')
        params = read_json(path)
        result = {'case': case.name, 'input': str(path), 'tasks': [], 'errors': [],
                  'policy_problems': [],
                  'configured_policy': {k: params.get(k) for k in (
                      'queue_exclusive', 'rotdpy_exclusive', 'slurm_exclude_nodes',
                      'queue_job_limit', 'qc_scratch_fallbacks', 'qc_scratch_min_mb')},
                  'branches': [], 'controllers': [p for p in processes
                      if Path(p['cwd']).is_relative_to(case)], 'worker_inputs': []}
        expected = set(params.get('slurm_exclude_nodes', []))
        for key in ('queue_exclusive', 'rotdpy_exclusive'):
            if params.get(key) != 1:
                result['policy_problems'].append(f'{path.name}: {key} is not explicitly enabled')
        try:
            result['pinned_policy'] = pinned_policy(root, case)
        except (OSError, ValueError) as exc:
            result['pinned_policy'] = {'status': 'unknown', 'detail': str(exc)}
        for worker in sorted(case.glob('[0-9]*')):
            if not worker.is_dir():
                continue
            node_state = worker/'.kinbot_slurm_nodes.json'
            worker_expected = expected | (set(read_json(node_state).get('excluded_nodes', [])) if node_state.exists() else set())
            worker_input = worker/(worker.name + '.json')
            if worker_input.exists():
                actual = read_json(worker_input)
                result['worker_inputs'].append({'file': str(worker_input),
                    'policy': {k: actual.get(k) for k in result['configured_policy']}})
                worker_expected.update(actual.get('slurm_exclude_nodes', []))
                for key in ('queue_exclusive', 'rotdpy_exclusive'):
                    if actual.get(key) == 0:
                        result['policy_problems'].append(f'{worker_input}: {key} explicitly disabled')
            branch = worker/'.kinbot_post_reaction_branches.json'
            if branch.exists():
                result['branches'].append({'worker': worker.name, **read_json(branch)})
            for pattern in GRAPH_PATTERNS:
                # Scan one worker at a time so its learned exclusions are respected.
                for graph in scan(worker, pattern.split('/', 1)[1], result['errors']):
                    try:
                        if not graph.resolve().is_relative_to(root):
                            raise ValueError('Graph resolves outside the run')
                        rows = task_rows(graph, historical.get(str(graph.resolve()), {}),
                            jobs, error, worker_expected, details, account.pw_name, account.pw_dir)
                        result['tasks'].extend(rows)
                    except (OSError, ValueError, RuntimeError, KeyError) as exc:
                        result['errors'].append(f'{graph}: {exc}')
        result['expected_exclusions'] = sorted(expected)
        known = {j['job_id'] for r in result['tasks'] for j in r['jobs']}
        result['other_allocations'] = [dict(j, scheduler_policy=scheduler_policy(j, details.get(j['job_id'], {}), expected))
            for j in jobs if Path(j['directory']).resolve().is_relative_to(case) and j['job_id'] not in known]
        report['cases'].append(result)
    observed = {r['graph'] for c in report['cases'] for r in c['tasks']}
    report['recovery_inventory'] = [
        {'graph': graph, 'completed': item.get('complete', []), 'retained': item.get('active', {}),
         'graph_seen_now': graph in observed} for graph, item in historical.items()]
    return report


def label(row, root):
    return str(Path(row['graph']).relative_to(root)) + '/' + row['task']


def brief_label(row):
    graph = Path(row['graph'])
    name = graph.parent.name if graph.name == 'composite' else graph.name
    if len(name) > 50:
        name = name[:24] + '...' + name[-20:]
    return name + '/' + row['task']


def render(report, *, full=False, limit=4):
    lines = [f"{report['time']} scheduler={'UNKNOWN' if report['scheduler_error'] else 'responding'}",
             f"Recovery={report['recovery_phase']} | {report['root']}"]
    for case in report['cases']:
        rows = case['tasks']; policy = case['configured_policy']
        counts = Counter(r['status'] for r in rows)
        lines += ['', case['case'] + ': ' + ', '.join(f'{k}={v}' for k, v in sorted(counts.items())),
                  '  Salvaged at recovery=' + str(sum(r['salvaged_at_recovery'] for r in rows)) +
                  '; new recorded complete=' + str(sum(r['status'] == 'complete' and not r['salvaged_at_recovery'] for r in rows)),
                  f"  Config: exclusive={policy['queue_exclusive']}; ROTD exclusive={policy['rotdpy_exclusive']}; cap={policy['queue_job_limit']}",
                  '  Exclusions: ' + ','.join(case['expected_exclusions'])]
        if case['branches']:
            lines.append('  Branches: ' + '; '.join(f"{k}={v.get('status')} slots={v.get('slots')}"
                for b in case['branches'] for k, v in b.get('branches', {}).items()))
        checks = [r for r in rows if r['status'] in ('allocated', 'ready', 'ready_to_stage')]
        lines.append('  Current script checks: ' + str(dict(Counter(r['script_policy']['status'] for r in checks))) +
                     '; scratch observations: ' + str(dict(Counter(r.get('scratch_policy', {}).get('status', 'unknown') for r in checks))))
        lines.append('  Live scheduler policy: ' + str(dict(Counter(
            r.get('scheduler_policy', {}).get('status', 'unknown') for r in rows if r['jobs']))) +
            '; pinned policy code=' + case['pinned_policy']['status'])
        failures = [r for r in checks if r['script_policy']['status'] == 'fail' or
                    r.get('scheduler_policy', {}).get('status') == 'fail' or
                    r.get('scratch_policy', {}).get('status') == 'fail']
        lines.append(f'  Task policy discrepancies={len(failures)}; input discrepancies={len(case["policy_problems"])}; read errors={len(case["errors"])}; local processes={len(case["controllers"])}')
        groups = [('LIVE', [r for r in rows if r['jobs']]),
                  ('READY (not promised launch order)', [r for r in rows if r['status'] in ('ready', 'ready_to_stage')])]
        if full:
            groups += [('SALVAGED', [r for r in rows if r['salvaged_at_recovery']]),
                       ('OTHER / DEPENDENCY ORDER', [r for r in rows if not r['jobs'] and not r['salvaged_at_recovery']
                           and r['status'] not in ('ready', 'ready_to_stage')])]
        for title, subset in groups:
            ordered = sorted(subset, key=lambda r: (r['dependency_level'], r['graph'], r['spec_position']))
            for row in ordered if full else ordered[:limit]:
                jobs = ','.join(j['job_id'] + ':' + j['state'] for j in row['jobs'])
                scratch = row.get('scratch_policy', {})
                name = label(row, report['root']) if full else brief_label(row)
                lines.append(f"  {title} {name} [{row['status']}] {jobs}" +
                    (f" scheduler_policy={row.get('scheduler_policy', {}).get('status')} scratch={scratch.get('status')} {scratch.get('scratch_root', '')}" +
                     (' RETAINED PRE-UPDATE JOB' if row['retained_job_at_recovery'] in {j['job_id'] for j in row['jobs']} else '') if jobs else ''))
                if full:
                    lines += [f"    molecule={row['molecule'].get('smiles', row['molecule'].get('symbols'))}; depth={row['dependency_level']}; resources={row['resources']}",
                              f"    geometry={row['geometry_source']}; waiting_on={row['waiting_on']}; fallback_on={row['fallback_on']}",
                              f"    script={row['script_policy']}; scheduler={row.get('scheduler_policy')}; scratch={scratch}"]
            if not full and len(ordered) > limit:
                lines.append(f'  {len(ordered)-limit} more {title} tasks in full report')
        if full:
            lines += ['  Other allocations: ' + json.dumps(case['other_allocations']),
                      '  Worker input policies: ' + json.dumps(case['worker_inputs']),
                      '  Pinned policy source: ' + json.dumps(case['pinned_policy'])]
        elif case['other_allocations']:
            lines.append(f'  Other live allocations={len(case["other_allocations"])} (outside materialized graphs; full report lists IDs/policies)')
        lines.extend('  INPUT POLICY: ' + problem for problem in case['policy_problems'])
        lines.extend('  READ ERROR: ' + e for e in case['errors'])
    if full:
        lines += ['', 'RECOVERY JOURNAL INVENTORY (verified at recovery, not rehashed by this audit)']
        for item in report['recovery_inventory']:
            lines.append(f"{item['graph']}: salvaged={item['completed']}; retained={item['retained']}; seen_now={item['graph_seen_now']}")
    lines += ['', ORDER_NOTE,
              'Salvaged counts come from the recovery journal; current completions are recorded, not freshly revalidated.',
              'Unstaged/pending tasks have no observed scratch yet. Retained pre-update jobs keep their old runtime.',
              'Only materialized task graphs are enumerated. Later VRC points, corrections, ROTD sampling and MESS are conditional future stages.',
              'No jobs, inputs, task states, logs or scratch directories were changed.']
    lines.extend('READ ERROR: ' + e for e in report['errors'])
    if report['scheduler_error']:
        lines.append('SCHEDULER UNKNOWN: ' + report['scheduler_error'])
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path, help='Existing report directory outside the run; files must not exist')
    parser.add_argument('--limit', type=int, default=3)
    args = parser.parse_args(argv)
    if args.limit < 0:
        parser.error('--limit must be nonnegative')
    if args.output:
        destination = args.output.resolve(strict=True)
        if not destination.is_dir() or destination.is_relative_to(args.root.resolve(strict=True)):
            parser.error('--output must be a directory outside the run')
        if any((destination/name).exists() for name in ('audit.txt', 'audit.json')):
            parser.error('Report files already exist; choose a new output directory')
    report = audit(args.root)
    if args.output:
        for name, text in [('audit.txt', render(report, full=True)),
                           ('audit.json', json.dumps(report, indent=2))]:
            with (destination/name).open('x') as stream:
                stream.write(text + '\n')
    print(render(report, limit=args.limit))
    if args.output:
        print(f'Full reports: {destination}/audit.txt and audit.json')


if __name__ == '__main__':
    main()
