"""Read-only outage/restart inventory. Never submit, reconcile, reset or cancel.

Recorded completion is distinct from verified reusable output. One scheduler
query covers the entire run; unavailable storage/scheduler state stays unknown.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fnmatch import fnmatchcase
import json
from pathlib import Path
import re
import sys


def scan(root, pattern, errors):
    """Like a small glob, but report inaccessible directories instead of hiding them."""
    paths = [root]
    parts = pattern.split('/')
    for index, part in enumerate(parts):
        matches = []
        for path in paths:
            try:
                matches.extend(child for child in path.iterdir()
                               if fnmatchcase(child.name, part) and
                               (index == len(parts)-1 or child.is_dir()))
            except OSError as error:
                errors.append(f'{path}: {error}')
        paths = matches
    return sorted(paths)


def read_json(path):
    if path.stat().st_size > 32 * 1024**2:
        raise ValueError(f'JSON exceeds 32 MiB: {path}')
    return json.loads(path.read_text())


def diagnostic(directory, prefix=''):
    """Read bounded current-file tails, not archives, integral files or WFUs."""
    found = []
    files = [p for p in directory.iterdir() if p.name.startswith(prefix) and p.is_file() and
             p.suffix in ('.out', '.log', '.err', '.stderr') and not p.name.endswith('.live.log')]
    for path in sorted(files, key=lambda p: p.stat().st_mtime, reverse=True)[:8]:
        with path.open('rb') as stream:
            stream.seek(max(0, path.stat().st_size - 65536))
            lines = stream.read().decode(errors='replace').splitlines()
        hits = [s.strip() for s in lines if re.search(
            r'^\s*\?|Input/output error|No space left|No lower point|QPErr|'
            r'Inconsistency found|New:.*Old:|NO CONVERGENCE|Convergence criterion not met|'
            r'Error termination|Frequency evaluation failed|frequency check failed|'
            r'MemoryError|Traceback|(?:Error|Exception):', s, re.I)]
        if hits:
            found.append({'file': str(path), 'lines': list(dict.fromkeys(hits))[-8:]})
        elif any(line.strip() for line in lines):
            found.append({'file': str(path), 'lines': lines[-6:], 'tail_only': True})
    return found


def graph_recap(graph, live, scheduler_error, verify):
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import geometry_invalidations
    directory, spec, state = dispatch._load(graph)
    invalid = geometry_invalidations(spec, state)
    rows = []
    for task in spec['tasks']:
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        work = directory / 'tasks' / ident
        job = str(entry.get('job_id', ''))
        row = {'graph': str(graph), 'task': ident, 'recorded_status': entry.get('status', 'waiting'),
               'job_id': job or None, 'geometry_invalid': ident in invalid,
               'disposition': 'waiting', 'reason': ''}
        execution = None
        owners = [j for j in live.values() if Path(j['directory']).resolve() == work.resolve()]
        if job and job in live and Path(live[job]['directory']).resolve() != work.resolve():
            row.update(disposition='unknown', reason='Recorded job ID belongs to a different working directory.')
        elif owners:
            job = str(owners[0]['job_id'])
            row['job_id'] = job
            row.update(disposition='active', scheduler_state=live[job]['state'],
                       reason='Retain this allocation; do not submit a duplicate.')
            if ident in invalid:
                row['reason'] = 'Active allocation has incompatible geometry; review before adoption.'
        elif job and scheduler_error and entry.get('status') in ('submitted', 'submitting'):
            row.update(disposition='unknown', reason='Scheduler unavailable; job ownership is unresolved.')
        else:
            if (work / 'execution.json').is_file():
                execution = read_json(work / 'execution.json')
            if ident in invalid and entry.get('status') not in ('waiting', 'staged', None):
                row.update(disposition='invalid_geometry', reason='Retry L3 before recomputing its dependent single points.')
            elif execution and execution.get('status') == 'executed':
                if verify:
                    try:
                        dispatch._verify_stage_files(directory, task, entry)
                        final_hash = dispatch._verify_execution(directory, task, entry, execution)
                        if entry.get('final_geometry_sha256') and final_hash != entry['final_geometry_sha256']:
                            raise ValueError('Accepted final geometry differs from execution')
                        if task['kind'] == 'external':
                            dispatch._validate_external_outputs(work, task)
                            if task.get('result_parser'):
                                from kinbot.anl.results import parse_result
                                parser = task['result_parser']
                                parse_result((work / parser['file']).read_text(errors='replace'), parser)
                        row['disposition'] = 'reusable'
                    except Exception as error:
                        row.update(disposition='review_required', reason=str(error))
                else:
                    row.update(disposition='recorded_success', reason='Run --verify to check hashes, geometry and native results.')
            elif execution and execution.get('status') == 'failed':
                row.update(disposition='failed', reason=execution.get('error', 'Execution failed'))
            elif entry.get('status') in ('failed', 'unavailable'):
                row.update(disposition='failed', reason=entry.get('unavailable_reason') or entry.get('error') or 'Native outcome needs review')
            elif entry.get('status') in ('complete', 'submitted', 'submitting'):
                row.update(disposition='unknown', reason='No live job or verified execution record; disappearance is not success.')
            elif entry.get('status') == 'skipped':
                row.update(disposition='skipped', reason=entry.get('skip_reason', ''))
            elif entry.get('status') == 'blocked':
                row.update(disposition='blocked', reason=entry.get('error', 'Dependency unavailable'))
            elif entry.get('status') == 'staged':
                row['disposition'] = 'ready'
        if row['disposition'] in ('failed', 'unknown', 'invalid_geometry', 'review_required') and work.is_dir():
            row['native_diagnostics'] = diagnostic(work)
        if execution:
            row['runtime'] = {k: v for k, v in execution.get('details', {}).items()
                              if k.startswith('scratch_') or k == 'native_work_directory'}
        if row['disposition'] in ('failed', 'unknown', 'review_required') and work.is_dir():
            row['local_recovery_locations'] = [read_json(path) for path in sorted(work.glob('*.runtime.json'))]
        rows.append(row)
    return rows


def recap(root, *, verify=False, scheduler=None):
    from kinbot.run_monitor import queue
    from kinbot.run_cleanup import _processes
    root = Path(root).expanduser().resolve(strict=True)
    # An inaccessible directory must raise instead of looking like zero tasks.
    list(root.iterdir())
    jobs, scheduler_error = queue() if scheduler is None else scheduler
    live = {str(j['job_id']): j for j in jobs}
    report = {'schema': 1, 'root': str(root), 'time': datetime.now(timezone.utc).isoformat(),
              'verified': verify, 'scheduler_error': scheduler_error, 'cases': [], 'errors': []}
    prepared = root / 'prepared.json'
    cases = ([Path(p).resolve().parent for p in read_json(prepared)['inputs']]
             if prepared.is_file() else [root])
    try:
        processes = _processes(root)
    except (OSError, RuntimeError) as error:
        processes = []
        report['process_error'] = str(error)
    known_jobs = set()
    for case in cases:
        if not case.is_relative_to(root):
            raise ValueError('Prepared case is outside the requested run')
        list(case.iterdir())
        result = {'case': case.name, 'directory': str(case), 'tasks': [], 'branches': [],
                  'controllers': [], 'worker_failures': [], 'worker_successes': [], 'l2_models': [], 'errors': []}
        result['live_processes'] = [p for p in processes if Path(p['cwd']).is_relative_to(case)]
        provenance = case/'resume_provenance.json'
        historical_exits = read_json(provenance).get('archived_failure_exits', {}) if provenance.is_file() else {}
        for pattern in ('pes_status.json', '[0-9]*/kinbot_status.json'):
            for path in scan(case, pattern, result['errors']):
                try:
                    result['controllers'].append({'file': str(path), **read_json(path)})
                except (OSError, ValueError) as error:
                    result['errors'].append(f'{path}: {error}')
        for path in scan(case, '[0-9]*/.kinbot_post_reaction_branches.json', result['errors']):
            try:
                result['branches'].append({'file': str(path), **read_json(path)})
            except (OSError, ValueError) as error:
                result['errors'].append(f'{path}: {error}')
        patterns = ('[0-9]*/anl_automatic/states/*/composite/workflow.json',
                    '[0-9]*/vrctst/molpro_scan/*/workflow.json')
        for pattern in patterns:
            for path in scan(case, pattern, result['errors']):
                try:
                    rows = graph_recap(path.parent, live, scheduler_error, verify)
                    result['tasks'].extend(rows)
                    known_jobs.update(r['job_id'] for r in rows if r.get('job_id'))
                except (OSError, ValueError, RuntimeError, KeyError) as error:
                    result['errors'].append(f'{path}: {error}')
        # Ordinary Gaussian/reference/pivot workers are not composite graphs.
        for pattern in ('[0-9]*/*.exitcode', '[0-9]*/vrctst/*.exitcode',
                        '[0-9]*/anl_automatic/reference_l2/*/*.exitcode'):
            for path in scan(case, pattern, result['errors']):
                try:
                    from kinbot.qc import _worker_exit_status
                    from kinbot.execution import file_sha256
                    if historical_exits.get(str(path.relative_to(case))) == file_sha256(path):
                        continue
                    code, explicit = _worker_exit_status(path)
                    if code != 0:
                        result['worker_failures'].append({'file': str(path),
                            'exitcode': code, 'generation_explicit': explicit,
                            'native_diagnostics': diagnostic(path.parent, path.stem)})
                    else:
                        result['worker_successes'].append({'file': str(path),
                            'disposition': 'recorded_success', 'generation_explicit': explicit})
                except (OSError, ValueError) as error:
                    result['errors'].append(f'{path}: {error}')
        for pattern in ('[0-9]*/thermochemistry_model.json',
                        '[0-9]*/anl_automatic/reference_l2/*/thermochemistry_model.json'):
            for path in scan(case, pattern, result['errors']):
                try:
                    model = read_json(path)
                    result['l2_models'].append({'file': str(path),
                        'species': list(model.get('species', {})),
                        'disposition': 'recorded_model'})
                except (OSError, ValueError) as error:
                    result['errors'].append(f'{path}: {error}')
        result['counts'] = dict(Counter(r['disposition'] for r in result['tasks']))
        report['cases'].append(result)
    owned = []
    for job in jobs:
        try:
            if Path(job['directory']).resolve().is_relative_to(root):
                owned.append({**job, 'matched_graph': str(job['job_id']) in known_jobs})
        except OSError as error:
            report['errors'].append(str(error))
    report['live_jobs'] = owned
    report['restart_ready'] = False  # Inventory never authorizes mutation or scientific acceptance.
    return report


def render(report, limit=12):
    lines = [f"{report['time']}  scheduler={'UNKNOWN' if report['scheduler_error'] else 'responding'}",
             f"Run: {report['root']}", f"Live allocations: {len(report['live_jobs'])}"]
    if report['scheduler_error']:
        lines.append('SCHEDULER: ' + report['scheduler_error'])
    if report.get('process_error'):
        lines.append('PROCESS CHECK: ' + report['process_error'])
    for case in report['cases']:
        lines.append(f"\n{case['case']}: " + ', '.join(f'{k}={v}' for k, v in case['counts'].items()))
        lines.append(f"  live local processes={len(case['live_processes'])}; worker failures={len(case['worker_failures'])}")
        lines.append(f"  recorded worker successes={len(case['worker_successes'])}; accepted L2 model files={len(case['l2_models'])}")
        for branch in case['branches']:
            for name, entry in branch.get('branches', {}).items():
                lines.append(f"  {name}: {entry.get('status')} {str(entry.get('error', '')).splitlines()[0][:180] if entry.get('error') else ''}")
        issues = [r for r in case['tasks'] if r['disposition'] in
                  ('failed', 'unknown', 'invalid_geometry', 'review_required')]
        for row in issues[:limit]:
            native = row.get('native_diagnostics', [])
            reason = ' | '.join(native[0]['lines']) if native else row['reason']
            lines.append(f"  {Path(row['graph']).parent.name}/{row['task']}: {row['disposition']} — {' '.join(reason.split())[:240]}")
        if len(issues) > limit:
            lines.append(f'  {len(issues)-limit} additional issues in JSON')
        for row in case['worker_failures'][:max(0, limit-len(issues))]:
            native = row.get('native_diagnostics', [])
            reason = ' | '.join(native[0]['lines']) if native else f"exit {row['exitcode']}"
            lines.append(f"  {Path(row['file']).stem}: worker failure — {' '.join(reason.split())[:240]}")
        for error in case['errors'][:3]:
            lines.append('  READ/VALIDATION ERROR: ' + error[:240])
    lines.append('\nNo jobs, logs or task states changed. Reusable means checked artifacts/geometry, not final scientific acceptance.')
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root')
    parser.add_argument('--verify', action='store_true', help='Read and verify declared artifact hashes and native results')
    parser.add_argument('--json', required=True, type=Path, help='Report destination outside the calculation tree')
    parser.add_argument('--limit', type=int, default=12)
    args = parser.parse_args(argv)
    try:
        if args.json.resolve().is_relative_to(Path(args.root).resolve()):
            raise ValueError('Write the report outside the run tree')
        report = recap(args.root, verify=args.verify)
        from kinbot.execution import atomic_json
        atomic_json(args.json, report)
        print(render(report, max(0, args.limit)))
        print(f'Full report: {args.json}')
    except (OSError, ValueError, RuntimeError) as error:
        print(f'Recap unavailable: {error}. No run files or jobs changed.', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
