"""Select native evidence for unresolved homolysis checks; never drive QC."""
import json
import os
from pathlib import Path
import re

RESOLVED = {'superseded_by_validated_saddle', 'parent_not_retained',
            'sampled_path_supports_provisional_pst'}
# These trees are not part of ordinary L1/L2 discovery. Retained ordinary
# restart/SCF/timeout archives, conf, hir and perm are deliberately traversed.
SKIP_DIRS = {'scratch', '__pycache__', '.git', '.venv', 'node_modules',
             'anl_automatic', 'vrctst', 'rotdPy', 'rotdpy', 'molpro'}


def workers(scopes):
    result = []
    for scope in scopes:
        path = scope/'homolysis_checks.json'
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text())
            pending = [c for c in data['checks'] if c.get('outcome') not in RESOLVED]
            if pending:
                result.append(scope)
        except (OSError, ValueError, KeyError, TypeError, AttributeError):
            # Capture the malformed ledger and its native evidence too.
            result.append(scope)
    return result


def sources(scopes, *, all_workers=False, progress=None):
    result = set()
    # A case already contains its workers. Walking both doubles NFS directory
    # traffic (and another inventory is deliberately taken after copying).
    roots = []
    for scope in sorted(set(scopes if all_workers else workers(scopes)),
                        key=lambda p: (len(p.parts), str(p))):
        if not any(scope.is_relative_to(parent) and not any(
                part in SKIP_DIRS or part.startswith(('.salvage_code.', '.fresh_code.'))
                for part in scope.relative_to(parent).parts) for parent in roots):
            roots.append(scope)
    for worker in roots:
        def failed(error):
            raise error  # An unreadable subtree cannot silently be omitted.
        for directory, dirs, files in os.walk(worker, followlinks=False, onerror=failed):
            parent = Path(directory)
            kept = []
            for name in dirs:
                path = parent/name
                if path.is_symlink():
                    result.add(path)  # Existing capture records refused symlinks.
                elif name not in SKIP_DIRS and not name.startswith(('.salvage_code.', '.fresh_code.')):
                    kept.append(name)
            dirs[:] = kept
            # Include complete ordinary files, including checkpoints, native
            # Hessians, trajectories, generation receipts and archived attempts.
            # Raw scratch wavefunctions are not needed to inspect convergence.
            result.update(parent/name for name in files
                          if not any(part in {'rwf', 'wfu'} for part in name.split('.')))
            if progress:
                progress(files_found=len(result), directory=str(parent))
    return result


def summarize(output, scopes):
    """Read only the exported copies; no QC readers, DB ingestion or retries."""
    lines = ['Inconclusive homolysis evidence capture',
             'Search outcomes below are recorded states, not export holds or retry recommendations.']
    records = []
    for scope in scopes:
        path=scope/'homolysis_checks.json'
        if not path.is_file():continue
        try:
            checks = json.loads(path.read_text())['checks']
            state_path=scope/'.kinbot_search_state.json'
            state=json.loads(state_path.read_text()) if state_path.is_file() else {}
            for check in checks:
                if check.get('outcome') in RESOLVED:
                    continue
                name = check.get('search', '')
                if not name or Path(name).name!=name:
                    raise ValueError('Missing or invalid search name in captured ledger.')
                entry = state.get('reactions', {}).get(name, {})
                item = dict(worker=str(path.parent.relative_to(output/'run')),
                            check=check, search_state=entry)
                scan = path.parent/(name+'.dissociation_scan.json')
                if scan.is_file():
                    raw = json.loads(scan.read_text())
                    item['scan'] = {k:v for k,v in raw.items() if k != 'points'}
                    item['point_count'] = len(raw.get('points', []))
                records.append(item)
                lines.extend(['', 'WORKER '+item['worker'],
                              'CHECK '+json.dumps(check, sort_keys=True),
                              'SEARCH '+json.dumps(entry, sort_keys=True),
                              'SCAN '+json.dumps(item.get('scan'), sort_keys=True),
                              'POINTS '+str(item.get('point_count', 0))])
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
            lines.extend(['', 'UNREADABLE_EXPORTED_RECORD '+str(path.relative_to(output)), str(error)])
    lines.extend(['', 'Native files are complete copies, not log tails.',
                  'See manifest.json for changed, unavailable or oversized files.',
                  'No search was reset, retried, accepted or rejected by this capture.'])
    return records, '\n'.join(lines)+'\n'


def scheduler(output, records):
    """Bounded read-only queries for recorded scan/TS/IRC allocations only."""
    from datetime import datetime, timedelta, timezone
    from kinbot.run_outage_audit import command
    ids = set()
    errors = []
    def collect(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in {'job_id', 'slurm_job_id', 'previous_job_id', 'submitted_job_id'}:
                    if re.fullmatch(r'\d+(?:_\d+)?', str(item)):
                        ids.add(str(item))
                collect(item)
        elif isinstance(value, list):
            for item in value:collect(item)
    for record in records:
        collect(record['search_state'])
        worker=output/'run'/record['worker']
        check=record['check']
        search=check.get('search','')
        if not search or Path(search).name!=search:continue
        for path in (worker/(search+'.dissociation_scan.json'),worker/'.kinbot_theory_jobs.json'):
            if not path.is_file():continue
            try:
                data=json.loads(path.read_text())
                if path.name=='.kinbot_theory_jobs.json':
                    prefixes=[n for n in [search,check.get('parent'),*check.get('saddles',[])] if n]
                    data=[entry for name,entry in data.get('jobs',{}).items()
                          if any(Path(name).name==n or Path(name).name.startswith(n+'_') for n in prefixes)]
                collect(data)
            except (OSError,ValueError,KeyError,TypeError) as error:
                errors.append(dict(path=str(path),error=str(error)))
    ordered=sorted(ids)
    start=(datetime.now(timezone.utc)-timedelta(days=90)).date().isoformat()
    result=dict(job_ids=ordered,since=start,read_errors=errors,accounting=[],
                note='Accounting can be unavailable or purged; that is not a failed-QC verdict.')
    for offset in range(0,len(ordered),100):
        print(f'Reading Slurm accounting for recorded jobs {offset+1}-{min(offset+100,len(ordered))}/{len(ordered)}.',flush=True)
        result['accounting'].append(command(['sacct','-n','-P','--duplicates','-S',start,
            '-j',','.join(ordered[offset:offset+100]),
            '--format=JobID,JobName%200,State,ExitCode,Submit,Start,End,Elapsed,NodeList,WorkDir%500']))
    return result
