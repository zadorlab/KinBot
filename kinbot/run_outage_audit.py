"""Capture bounded, read-only outage evidence without running any owner/recovery.

The archive contains metadata, small inputs/trajectories and bounded output
tails. Large checkpoint files and databases are inventoried, NOT backed up.
It is a sequential observation of a live tree, not a consistent filesystem
snapshot or a scientific acceptance/restart decision.
"""
import argparse
from collections import Counter
from datetime import datetime, timedelta, timezone
import getpass
import hashlib
import io
import json
import os
from pathlib import Path
import re
import socket
import stat
import subprocess
import tarfile
import uuid


MIB = 1024 * 1024
SKIP = {'.git', '.venv', '__pycache__', 'scratch', 'node_modules'}
CHECKPOINTS = {'fort.16', 'fort.17', 'JOBARC', 'JAINDX', 'OPTARC', 'FCMINT',
               'FCMFINAL', 'GRD', 'HESSIAN'}
INPUTS = {'.inp', '.com', '.xyz', '.slurm', '.sh'}
LOGS = {'.out', '.log', '.err', '.stdout', '.stderr'}
JOB_ID = re.compile(r'\d+(?:_\d+)?\Z')


def utc():
    return datetime.now(timezone.utc).isoformat()


def command(argv):
    """Only callers below construct commands; never execute captured content."""
    try:
        value = subprocess.run(argv, capture_output=True, text=True, timeout=20)
        return dict(argv=argv, returncode=value.returncode,
                    status='observed' if value.returncode == 0 else 'unknown',
                    stdout=value.stdout[-8*MIB:], stderr=value.stderr[-65536:])
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(argv=argv, status='unknown', error=str(error))


def scheduler(root, jobs, start):
    queue = command(['squeue', '-h', '-u', getpass.getuser(), '-o',
                     '%i|%T|%N|%M|%Z'])
    result = dict(queue=queue, scoped_jobs=[], details=[], accounting=[])
    if queue['status'] == 'observed':
        for line in queue['stdout'].splitlines():
            parts = line.split('|', 4)
            if len(parts) != 5:
                queue['status'] = 'unknown'
                queue['parse_error'] = 'Unrecognized squeue row.'
                continue
            ident, state, nodes, elapsed, directory = map(str.strip, parts)
            # Lexical containment: do not follow a remote WorkDir symlink.
            if Path(directory).is_absolute() and Path(directory).is_relative_to(root):
                result['scoped_jobs'].append(dict(job_id=ident, state=state,
                    nodes=nodes, elapsed=elapsed, directory=directory))
    ids = list(dict.fromkeys([*jobs, *(j['job_id'] for j in result['scoped_jobs'])]))
    valid = [j for j in ids if JOB_ID.fullmatch(j)]
    result['omitted_ids'] = [j for j in ids if j not in valid] + valid[100:]
    valid = valid[:100]
    # Explicit IDs include jobs already absent from squeue/scontrol. Duplicate
    # accounting records can distinguish several starts under one job number.
    for offset in range(0, len(valid), 50):
        result['accounting'].append(command(['sacct', '-n', '-P', '--duplicates',
            '-u', getpass.getuser(), '-S', start, '-j', ','.join(valid[offset:offset+50]),
            '--format=JobID,JobName,State,ExitCode,Submit,Start,End,Elapsed,NodeList,WorkDir']))
    # One bounded request per explicitly named job; no user-wide scontrol scan.
    for ident in jobs:
        result['details'].append(command(['scontrol', 'show', 'job', '-o', ident]))
    return result


def inventory(root, issues, maximum):
    paths = []
    def error(value):
        issues.append(dict(path=value.filename, error=str(value)))
    for directory, dirs, files in os.walk(root, followlinks=False, onerror=error):
        parent = Path(directory)
        kept = []
        for name in sorted(dirs):
            path = parent/name
            if path.is_symlink():
                issues.append(dict(path=str(path.relative_to(root)), excluded='symlink'))
            elif name not in SKIP and not name.startswith(('.salvage_code.', '.fresh_code.')):
                kept.append(name)
        dirs[:] = kept
        for name in sorted(files):
            path = parent/name
            try:
                info = path.lstat()
                if not stat.S_ISREG(info.st_mode):
                    issues.append(dict(path=str(path.relative_to(root)), excluded='not_regular'))
                    continue
                paths.append((path, info))
            except OSError as exc:
                error(exc)
            if len(paths) >= maximum:
                issues.append(dict(error='File inventory limit reached; capture is incomplete.'))
                return paths
    return paths


def read_file(path, root, limit, *, tail=False):
    """Reject symlinks, including ancestor substitutions, and detect writes."""
    if path.resolve(strict=True) != path or not path.is_relative_to(root):
        raise ValueError('Refusing a symlink or off-root source.')
    fd = os.open(path, os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0))
    with os.fdopen(fd, 'rb') as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError('Source is not a regular file.')
        if before.st_size > limit and not tail:
            raise ValueError(f'Full-file capture exceeds {limit} bytes.')
        offset = max(0, before.st_size-limit) if tail else 0
        stream.seek(offset)
        data = stream.read(limit)
        after = os.fstat(stream.fileno())
    current = path.stat(follow_symlinks=False)
    signature = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    return data, dict(offset=offset, captured_bytes=len(data), bytes=before.st_size,
        mtime_ns=before.st_mtime_ns,
        changed=signature(before) != signature(after) or signature(before) != signature(current),
        sha256=hashlib.sha256(data).hexdigest())


def add_bytes(archive, name, data):
    item = tarfile.TarInfo(name)
    item.size, item.mode = len(data), 0o600
    archive.addfile(item, io.BytesIO(data))


def scratch_refs(value, origin, output):
    if isinstance(value, dict):
        fields = {k: value[k] for k in ('scratch_host', 'scratch_directory',
                  'native_work_directory', 'execution_node', 'slurm_job_id') if k in value}
        if 'scratch_directory' in fields or 'native_work_directory' in fields:
            output.append(dict(source=origin, **fields, existence='not_checked_on_node'))
        for child in value.values():
            scratch_refs(child, origin, output)
    elif isinstance(value, list):
        for child in value:
            scratch_refs(child, origin, output)


def describe(root, records, queue):
    graphs, campaigns, phases = [], [], []
    for name, value in records.items():
        if not isinstance(value, dict):
            continue
        path = Path(name)
        if path.name == 'state.json' and str(path.parent/'workflow.json') in records:
            workflow = records[str(path.parent/'workflow.json')]
            if not isinstance(workflow, dict):
                continue
            definitions = {t.get('id'): t for t in workflow.get('tasks', []) if isinstance(t, dict)}
            entries = value.get('tasks', {})
            if not isinstance(entries, dict):
                continue
            tasks = []
            for ident, entry in entries.items():
                if not isinstance(entry, dict):
                    continue
                definition = definitions.get(ident, {})
                directory = path.parent/'tasks'/ident
                execution = records.get(str(directory/'execution.json'), {})
                owners = [j for j in queue.get('scoped_jobs', [])
                          if j['directory'] == str(root/directory)]
                tasks.append(dict(id=ident, ledger=entry, kind=definition.get('kind'),
                    backend=definition.get('backend'), profile=definition.get('profile'),
                    execution_status=execution.get('status') if isinstance(execution, dict) else None,
                    observed_allocations=owners,
                    acceptance='not_revalidated_by_capture',
                    restart='requires_native_evidence_and_owner_reconciliation'))
            graphs.append(dict(path=str(path.parent),
                historical=any(p in ('attempts', 'archive', 'archives') for p in path.parts),
                counts=dict(Counter(e.get('status', 'unknown') for e in entries.values() if isinstance(e, dict))),
                tasks=tasks))
        elif path.name == 'review_launch.json':
            campaigns.append(dict(path=name, record=value))
        elif path.name == 'pes_exploration.json':
            from kinbot.pes_phase import process_identity, same_process_identity
            owner = value.get('owner_identity')
            current = process_identity(value.get('owner_pid')) if owner else None
            workers = []
            for ident, worker in value.get('workers', {}).items():
                saved = worker.get('identity')
                observed = process_identity(worker.get('pid')) if saved else None
                workers.append(dict(worker=ident, registered=worker, current_identity=observed,
                    identity_matches=same_process_identity(saved, observed)))
            phases.append(dict(path=name, phase=value.get('phase'), token=value.get('token'),
                owner_identity=owner, current_identity=current,
                owner_identity_matches=same_process_identity(owner, current), workers=workers))
    return dict(graphs=graphs, pes=phases, review_campaigns=campaigns)


def capture(root, output_parent, *, jobs=(), since=None, query_scheduler=True,
            max_files=300000, metadata_mb=192, logs_mb=192):
    root = Path(root).expanduser().resolve(strict=True)
    output_parent = Path(output_parent).expanduser().resolve()
    if not root.is_dir() or root == Path('/') or root == Path.home().resolve():
        raise ValueError('Expected one calculation run directory.')
    if output_parent == root or output_parent.is_relative_to(root):
        raise ValueError('Write the audit outside the run being observed.')
    if any(not JOB_ID.fullmatch(j) for j in jobs) or len(jobs) > 20:
        raise ValueError('Supply at most 20 numeric job IDs (array elements allowed).')
    since = since or (datetime.now()-timedelta(days=14)).strftime('%Y-%m-%d')
    report = dict(schema=1, started_utc=utc(), root=str(root), host=socket.gethostname(),
        issues=[], files=[], scratch_locations=[],
        limitations=['Sequential live observations, not an atomic snapshot.',
            'Checkpoint/database/pickle entries are metadata only; no native checkpoint backup.',
            'No accepted result is revalidated; no restart is applied.',
            'Missing scheduler evidence is unknown, not proof a job ended.'])
    report['scheduler'] = (scheduler(root, list(jobs), since) if query_scheduler
                           else dict(queue={'status': 'not_queried'}, scoped_jobs=[]))
    try:
        from kinbot.run_cleanup import _processes
        report['processes'] = dict(status='observed', records=_processes(root))
    except (OSError, RuntimeError) as error:
        report['processes'] = dict(status='unknown', error=str(error))
    paths = inventory(root, report['issues'], max_files)
    # Metadata/trajectory evidence is budgeted independently from verbose logs.
    priority = lambda pair: (0 if pair[0].suffix == '.json' else
                             1 if pair[0].suffix == '.traj' else
                             2 if pair[0].suffix in INPUTS else 3, str(pair[0]))
    output_parent.mkdir(parents=True, exist_ok=True)
    target = output_parent/('outage-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    target.mkdir(mode=0o700)
    archive_path = target/'outage-evidence.tar.gz'
    records, spent = {}, dict(metadata=0, logs=0)
    budgets = dict(metadata=int(metadata_mb*MIB), logs=int(logs_mb*MIB))
    with tarfile.open(archive_path, 'w:gz') as archive:
        for path, info in sorted(paths, key=priority):
            rel = str(path.relative_to(root))
            item = dict(path=rel, bytes=info.st_size, mtime_ns=info.st_mtime_ns)
            suffix = path.suffix.lower()
            checkpoint = suffix in {'.wfu', '.chk', '.fchk', '.rwf', '.gbw'} or path.name in CHECKPOINTS or path.name.endswith('.wfu.failed')
            item['kind'] = 'checkpoint_metadata' if checkpoint else 'file_metadata'
            report['files'].append(item)
            if checkpoint or suffix in {'.db', '.pickle', '.pkl', '.tar', '.gz', '.zip'}:
                continue
            if suffix not in {'.json', '.traj', *INPUTS, *LOGS} and path.name not in {'MINP', 'ZMAT', 'EXIT'}:
                continue
            is_log = suffix in LOGS
            budget = 'logs' if is_log else 'metadata'
            limit = 65536 if is_log else (32*MIB if suffix == '.traj' else 8*MIB)
            amount = min(info.st_size, limit) if is_log else info.st_size
            if amount > limit or spent[budget]+amount > budgets[budget]:
                item['omitted'] = 'size_or_capture_budget'
                continue
            try:
                data, metadata = read_file(path, root, limit, tail=is_log)
                if spent[budget]+len(data) > budgets[budget]:
                    item['omitted'] = 'capture_budget_after_growth'
                    continue
                spent[budget] += len(data)
                member = 'files/'+rel+('.tail' if metadata['offset'] else '')
                add_bytes(archive, member, data)
                item.update(metadata, member=member, kind='tail' if metadata['offset'] else 'full_file')
                if suffix == '.json':
                    try:
                        value = json.loads(data)
                        scratch_refs(value, rel, report['scratch_locations'])
                        if not metadata['changed']:
                            records[rel] = value
                    except (ValueError, RecursionError) as error:
                        report['issues'].append(dict(path=rel, error=str(error)))
            except (OSError, ValueError) as error:
                report['issues'].append(dict(path=rel, error=str(error)))
        report.update(describe(root, records, report['scheduler']))
        report['finished_utc'] = utc()
        report['captured_bytes'] = spent
        report['omitted_content_count'] = sum('omitted' in f for f in report['files'])
        report['changed_files'] = [f['path'] for f in report['files'] if f.get('changed')]
        data = json.dumps(report, indent=2).encode()
        add_bytes(archive, 'report.json', data)
    (target/'report.json').write_bytes(data)
    digest = hashlib.sha256()
    with archive_path.open('rb') as stream:
        for block in iter(lambda: stream.read(MIB), b''):
            digest.update(block)
    (target/'SHA256SUMS').write_text(f'{digest.hexdigest()}  {archive_path.name}\n')
    summary = dict(directory=str(target), archive=str(archive_path), sha256=digest.hexdigest(),
        files=len(report['files']), issues=len(report['issues']),
        omitted_content=report['omitted_content_count'], changed_files=len(report['changed_files']),
        graphs=len(report['graphs']), scratch_locations=len(report['scratch_locations']),
        scheduler=report['scheduler']['queue']['status'],
        checkpoint_backup=False, restart_applied=False)
    (target/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--output-parent', type=Path, required=True)
    parser.add_argument('--job', action='append', default=[])
    parser.add_argument('--since')
    parser.add_argument('--no-scheduler', action='store_true')
    args = parser.parse_args(argv)
    print(json.dumps(capture(args.run, args.output_parent, jobs=args.job,
                            since=args.since, query_scheduler=not args.no_scheduler), indent=2))


if __name__ == '__main__':
    main()
