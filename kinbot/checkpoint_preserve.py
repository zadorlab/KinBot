"""Copy named checkpoint artifacts from recorded scratch on the current node.

Run as a new step within a STOP-signalled allocation, using a run_outage_audit
report. No calculation is resumed, signalled, accepted or deleted. Unreferenced
scratch directories are listed for later investigation, never assumed to belong
to this run. This is checkpoint salvage, not an integral-file/full-scratch backup.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import re
import socket
import stat
import subprocess
import sys
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed

from kinbot.run_outage_audit import CHECKPOINTS, MIB


def file_hash(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(8*MIB), b''):
            result.update(data)
    return result.hexdigest()


def host_matches(left, right):
    left, right = str(left).rstrip('.').lower(), str(right).rstrip('.').lower()
    return left == right or (('.' not in left or '.' not in right)
                              and left.split('.')[0] == right.split('.')[0])


def checkpoint(path):
    return (path.suffix.lower() in {'.wfu', '.chk', '.fchk', '.rwf', '.gbw',
                                  '.traj', '.xyz', '.xml', '.molden', '.grad',
                                  '.inp', '.com', '.json'}
            or path.name in CHECKPOINTS | {'MINP', 'ZMAT', 'EXIT'}
            or path.name.endswith('.wfu.failed'))


def preserve(report_path, output_parent, *, max_gib=128, include_unattributed=False):
    report_path = Path(report_path).resolve(strict=True)
    report = json.loads(report_path.read_text())
    run = Path(report['root']).resolve(strict=True)
    output_parent = Path(output_parent).expanduser().resolve()
    if output_parent == run or output_parent.is_relative_to(run):
        raise ValueError('Checkpoint preservation must be outside the source run.')
    if not math.isfinite(max_gib) or max_gib <= 0:
        raise ValueError('The copy budget must be positive.')
    host = socket.gethostname()
    locations = {}
    for item in report.get('scratch_locations', []):
        if not host_matches(item.get('scratch_host', item.get('execution_node', '')), host):
            continue
        value = item.get('scratch_directory')
        if not value and item.get('native_work_directory'):
            value = str(Path(item['native_work_directory']).parent)
        if not value:
            continue
        path = Path(value)
        # A runtime record is evidence, not authority to copy arbitrary files.
        if not path.is_absolute() or not path.name.startswith('kinbot-'):
            continue
        locations.setdefault(str(path), []).append(item['source'])
    unattributed = []
    for parent in sorted({Path(p).parent for p in locations}):
        if not parent.is_dir() or parent.resolve() != parent:
            continue
        for child in sorted(parent.glob('kinbot-*')):
            if child.is_symlink() or str(child) in locations:
                continue
            try:
                info = child.stat()
                if child.is_dir() and info.st_uid == os.getuid():
                    unattributed.append(dict(path=str(child), mtime_ns=info.st_mtime_ns,
                        attribution='unknown; do not reuse without run/attempt provenance'))
                    if include_unattributed:
                        locations[str(child)] = []
            except OSError:
                continue
    output_parent.mkdir(parents=True, exist_ok=True)
    dest = output_parent/(host.split('.')[0]+'-'+uuid.uuid4().hex[:8])
    dest.mkdir(mode=0o700)
    result = dict(schema=1, host=host, slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        source_report=str(report_path), source_report_sha256=file_hash(report_path),
        started_utc=datetime.now(timezone.utc).isoformat(), source_run=str(run),
        files=[], issues=[], missing_scratch=[], unreferenced_scratch=unattributed,
        limitations=['Named checkpoints and diagnostics only; integral scratch is not backed up.',
            'A verified byte copy does not prove the native program can restart it.',
            'Old scratch on other nodes requires a separate node capture.'])
    try:
        from kinbot.run_cleanup import _processes
        result['run_processes_on_node'] = _processes(run)
    except (OSError, RuntimeError) as error:
        result['process_inventory_error'] = str(error)
    receipt = dest/'preservation.json'
    def save():
        temporary = dest/'preservation.json.tmp'
        temporary.write_text(json.dumps(result, indent=2)+'\n')
        temporary.replace(receipt)
    save()
    spent, count = 0, 0
    for value, sources in sorted(locations.items()):
        path = Path(value)
        if not path.exists():
            result['missing_scratch'].append(value)
            continue
        if (path.resolve() != path or not path.is_dir() or path.stat().st_uid != os.getuid()
                or path == Path.home().resolve() or run.is_relative_to(path)
                or output_parent.is_relative_to(path)):
            result['issues'].append(dict(path=value, error='Unsafe scratch directory; not copied.'))
            continue
        key = hashlib.sha256(value.encode()).hexdigest()[:16]
        def walk_error(error):
            result['issues'].append(dict(path=error.filename, error=str(error)))
        for parent, dirs, files in os.walk(path, followlinks=False, onerror=walk_error):
            dirs[:] = sorted(d for d in dirs if not (Path(parent)/d).is_symlink())
            for name in sorted(files):
                source = Path(parent)/name
                count += 1
                if count > 100000:
                    raise RuntimeError(f'File count limit reached; partial receipt: {receipt}')
                try:
                    info = source.lstat()
                    if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid():
                        continue
                    log = source.suffix.lower() in {'.out', '.log', '.stdout', '.stderr', '.err'}
                    if not checkpoint(source) and not log:
                        continue
                    amount = min(info.st_size, MIB) if log else info.st_size
                    if spent+amount > max_gib*1024*MIB:
                        result['issues'].append(dict(path=str(source), error='Copy budget exceeded; not copied.'))
                        continue
                    if shutil.disk_usage(dest).free < amount+1024*MIB:
                        raise OSError('Insufficient destination free space (1 GiB reserve).')
                    # Recheck every ancestor immediately before opening.
                    if source.resolve(strict=True) != source:
                        raise OSError('Symlink source is not eligible.')
                    target = dest/'files'/key/source.relative_to(path)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    digest = hashlib.sha256()
                    fd = os.open(source, os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0))
                    with os.fdopen(fd, 'rb') as src, target.open('xb') as out:
                        before = os.fstat(src.fileno())
                        if not stat.S_ISREG(before.st_mode):
                            raise OSError('Source changed to a nonregular file.')
                        offset = max(0, before.st_size-MIB) if log else 0
                        size = before.st_size-offset
                        if spent+size > max_gib*1024*MIB:
                            raise OSError('File grew past copy budget.')
                        src.seek(offset)
                        remaining = size
                        while remaining:
                            data = src.read(min(8*MIB, remaining))
                            if not data:
                                raise OSError('Source shortened during copy.')
                            out.write(data); digest.update(data); remaining -= len(data)
                        out.flush(); os.fsync(out.fileno())
                        after = os.fstat(src.fileno())
                    current = source.stat(follow_symlinks=False)
                    signature = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
                    stable = signature(before) == signature(after) == signature(current)
                    copied_hash = file_hash(target)
                    result['files'].append(dict(source=str(source), target=str(target),
                        runtime_sources=sources, attribution='recorded' if sources else 'unknown',
                        bytes=size, offset=offset,
                        mtime_ns=before.st_mtime_ns, sha256=digest.hexdigest(),
                        verified_copy=copied_hash == digest.hexdigest(), source_stable=stable))
                    spent += size
                    if not stable or copied_hash != digest.hexdigest():
                        result['issues'].append(dict(path=str(source), error='Copy verification or source stability failed.'))
                    save()
                except OSError as error:
                    result['issues'].append(dict(path=str(source), error=str(error)))
                    save()
    result.update(finished_utc=datetime.now(timezone.utc).isoformat(), copied_bytes=spent,
                  status='review_required', selected_scratch_count=len(locations))
    save()
    return dict(receipt=str(receipt), host=host, files=len(result['files']),
        copied_bytes=spent, issues=len(result['issues']), missing_scratch=len(result['missing_scratch']),
        unreferenced_scratch=len(result['unreferenced_scratch']),
        status='review_required', cancellation_authorized_by_receipt=False)


def allocations(report_path, freeze_path, output_parent, *, max_gib=128):
    """Use existing frozen allocations for copies; never submit or cancel jobs."""
    from kinbot.run_cleanup import _jobs, _processes
    report_path, freeze_path = Path(report_path).resolve(), Path(freeze_path).resolve()
    report, freeze = (json.loads(p.read_text()) for p in (report_path, freeze_path))
    root = Path(report['root']).resolve(strict=True)
    if Path(freeze['root']).resolve() != root:
        raise ValueError('Freeze receipt and report name different runs.')
    live = _processes(root)
    if any(p['state'] not in ('T', 't') for p in live):
        raise RuntimeError('Run-scoped local processes are not all stopped; no copy step launched.')
    stopped = {str(j) for a in freeze.get('actions', [])
               if a.get('action') == 'freeze_running' and a.get('returncode') == 0
               for j in a['jobs']}
    jobs = _jobs(root)  # A failed scheduler query must stop this operation.
    by_node, omitted = {}, []
    for job in jobs:
        node = job.get('reason', '')
        if (job['job_id'] not in stopped or job['state'] != 'RUNNING'
                or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', node)):
            omitted.append(job)
            continue
        by_node.setdefault(node, job)
    if not by_node:
        raise RuntimeError('No currently running, STOP-signalled single-node allocation found.')
    output_parent = Path(output_parent).expanduser().resolve()
    if output_parent == root or output_parent.is_relative_to(root):
        raise ValueError('Copy destination must be outside the source run.')
    output_parent.mkdir(parents=True, exist_ok=True)
    session = output_parent/('copy-steps-'+uuid.uuid4().hex[:8]); session.mkdir(mode=0o700)
    result = dict(root=str(root), freeze_receipt=str(freeze_path),
                  report=str(report_path), steps=[], omitted_jobs=omitted,
                  status='copying', cancels_jobs=False)
    def save():
        (session/'steps.json').write_text(json.dumps(result, indent=2)+'\n')
    save()
    def one(node, job):
        command = ['srun', '--jobid='+job['job_id'], '--overlap', '--exact',
            '--nodes=1', '--ntasks=1', '--cpus-per-task=1', '--immediate=30',
            sys.executable, '-B', '-m', 'kinbot.checkpoint_preserve', str(report_path),
            '--output-parent', str(output_parent/'nodes'), '--max-gib', str(max_gib),
            '--include-unattributed']
        stdout, stderr = session/(node+'.stdout'), session/(node+'.stderr')
        print('COPY_START', node, job['job_id'], flush=True)
        try:
            with stdout.open('x') as out, stderr.open('x') as err:
                value = subprocess.run(command, stdout=out, stderr=err, check=False)
            return dict(node=node, job=job, argv=command, returncode=value.returncode,
                        stdout=str(stdout), stderr=str(stderr))
        except OSError as error:
            return dict(node=node, job=job, argv=command, error=str(error))
    # Two simultaneous copies bound shared-filesystem pressure. They consume
    # existing allocations, never request an additional node or chemistry task.
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(one, node, job) for node, job in sorted(by_node.items())]
        for future in as_completed(futures):
            item = future.result(); result['steps'].append(item); save()
            print('COPY_RESULT', json.dumps(item), flush=True)
    result['status'] = 'review_required'; save()
    return dict(receipt=str(session/'steps.json'), steps=len(result['steps']),
                failed_steps=sum(s.get('returncode') != 0 for s in result['steps']),
                omitted_jobs=len(omitted), status='review_required', cancels_jobs=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--output-parent', required=True, type=Path)
    parser.add_argument('--max-gib', type=float, default=128)
    parser.add_argument('--include-unattributed', action='store_true',
        help='Also copy named checkpoints from owned sibling kinbot-* workspaces; keep attribution unknown.')
    parser.add_argument('--allocations', type=Path, metavar='FREEZE_RECEIPT',
        help='On the controller host, run copy steps in allocations recorded as STOP-signalled.')
    args = parser.parse_args(argv)
    result = (allocations(args.report, args.allocations, args.output_parent, max_gib=args.max_gib)
              if args.allocations else preserve(args.report, args.output_parent,
                  max_gib=args.max_gib, include_unattributed=args.include_unattributed))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
