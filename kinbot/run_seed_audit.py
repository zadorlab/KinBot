"""Export seed/restart evidence without changing or driving the source run.

Capture shallow PES/search ledgers, full SQLite database files and their sidecars,
graph state, and explicitly selected interrupted tasks. This is an audit bundle,
not a complete run backup or a declaration that any result is reusable.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import tarfile
import time
import traceback
import uuid

from kinbot.execution import atomic_json, file_sha256
from kinbot.run_outage_audit import read_file

MIB = 1024**2
METADATA = {'.kinbot_search_state.json', '.kinbot_theory_jobs.json',
            '.kinbot_slurm_jobs.json', '.kinbot_slurm_nodes.json', 'input.json',
            'pes_exploration.json', 'accepted_l2_network.json', 'kinbot_monitor.out',
            'network_request.json', 'network_revisit.json'}
TASK_NAMES = {'task.json', 'execution.json', 'job.slurm', 'geometry.xyz', 'final.xyz',
              'continuation.xyz', 'optimization.log', 'optimization.traj',
              'optimization_budget.json', 'optimization_resume.json',
              'slurm.stdout', 'slurm.stderr'}
TEXT_SUFFIXES = {'.inp', '.com', '.out', '.log', '.xml', '.xyz', '.stdout', '.stderr', '.err'}


class Progress:
    """Flush phase changes and bounded progress to a terminal or detached log."""
    def __init__(self, callback=None):
        self.callback = callback
        self.phase = None
        self.last = float('-inf')

    def __call__(self, phase, **details):
        now = time.monotonic()
        if phase == self.phase and now - self.last < 5:
            return
        self.phase, self.last = phase, now
        record = dict(updated_utc=datetime.now(timezone.utc).isoformat(),
                      phase=phase, **details)
        print('CAPTURE_PROGRESS '+json.dumps(record, sort_keys=True), flush=True)
        if self.callback:
            self.callback(record)


def _scopes(root):
    """Known PES layout only; never recurse through millions of sample files."""
    cases = {root} if (root/'input.json').is_file() else {
        p.parent for p in root.glob('*/input.json') if not p.is_symlink()}
    scopes = set(cases)
    for case in cases:
        for child in case.iterdir():
            if child.is_dir() and not child.is_symlink() and (
                    (child/'kinbot.db').is_file() or (child/(child.name+'.json')).is_file()):
                scopes.add(child)
    return sorted(scopes)


def sources(root, tasks=(), evidence=None, *, mess_handoff=False, homolysis_native=False,
            ordinary_native=False, progress=None):
    paths = set()
    for name in ('prepared.json', '.kinbot_resource_policy.json', '.kinbot_live_control.json'):
        if (root/name).is_file():paths.add(root/name)
    for scope in _scopes(root):
        if progress:progress(directory=str(scope), files_found=len(paths))
        for name in METADATA | {scope.name+'.json', 'kinbot.db', 'kinbot.db-wal', 'kinbot.db-shm'}:
            if (scope/name).is_file():paths.add(scope/name)
        if mess_handoff:
            # Complete, shallow assembly handoff: no QC/sampling tree traversal.
            # Full databases retain atom mapping, selected observations and all
            # conformer records needed to replay final PES energy selection.
            for pattern in ('*.json', '*.mess', 'summary_*.out'):
                paths.update(scope.glob(pattern))
            for name in ('chemids', 'pes.log', 'driver.log', 'kinbot.log', 'kinbot.err',
                         'pes_summary.txt', 'reactionList.log'):
                if (scope/name).is_file():paths.add(scope/name)
            for folder in ('me', 'xyz'):
                paths.update(p for p in (scope/folder).glob('*') if p.is_file())
            continue
        paths.update(scope.glob('vrctst/orbital_review/*/request.json'))
        graphs = {*scope.glob('anl_automatic/states/*/composite/workflow.json'),
                  *scope.glob('vrctst/*/*/workflow.json')}
        for workflow in graphs:
            graph = workflow.parent
            paths.update(graph.glob('*.json'))
            for task in (graph/'tasks').iterdir() if (graph/'tasks').is_dir() else ():
                paths.update(p for name in ('task.json', 'execution.json')
                             if (p := task/name).is_file())
    for relative in tasks:
        path = root/relative
        if Path(relative).is_absolute() or not path.resolve().is_relative_to(root):
            raise ValueError('Selected task must be relative to the source run.')
        if not (path/'task.json').is_file():raise FileNotFoundError(path/'task.json')
        # These are explicitly selected small task trees, excluding binary
        # checkpoint data. No general recursive run traversal is performed.
        for parent, dirs, files in os.walk(path, followlinks=False):
            dirs[:] = [d for d in dirs if d not in ('scratch', '__pycache__')
                       and not (Path(parent)/d).is_symlink()]
            for name in files:
                p=Path(parent)/name
                if name in TASK_NAMES or p.suffix in TEXT_SUFFIXES or p.suffix=='.json':paths.add(p)
    if homolysis_native or ordinary_native:
        from kinbot.homolysis_evidence import sources as native_sources
        paths.update(native_sources(_scopes(root), all_workers=ordinary_native, progress=progress))
        for scope in {root, *[p.parent for p in root.glob('*/input.json')]}:
            paths.update(scope.glob('*.json'))
            for name in ('pes.log', 'driver.log', 'pes_summary.txt', 'chemids'):
                if (scope/name).is_file():paths.add(scope/name)
    result = [(p, root, Path('run')/p.relative_to(root)) for p in sorted(paths)]
    if evidence is not None:
        extra = {p for pattern in ('retirement.json', 'outage-*/report.json',
            'checkpoints/copy-steps*/*.json', 'checkpoints/copy-steps*/*.stdout',
            'checkpoints/copy-steps*/*.stderr', 'checkpoints/nodes/*/preservation.json')
                 for p in evidence.glob(pattern)}
        for receipt in evidence.glob('checkpoints/nodes/*/preservation.json'):
            data=json.loads(receipt.read_text())
            for entry in data.get('files',[]):
                path=Path(entry['target'])
                if path.suffix in TEXT_SUFFIXES or path.name in ('MINP','ZMAT','EXIT'):
                    if not path.resolve().is_relative_to(evidence):
                        raise ValueError('Checkpoint evidence target escaped its supplied directory.')
                    extra.add(path)
        result.extend((p,evidence,Path('evidence')/p.relative_to(evidence)) for p in sorted(extra))
    return result


def _capture_paths(root, output_parent, evidence, max_mb):
    root=Path(root).resolve(strict=True);output_parent=Path(output_parent).resolve()
    evidence=Path(evidence).resolve(strict=True) if evidence else None
    if output_parent==root or output_parent.is_relative_to(root):
        raise ValueError('Write audit evidence outside the source run.')
    if evidence and (output_parent==evidence or output_parent.is_relative_to(evidence)):
        raise ValueError('Write this audit outside the supplied evidence tree.')
    if type(max_mb) is not int or max_mb<1:raise ValueError('max_mb must be positive.')
    return root, output_parent, evidence


def _verify_archive(archive, expected, progress):
    expected = dict(expected)
    total = len(expected)
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            if not member.isfile():continue
            with tar.extractfile(member) as stream:
                digest = hashlib.sha256()
                for block in iter(lambda: stream.read(MIB), b''):
                    digest.update(block)
            if digest.hexdigest() != expected.pop(member.name, None):
                raise RuntimeError('Archive verification mismatch: '+member.name)
            progress('verifying_archive', checked=total-len(expected), total=total)
    if expected:raise RuntimeError('Archive is missing captured evidence.')


def capture(root, output_parent, *, tasks=(), evidence=None, max_mb=2048, mess_handoff=False,
            homolysis_native=False, ordinary_native=False, progress_callback=None):
    root,output_parent,evidence = _capture_paths(root,output_parent,evidence,max_mb)
    progress = Progress(progress_callback)
    output=output_parent/('seed-audit-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    output.mkdir(parents=True)
    print('DESTINATION '+str(output), flush=True)
    progress('inventory', output_directory=str(output))
    homolysis_native = homolysis_native or ordinary_native
    selected=sources(root,tasks,evidence,mess_handoff=mess_handoff,homolysis_native=homolysis_native,
                     ordinary_native=ordinary_native,progress=lambda **d:progress('inventory',**d))
    if homolysis_native:
        print(f'Capturing {len(selected)} native/metadata files; complete files, not tails.',flush=True)
    def handoff_directories():
        scopes=_scopes(root)
        return {p for scope in scopes for p in (scope,*scope.iterdir())
                if p.is_dir() and not p.is_symlink()}
    directories=handoff_directories() if mess_handoff else set()
    if not selected:raise ValueError('No recognized seed evidence in source run.')
    entries=[];issues=[];used=0
    progress('copying', copied=0, total=len(selected))
    for index,(source,base,relative) in enumerate(selected,1):
        try:
            data,info=read_file(source,base,min(256*MIB,max_mb*MIB-used))
            if info['changed']:raise RuntimeError('Source changed during capture.')
            target=output/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
            entries.append(dict(source=str(source),path=str(relative),**info))
            used+=len(data)
        except (OSError,ValueError,RuntimeError) as error:
            issues.append(dict(source=str(source),error=str(error)))
        progress('copying', processed=index, copied=len(entries), total=len(selected),
                 bytes=used, issues=len(issues))
    # Includes DB/WAL sidecars: a change anywhere invalidates the captured
    # group. An actively changing run must be quiesced before seed import.
    progress('rechecking_source_hashes', checked=0, total=len(entries), bytes=used, issues=len(issues))
    for index,entry in enumerate(entries,1):
        source=Path(entry['source'])
        try:
            if source.stat().st_mtime_ns!=entry['mtime_ns'] or file_sha256(source)!=entry['sha256']:
                issues.append(dict(source=str(source),error='Source changed before capture completed.'))
        except OSError as error:issues.append(dict(source=str(source),error=str(error)))
        progress('rechecking_source_hashes', checked=index, total=len(entries), issues=len(issues))
    if mess_handoff or homolysis_native:
        # A new/replaced worker or late-created artifact must not disappear
        # from an otherwise successful capture of the old file inventory.
        progress('rechecking_inventory')
        current=sources(root,tasks,evidence,mess_handoff=mess_handoff,homolysis_native=homolysis_native,
                        ordinary_native=ordinary_native,
                        progress=lambda **d:progress('rechecking_inventory',**d))
        if {str(p) for p,_,_ in current}!={str(p) for p,_,_ in selected}:
            issues.append(dict(source=str(root),error='Evidence inventory changed during capture.'))
        if mess_handoff and handoff_directories()!=directories:
            issues.append(dict(source=str(root),error='MESS handoff directory inventory changed during capture.'))
    manifest=dict(schema=1,root=str(root),files=entries,issues=issues,
                  directories=[str(Path('run')/p.relative_to(root)) for p in sorted(directories)],
                  status='captured' if not issues else 'incomplete',captured_bytes=used,
                  purpose=('ordinary_native_audit' if ordinary_native else
                           'homolysis_native_audit' if homolysis_native else
                           'mess_handoff_replay' if mess_handoff else 'seed_and_restart_audit_only'),
                  native_checkpoint_backup=False,
                  consistent_live_snapshot=False,submits_jobs=False,modifies_source=False)
    atomic_json(output/'manifest.json',manifest)
    extras=[]
    if homolysis_native:
        progress('summarizing_copies')
        from kinbot.homolysis_evidence import summarize, scheduler, SKIP_DIRS
        records,summary=summarize(output,_scopes(output/'run'))
        atomic_json(output/'homolysis_inventory.json',dict(holds=records,
            note='The legacy holds key contains inconclusive evidence records, not export-blocking decisions.',
            excluded_directories=sorted(SKIP_DIRS),excluded_file_components=['rwf','wfu']))
        (output/'summary.txt').write_text(summary)
        progress('scheduler_accounting', inconclusive_checks=len(records))
        atomic_json(output/'scheduler.json',scheduler(output,records))
        extras=['homolysis_inventory.json','summary.txt','scheduler.json']
    archive=output/'seed-evidence.tar.gz'
    partial=archive.with_name(archive.name+'.partial')
    progress('compressing_archive', added=0, total=len(entries), compression_level=1)
    # Level one is already used by the native orbital exporter. Checksums,
    # source rechecks and complete file contents are unchanged by compression.
    with tarfile.open(partial,'w:gz',compresslevel=1) as tar:
        tar.add(output/'manifest.json',arcname='manifest.json',recursive=False)
        for name in extras:tar.add(output/name,arcname=name,recursive=False)
        for name in manifest['directories']:
            path=output/name;path.mkdir(parents=True,exist_ok=True)
            tar.add(path,arcname=name,recursive=False)
        for index,item in enumerate(entries,1):
            tar.add(output/item['path'],arcname=item['path'],recursive=False)
            progress('compressing_archive', added=index, total=len(entries))
    with partial.open('rb') as stream:os.fsync(stream.fileno())
    expected={item['path']:item['sha256'] for item in entries}
    for name in ['manifest.json',*extras]:expected[name]=file_sha256(output/name)
    progress('verifying_archive', checked=0, total=len(expected))
    _verify_archive(partial,expected,progress)
    progress('hashing_archive', archive_bytes=partial.stat().st_size)
    checksum=file_sha256(partial)
    partial.rename(archive)
    result=dict(directory=str(output),archive=str(archive),sha256=checksum,
                files=len(entries),bytes=used,issues=len(issues),status=manifest['status'])
    if homolysis_native:
        result.update(summary=str(output/'summary.txt'),inconclusive_checks=len(records))
    atomic_json(output/'result.json',result)
    progress('finished', **result)
    return result


def detach(a):
    """Launch only the read-only collector, independent of the SSH terminal."""
    root,parent,evidence = _capture_paths(a.root,a.output_parent,a.evidence,a.max_mb)
    job=parent/('seed-capture-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    job.mkdir(parents=True)
    options=dict(root=str(root),output_parent=str(job),evidence=str(evidence) if evidence else None,
                 tasks=a.task,max_mb=a.max_mb,mess_handoff=a.mess_handoff,
                 homolysis_native=a.homolysis_native,ordinary_native=a.ordinary_native)
    atomic_json(job/'request.json',options)
    atomic_json(job/'status.json',dict(state='starting',updated_utc=datetime.now(timezone.utc).isoformat()))
    args=[sys.executable,'-u','-B','-m','kinbot.run_seed_audit',str(root),
          '--output-parent',str(job),'--max-mb',str(a.max_mb),'--job-directory',str(job)]
    for option in ('mess_handoff','homolysis_native','ordinary_native'):
        if getattr(a,option):args.append('--'+option.replace('_','-'))
    if evidence:args.extend(['--evidence',str(evidence)])
    for task in a.task:args.extend(['--task',task])
    env=dict(os.environ)
    package_root=str(Path(__file__).resolve().parent.parent)
    env['PYTHONPATH']=os.pathsep.join(filter(None,[package_root,env.get('PYTHONPATH')]))
    try:
        with (job/'collector.log').open('xb') as log:
            log.write(b'Starting detached read-only evidence capture. No QC jobs or controllers are launched.\n')
            log.flush()
            process=subprocess.Popen(args,cwd=package_root,env=env,stdin=subprocess.DEVNULL,
                                     stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                                     close_fds=True)
        # Separate from status: a quick child can finish before Popen returns.
        launch=dict(job_directory=str(job),pid=process.pid,host=socket.gethostname(),
                    log=str(job/'collector.log'),status_file=str(job/'status.json'),
                    launched_utc=datetime.now(timezone.utc).isoformat(),submits_qc_jobs=False)
    except Exception as error:
        atomic_json(job/'status.json',dict(state='launch_failed',error=str(error)))
        raise
    atomic_json(job/'launch.json',launch)
    return launch


def show_status(job):
    job=Path(job).resolve(strict=True)
    # Also accept the export parent so a new login needs no shell variables.
    if not (job/'status.json').is_file():
        jobs=list(job.glob('seed-capture-*/status.json'))
        if not jobs:raise ValueError('No detached capture status in '+str(job))
        job=max(jobs,key=lambda p:(p.parent/'request.json').stat().st_mtime_ns).parent
    print('CAPTURE_JOB '+str(job))
    print((job/'status.json').read_text())
    log=job/'collector.log'
    if log.is_file():
        with log.open('rb') as stream:
            stream.seek(max(0,log.stat().st_size-16000))
            print('\n'.join(stream.read().decode(errors='replace').splitlines()[-12:]))
    print('Status is saved evidence; a stale running record does not prove a live collector.')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path,nargs='?')
    parser.add_argument('--output-parent',type=Path)
    parser.add_argument('--evidence',type=Path)
    parser.add_argument('--task',action='append',default=[])
    parser.add_argument('--max-mb',type=int,default=2048)
    parser.add_argument('--mess-handoff',action='store_true',
                        help='Include complete MESS artifacts, summaries and PES assembly metadata.')
    parser.add_argument('--homolysis-native',action='store_true',
                        help='Include full ordinary QC and archived attempts for workers with unresolved homolysis checks.')
    parser.add_argument('--ordinary-native',action='store_true',
                        help='Include ordinary native QC for all PES/reverse workers, including HIR, conformers and MESS artifacts.')
    parser.add_argument('--detach',action='store_true',
                        help='Run the read-only collector in its own session; save progress and result after logout.')
    parser.add_argument('--status',type=Path,
                        help='Read status and log tail for a detached capture directory, or the latest under an export parent.')
    parser.add_argument('--job-directory',type=Path,help=argparse.SUPPRESS)
    a=parser.parse_args(argv)
    if a.status:
        if a.root or a.detach or a.job_directory:parser.error('--status is a read-only standalone command.')
        show_status(a.status)
        return
    if a.root is None or a.output_parent is None:parser.error('root and --output-parent are required for capture.')
    if a.detach:
        if a.job_directory:parser.error('--detach cannot be nested.')
        print(json.dumps(detach(a),indent=2),flush=True)
        return
    job=a.job_directory
    state={}
    callback=None
    if job:
        _capture_paths(a.root,job,a.evidence,a.max_mb)
        if job.resolve()!=a.output_parent.resolve():
            parser.error('Detached job directory must equal its output parent.')
        # New session plus redirected descriptors detaches the terminal. Ignore
        # SIGHUP too; this does not promise survival of SIGKILL or a host reboot.
        signal.signal(signal.SIGHUP,signal.SIG_IGN)
        state=dict(state='running',pid=os.getpid(),host=socket.gethostname(),
                   session_id=os.getsid(0),started_utc=datetime.now(timezone.utc).isoformat())
        atomic_json(job/'status.json',state)
        def callback(record):
            state.update(record)
            atomic_json(job/'status.json',state)
    if a.mess_handoff:
        print('Capturing MESS handoff: reading source files; no jobs or source changes.',flush=True)
    if a.homolysis_native or a.ordinary_native:
        print('Capturing native evidence: read-only; no retries or controller changes.',flush=True)
    try:
        result=capture(a.root,a.output_parent,tasks=a.task,evidence=a.evidence,
                       max_mb=a.max_mb,mess_handoff=a.mess_handoff,homolysis_native=a.homolysis_native,
                       ordinary_native=a.ordinary_native,progress_callback=callback)
        if job:
            atomic_json(job/'result.json',result)
            state.update(state='complete' if result['status']=='captured' else 'incomplete',
                         result=result,finished_utc=datetime.now(timezone.utc).isoformat())
            atomic_json(job/'status.json',state)
    except BaseException as error:
        if job:
            state.update(state='failed',error=f'{type(error).__name__}: {error}',
                         finished_utc=datetime.now(timezone.utc).isoformat())
            atomic_json(job/'status.json',state)
        traceback.print_exc()
        raise SystemExit(1)
    print(json.dumps(result,indent=2),flush=True)
    if (a.mess_handoff or a.homolysis_native or a.ordinary_native) and result['status']!='captured':
        raise SystemExit(1)


if __name__=='__main__':main()
