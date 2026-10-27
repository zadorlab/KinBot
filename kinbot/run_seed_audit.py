"""Export seed/restart evidence without changing or driving the source run.

Capture shallow PES/search ledgers, full SQLite database files and their sidecars,
graph state, and explicitly selected interrupted tasks. This is an audit bundle,
not a complete run backup or a declaration that any result is reusable.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tarfile
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


def sources(root, tasks=(), evidence=None):
    paths = set()
    for name in ('prepared.json', '.kinbot_resource_policy.json', '.kinbot_live_control.json'):
        if (root/name).is_file():paths.add(root/name)
    for scope in _scopes(root):
        for name in METADATA | {scope.name+'.json', 'kinbot.db', 'kinbot.db-wal', 'kinbot.db-shm'}:
            if (scope/name).is_file():paths.add(scope/name)
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


def capture(root, output_parent, *, tasks=(), evidence=None, max_mb=2048):
    root=Path(root).resolve(strict=True);output_parent=Path(output_parent).resolve()
    evidence=Path(evidence).resolve(strict=True) if evidence else None
    if output_parent==root or output_parent.is_relative_to(root):
        raise ValueError('Write audit evidence outside the source run.')
    if evidence and (output_parent==evidence or output_parent.is_relative_to(evidence)):
        raise ValueError('Write this audit outside the supplied evidence tree.')
    if type(max_mb) is not int or max_mb<1:raise ValueError('max_mb must be positive.')
    selected=sources(root,tasks,evidence)
    if not selected:raise ValueError('No recognized seed evidence in source run.')
    output=output_parent/('seed-audit-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8])
    output.mkdir(parents=True)
    entries=[];issues=[];used=0
    for source,base,relative in selected:
        try:
            data,info=read_file(source,base,min(256*MIB,max_mb*MIB-used))
            if info['changed']:raise RuntimeError('Source changed during capture.')
            target=output/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
            entries.append(dict(source=str(source),path=str(relative),**info))
            used+=len(data)
        except (OSError,ValueError,RuntimeError) as error:
            issues.append(dict(source=str(source),error=str(error)))
    # Includes DB/WAL sidecars: a change anywhere invalidates the captured
    # group. An actively changing run must be quiesced before seed import.
    for entry in entries:
        source=Path(entry['source'])
        try:
            if source.stat().st_mtime_ns!=entry['mtime_ns'] or file_sha256(source)!=entry['sha256']:
                issues.append(dict(source=str(source),error='Source changed before capture completed.'))
        except OSError as error:issues.append(dict(source=str(source),error=str(error)))
    manifest=dict(schema=1,root=str(root),files=entries,issues=issues,
                  status='captured' if not issues else 'incomplete',captured_bytes=used,
                  purpose='seed_and_restart_audit_only',native_checkpoint_backup=False,
                  consistent_live_snapshot=False,submits_jobs=False,modifies_source=False)
    atomic_json(output/'manifest.json',manifest)
    archive=output/'seed-evidence.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:
        tar.add(output/'manifest.json',arcname='manifest.json',recursive=False)
        for item in entries:tar.add(output/item['path'],arcname=item['path'],recursive=False)
    return dict(directory=str(output),archive=str(archive),sha256=file_sha256(archive),
                files=len(entries),bytes=used,issues=len(issues),status=manifest['status'])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--output-parent',type=Path,required=True)
    parser.add_argument('--evidence',type=Path)
    parser.add_argument('--task',action='append',default=[])
    parser.add_argument('--max-mb',type=int,default=2048)
    a=parser.parse_args(argv)
    print(json.dumps(capture(a.root,a.output_parent,tasks=a.task,evidence=a.evidence,max_mb=a.max_mb),indent=2))


if __name__=='__main__':main()
