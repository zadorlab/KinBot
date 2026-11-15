"""Stage verified dependency checkpoints into Gaussian's isolated workspace.

The logical deck and workflow stay unchanged. Only the executed copy uses a
local OLDCHK name; its bytes and the imported checkpoint are retained as
execution artifacts. Arbitrary external files are never imported.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

from kinbot.execution import file_sha256

OLDCHK = re.compile(r'(?im)^(\s*%oldchk\s*=\s*)(\S+)(\s*)$')
REASON = 'gaussian_dependency_checkpoint_v1: retry verified pre-execution path rejection'


def plan(directory, input_name, text):
    """Return a rewritten execution copy and hash-bound dependency imports."""
    matches = [m for m in OLDCHK.finditer(text) if Path(m[2]).is_absolute()]
    if not matches:
        return text, []
    from kinbot.anl import dispatch
    directory = Path(directory).resolve()
    if directory.parent.name != 'tasks' or not (directory/'task.json').is_file():
        raise ValueError(f'Native QC file must be a relative basename: {matches[0][2]!r}')
    graph, spec, state = dispatch._load(directory.parent.parent)
    task = next(t for t in spec['tasks'] if t['id'] == directory.name)
    if (task.get('kind') != 'external' or task.get('backend') != 'gaussian'
            or task.get('stdin') != input_name or task['input_name'] != input_name):
        raise ValueError('Checkpoint import requires a declared Gaussian input.')
    dispatch._verify_stage_files(graph, task, state['tasks'][task['id']])
    if (directory/input_name).read_text() != text:
        raise ValueError('Checkpoint input changed during staging.')
    bindings = {}
    for match in matches:
        source = Path(match[2])
        # Require the exact published path, not a symlink or a path outside
        # this graph. The producer must be an explicit dependency.
        if (source != source.resolve() or source.parent.parent != graph/'tasks'
                or source.parent.name not in dispatch._task_dependencies(task)
                or source.suffix.lower() != '.chk'):
            raise ValueError('OLDCHK is not a checkpoint of a declared dependency.')
        parent = next(t for t in spec['tasks'] if t['id'] == source.parent.name)
        entry = state['tasks'].get(parent['id'], {})
        if (entry.get('status') != 'complete' or parent.get('backend') != 'gaussian'
                or source.name not in parent.get('required_outputs', ())):
            raise ValueError('OLDCHK producer has no accepted checkpoint output.')
        dispatch._verify_stage_files(graph, parent, entry)
        execution = source.parent/'execution.json'
        outcome = json.loads(execution.read_text())
        if outcome.get('status') != 'executed':
            raise ValueError('OLDCHK producer has no successful execution.')
        dispatch._verify_execution(graph, parent, entry, outcome)
        digest = outcome['artifacts'][source.name]
        local = 'kinbot_parent_'+hashlib.sha256(str(source).encode()).hexdigest()[:20]+'.chk'
        if local in text or (directory/local).exists():
            # A newly staged retry directory is empty of imported files. Do
            # not silently use a previous workspace's copy.
            raise ValueError('Checkpoint import filename already exists.')
        bindings[str(source)] = dict(source=str(source), source_task=parent['id'],
            sha256=digest, source_execution_sha256=file_sha256(execution), local=local)
    rendered = OLDCHK.sub(lambda m: m[1]+bindings[m[2]]['local']+m[3]
                         if m[2] in bindings else m[0], text)
    return rendered, list(bindings.values())


def retry_requests(root, *, apply=False):
    """Queue one repair through the live owner, never edit its state or jobs."""
    from kinbot.anl import dispatch
    from kinbot import live_control as control
    root = Path(root).resolve(strict=True)
    graphs = [root] if (root/'workflow.json').is_file() else sorted(
        p.parent for p in root.glob('*/workflow.json'))
    reports = []
    for graph in graphs:
        _, spec, state = dispatch._load(graph)
        marker = control.find(graph)
        for task in spec['tasks']:
            ident = task['id']; entry = state['tasks'].get(ident, {})
            if entry.get('status') != 'failed' or task.get('backend') != 'gaussian':
                continue
            directory = graph/'tasks'/ident
            outcome_path = directory/'execution.json'
            if not outcome_path.is_file():
                continue
            outcome = json.loads(outcome_path.read_text())
            text = (directory/task['input_name']).read_text()
            sources = [m[2] for m in OLDCHK.finditer(text) if Path(m[2]).is_absolute()]
            if (outcome.get('status') != 'failed' or outcome.get('failure_phase') != 'validation'
                    or outcome.get('error') not in {
                        f'Native QC file must be a relative basename: {s!r}' for s in sources}):
                continue
            report = dict(target=str(directory), status='eligible'); reports.append(report)
            try:
                dispatch._verify_execution(graph, task, entry, outcome)
                if any((directory/task.get(k, default)).is_file()
                       and (directory/task.get(k, default)).stat().st_size
                       for k, default in (('stdout', 'stdout.txt'), ('stderr', 'stderr.txt'))):
                    raise ValueError('Native output exists; inspect it before retrying.')
                _, bindings = plan(directory, task['input_name'], text)
                report['checkpoints'] = bindings
                if any(r.get('reason') == REASON for r in state.get('manual_retries', {}).get(ident, [])):
                    report['status'] = 'exhausted'; continue
                if marker is None:
                    raise ValueError('Enable live controls at the run root first.')
                previous = [json.loads(p.read_text()) for p in
                            (marker.parent/'.kinbot_controls/requests').glob('*.json')]
                if any(r.get('target') == str(directory) and
                       (r.get('status') in ('queued', 'applying') or
                        (r.get('reason') == REASON and r.get('status') == 'applied')) for r in previous):
                    report['status'] = 'already_requested'; continue
                expected = {k: entry.get(k) for k in ('attempt', 'job_id', 'task_sha256', 'status')}
                if apply:
                    request = control.request(marker.parent, directory, reason=REASON, expected=expected)
                    report.update(status='queued', request=request['id'])
            except (OSError, ValueError, RuntimeError, KeyError) as error:
                report.update(status='blocked', error=str(error))
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path, help='One graph or its immediate containing directory.')
    parser.add_argument('--apply', action='store_true', help='Queue verified repairs for the owning controllers.')
    args = parser.parse_args()
    print(json.dumps(retry_requests(args.root, apply=args.apply), indent=2))


if __name__ == '__main__':
    main()
