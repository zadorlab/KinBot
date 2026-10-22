"""Bounded read-only native failure evidence, including archived attempts.

Usage: python -m kinbot.native_failure_report TASK_DIRECTORY --label STEP
No scheduler actions, calculations, or file writes are performed.
"""
import argparse
import json
from pathlib import Path
import re


ERROR = re.compile(r'\?|\b(?:error|fatal|failed|abort|fehler)\b|'
                   r'not converg|no convergence|insufficient|not enough|cannot|'
                   r'killed by signal|out.of.memory', re.I)


def excerpt(path, *, limit=16*1024*1024):
    """Keep the tail and earlier errors in a bounded window, with offsets."""
    with Path(path).open('rb') as stream:
        stream.seek(0, 2)
        size = stream.tell()
        start = max(0, size-limit)
        stream.seek(start)
        data = stream.read(limit)
    lines = data.decode(errors='replace').splitlines()
    indices = set()
    for i, line in enumerate(lines):
        if ERROR.search(line):
            indices.update(range(max(0, i-2), min(len(lines), i+4)))
    selected = sorted(indices)
    return dict(size_bytes=size, scan_start_byte=start,
                scan_truncated=start > 0, error_context_truncated=len(selected) > 120,
                error_context='\n'.join(lines[i][:1000] for i in selected[-120:]),
                tail=data[-65536:].decode(errors='replace'))


def snapshot(task, label):
    task = Path(task).resolve(strict=True)
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', label):
        raise ValueError('Label must be a basename.')
    graph = task.parent.parent
    if task.parent.name != 'tasks' or not (graph/'workflow.json').is_file():
        raise ValueError('Expected a dispatcher task directory.')
    def obj(path):
        return json.loads(path.read_text()) if path.is_file() else {}
    result = dict(task=str(task), current=obj(graph/'state.json').get('tasks', {}).get(task.name),
                  resources=obj(task/'task.json').get('task', {}).get('resources'), attempts=[])
    # An explicit retry can move the files while the operator is reading them.
    # Report both active and archived attempts rather than interpreting a
    # missing active file as missing evidence or a new failure.
    receipts = sorted(graph.glob('attempts/**/'+label+'.failure.json'))
    receipts += list(task.glob(label+'.failure.json'))
    result['attempt_inventory_truncated'] = len(receipts) > 8
    for receipt in receipts[-8:]:
        if not receipt.resolve().is_relative_to(graph):
            continue
        try:
            failure = obj(receipt)
            entry = dict(path=str(receipt), execution=obj(receipt.parent/'execution.json'),
                         host=failure.get('host'), returncode=failure.get('returncode'),
                         capture_errors=failure.get('capture_errors'), files={})
            entry['execution'] = {k:entry['execution'].get(k) for k in
                                  ('status','slurm_job_id','error')}
            for suffix in ('.out','.log','.stdout','.stderr'):
                path = receipt.parent/(label+suffix)
                if path.is_file() and path.resolve().is_relative_to(graph):
                    evidence = excerpt(path)
                    evidence['tail'] = evidence['tail'][-2200:]
                    entry['files'][path.name] = evidence
            # Rank-local errors may survive only in the scratch receipt even
            # when the main output file is available. Keep distinct evidence.
            fingerprints = {(e.get('tail',''),e.get('error_context',''))
                            for e in entry['files'].values()}
            for name, evidence in failure.get('files', {}).items():
                saved = {**evidence, 'tail': evidence.get('tail','')[-2200:]}
                fingerprint = (saved['tail'], saved.get('error_context',''))
                if fingerprint not in fingerprints:
                    entry['files'][name] = saved
                    fingerprints.add(fingerprint)
            result['attempts'].append(entry)
        except (OSError, ValueError) as error:
            result['attempts'].append(dict(path=str(receipt), read_error=str(error)))
    execution = obj(task/'execution.json')
    result['current_execution'] = {k:execution.get(k) for k in ('status','slurm_job_id','error')}
    log = task/'optimization.log'
    if log.is_file():
        result['current_optimizer'] = excerpt(log, limit=2200)['tail']
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task', type=Path)
    parser.add_argument('--label', required=True)
    args = parser.parse_args()
    print(json.dumps(snapshot(args.task, args.label), indent=2))


if __name__ == '__main__':
    main()
