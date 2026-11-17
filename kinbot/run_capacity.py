"""Size a one-node allocation pool from unique, eligible Slurm nodes.

This changes only the shared admission ceiling. Slurm remains responsible for
availability, other users, reservations and account limits. Existing running
allocations and every other policy setting are preserved.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess

from kinbot import run_policy


def node_inventory(text, excluded, cores, memory_mb):
    nodes = {}
    for line in text.splitlines():
        if not line.strip(): continue
        fields = line.strip().split('|')
        if (len(fields) != 5 or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', fields[0])
                or not fields[2].isdigit() or not fields[3].isdigit()
                or fields[4] not in ('up', 'down', 'inact', 'drain')):
            raise ValueError('Malformed per-node Slurm inventory.')
        name, state, cpu, memory, available = fields
        record = dict(state=state.lower(), cores=int(cpu), memory_mb=int(memory), partition=available)
        if name in nodes and nodes[name] != record:
            raise ValueError('Conflicting node records; repeat the Slurm inventory.')
        nodes[name] = record
    if not nodes: raise ValueError('Empty Slurm node inventory.')
    eligible, idle, rejected = [], [], {}
    for name, row in sorted(nodes.items()):
        if name in excluded: reason = 'excluded'
        elif row['partition'] != 'up': reason = 'partition_unavailable'
        elif row['state'] not in ('idle', 'allocated', 'mixed', 'completing'): reason = row['state']
        elif row['cores'] < cores or row['memory_mb'] < memory_mb: reason = 'resources'
        else:
            eligible.append(name)
            if row['state'] == 'idle': idle.append(name)
            continue
        rejected[name] = reason
    return dict(eligible=eligible, idle=idle, rejected=rejected, unique_nodes=len(nodes))


def fit(directory, *, partition, cores, memory_mb, apply=False):
    if (type(cores) is not int or cores < 1 or type(memory_mb) is not int or memory_mb < 1
            or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', partition)):
        raise ValueError('Use one named partition and positive per-job resource requirements.')
    path = run_policy.find(directory)
    if path is None: raise ValueError('Enable a shared policy before fitting its capacity.')
    previous = run_policy._read(path)
    query = subprocess.run(['sinfo', '-N', '-h', '-p', partition, '-o', '%N|%T|%c|%m|%a'],
                           text=True, capture_output=True, check=True, timeout=30)
    report = node_inventory(query.stdout, previous['excluded_nodes'], cores, memory_mb)
    ceiling = len(report['eligible'])
    if not ceiling: raise ValueError('No healthy, nonexcluded nodes fit these jobs; policy unchanged.')
    report.update(root=str(path.parent), partition=partition, cores=cores, memory_mb=memory_mb,
        previous_max_jobs=previous['max_jobs'], max_jobs=ceiling, applied=False,
        excluded_nodes=previous['excluded_nodes'])
    if apply and previous['max_jobs'] != ceiling:
        result = run_policy.update(path.parent, max_jobs=ceiling, expected_policy=previous)
        report.update(applied=True, revision=result['revision'])
    elif run_policy._read(path) != previous:
        raise RuntimeError('Shared policy changed during capacity review; repeat the inventory.')
    report['effective_max_jobs'] = run_policy.read(path.parent)['max_jobs']
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--partition', required=True)
    parser.add_argument('--cores', required=True, type=int)
    parser.add_argument('--memory-mb', required=True, type=int)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    print(json.dumps(fit(args.root, partition=args.partition, cores=args.cores,
        memory_mb=args.memory_mb, apply=args.apply), indent=2))


if __name__ == '__main__': main()
