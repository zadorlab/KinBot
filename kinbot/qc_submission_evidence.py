"""Allocation completion evidence for ordinary KinBot QC submissions.

Releasing a scheduler slot is independent of validating the native result. An
unsuccessful saddle search releases its slot just as a successful minimum does.
The caller still requires a healthy queue snapshot and an explicit purged job
ID reply; missing files or an unconfirmed sbatch never release capacity.
"""
import hashlib
import json
from pathlib import Path
import re
import shlex

from kinbot.execution import file_sha256
from kinbot.slurm_retry import completion_evidence


def terminal_attempt(root, claim):
    root = Path(root).resolve()
    directory = Path(claim.get('directory', '/')).resolve()
    job = str(claim.get('job_id', ''))
    command = claim.get('command', [])
    if (not job.isdecimal() or claim.get('status') != 'submitted'
            or not directory.is_relative_to(root) or len(command) < 2
            or command[0] != 'sbatch'
            or any(x not in ('--parsable', '--export=ALL') for x in command[1:-1])):
        return None
    receipt = root/'.kinbot_policy_submissions'/(job+'.json')
    if json.loads(receipt.read_text()) != claim:
        return None
    script = Path(command[-1])
    script = (directory/script).resolve()
    if not script.is_relative_to(directory) or script.suffix != '.sbatch':
        return None
    # Parse only the generated wrapper's arguments, never evaluate shell code.
    # Timeout-recovery scripts may have a different filename from the QC job.
    guards = []
    for line in script.read_text().splitlines():
        if 'kinbot.slurm_retry' not in line or not line.lstrip().startswith('exec '):
            continue
        args = shlex.split(line)
        if args[2:5] == ['-m', 'kinbot.slurm_retry', 'run']:
            guards.append(args[5:args.index('--')])
    if len(guards) != 1:
        return None
    args = guards[0]
    def option(name):
        return args[args.index(name)+1] if args.count(name) == 1 else None
    identity = option('--identity')
    if (not identity or identity != option('--qc-job')
            or option('--guard-token') != hashlib.sha256(identity.encode()).hexdigest()[:20]
            or option('--exit-marker') != identity+'.exitcode'):
        return None
    marker = (directory/(identity+'.exitcode')).resolve()
    if not marker.is_relative_to(directory):
        return None
    proof = completion_evidence(identity, job, directory)
    if proof:
        source = Path(proof['path'])
        status = proof['status']
    else:
        # Old wrappers and a missing final journal can still have an atomic,
        # allocation-tagged exit marker. A bare legacy exit code is insufficient.
        fields = marker.read_text().split()
        if (len(fields) != 2 or fields[0] != job
                or not re.fullmatch(r'-?[0-9]+', fields[1])):
            return None
        source, status = marker, 'complete' if int(fields[1]) == 0 else 'failed'
        # If an allocation journal exists but is running/malformed, do not
        # fall back to a marker left by an earlier incarnation of a requeue.
        folder = directory/'.kinbot_slurm_io_retries'/option('--guard-token')/job
        if folder.exists():
            return None
    return dict(current={}, receipt=receipt, execution=source,
                execution_status=status)
