"""Allocation identity for shared-pool orbital reconnaissance attempts."""
import json
from pathlib import Path

from kinbot.execution import file_sha256


def terminal_attempt(root, claim):
    """Prove a terminal native attempt; queue absence alone proves nothing.

    The caller additionally requires an explicit purged-ID scheduler reply.
    A matching failed execution releases capacity, not scientific acceptance.
    """
    directory = Path(claim.get('directory', '/')).resolve()
    ident = str(claim.get('job_id', ''))
    if (not ident.isdecimal() or claim.get('status') != 'submitted'
            or not directory.is_relative_to(root)):
        return None
    receipt = root/'.kinbot_policy_submissions'/(ident+'.json')
    if json.loads(receipt.read_text()) != claim:
        return None
    manifest = json.loads((directory.parent/'reconnaissance.json').read_text())
    if manifest.get('shared_policy_file') != str(root/'.kinbot_resource_policy.json'):
        return None
    matches = [j for j in manifest['jobs'] if Path(j['directory']).resolve() == directory]
    if len(matches) != 1:
        return None
    job = matches[0]
    if claim.get('command') != ['sbatch', '--export=ALL', '--parsable', job['name']+'.sbatch']:
        return None
    attempt, proof = directory, job
    if str(job.get('job_id')) != ident:
        history = [a for a in job.get('attempt_history', []) if str(a.get('job_id')) == ident]
        if len(history) != 1:
            return None
        proof = history[0]
        attempt = Path(proof['archive']).resolve()
        if not attempt.is_relative_to(directory/'attempts'):
            return None
    execution_path = attempt/'execution.json'
    execution = json.loads(execution_path.read_text())
    if (str(execution.get('slurm_job_id')) != ident
            or execution.get('status') not in ('failed', 'complete')
            or execution.get('input_sha256') != proof['input_sha256']
            or file_sha256(attempt/(job['name']+'.inp')) != proof['input_sha256']):
        return None
    return dict(current=job, receipt=receipt, execution=execution_path,
                execution_status=execution['status'])
