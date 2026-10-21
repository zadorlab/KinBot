"""Generation-fenced timeout retries for ordinary L1/L2 QC workers.

No chemistry or optimizer input is regenerated here. A clean retry executes
the same Python worker with a longer Slurm budget, through shared admission.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import threading
import time

from kinbot import run_policy, timeout_retry
from kinbot.execution import atomic_json, file_sha256

_LOCK = threading.RLock()
_GRACE = 60.
logger = logging.getLogger('KinBot')


def _marker(job):
    token = hashlib.sha256(str(job).encode()).hexdigest()[:24]
    return Path.cwd()/'.kinbot_timeout_retries'/token/'current.json'


@contextmanager
def _locked(marker):
    marker.parent.mkdir(parents=True, exist_ok=True)
    with _LOCK, (marker.parent/'lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def _belongs(record, job_id):
    return str(job_id) in [record['original_job_id'], *record.get('replacement_ids', [])]


def _inputs_match(record):
    for name, expected in record['inputs'].items():
        if file_sha256(name) != expected:
            raise RuntimeError(f'Timeout retry input changed: {name}; review before submitting.')


def _finish(qc, job, marker, record):
    """Resume archival/submission; never issue sbatch twice after uncertainty."""
    from kinbot.anl.dispatch import SchedulerUnavailableError
    from kinbot.qc import _slurm_job_active

    _inputs_match(record)
    attempt = record['attempts'][-1]
    old_job = attempt['job_id']
    if record['phase'] == 'preparing':
        if _slurm_job_active(old_job) is not False:
            return 'running'
        archive = Path(attempt['archive'])
        for source, target in record['moves'].items():
            source, target = Path(source), Path(target)
            if source.exists():
                if target.exists():
                    raise RuntimeError(f'Timeout artifact reappeared: {source}')
                target.parent.mkdir(parents=True, exist_ok=True)
                source.rename(target)
        record['phase'] = 'restoring'
        atomic_json(marker, record)
    if record['phase'] == 'restoring':
        archive = Path(attempt['archive'])
        for target, seed in record['seeds'].items():
            saved = Path(seed['file'])
            if file_sha256(saved) != seed['sha256']:
                raise RuntimeError('Timeout retry pre-allocation seed changed.')
            Path(target).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(saved, target)
        # Leave old successful rows as evidence, but do not let them satisfy
        # the new generation if that worker exits without a new result.
        rows = list(qc.db.select(name=job))
        if rows and rows[-1].data.get('status') != 0:
            qc.db.write(rows[-1].toatoms(), name=job,
                        data={'status': 0, 'timeout_retry_of': old_job})
        original = (archive/'batch.original').read_text()
        fallback=record.get('fallback', {})
        if record.get('numerical_recovery') or record.get('resource_overrides') or fallback:
            from kinbot.manual_fallback import worker_options, worker_script
            from kinbot.qc_live import worker_resources, script_resources
            worker=Path(record['effective_worker'])
            revised=worker.read_text()
            if record.get('numerical_recovery'):
                from kinbot.qc_optimizer_recovery import render
                revised=render(revised,record['numerical_recovery'],job)
            if record.get('resource_overrides'):
                revised=worker_resources(revised,record['backend'],record['resource_overrides'])
                original=script_resources(original,record['resource_overrides'])
            if fallback:
                revised=worker_options(revised,fallback)
            replacement_worker=archive/'retry_worker.py'
            replacement_worker.write_text(revised)
            record['inputs'][str(replacement_worker)]=file_sha256(replacement_worker)
            original=worker_script(original,record['worker_command'],replacement_worker,fallback)
        # Slurm stops reading directives at the first executable line; the
        # new final directive must precede that line and override old -t.
        lines = original.splitlines(keepends=True)
        index = 0
        while index < len(lines) and (not lines[index].strip() or lines[index].lstrip().startswith('#')):
            index += 1
        lines.insert(index, '#SBATCH --time=' + attempt['walltime'] + '\n')
        lines.insert(index + 1, '#SBATCH --partition=' + attempt['partition'] + '\n')
        lines.insert(index + 2, 'export KINBOT_SLURM_IO_RETRY_SEED=' + record['random_seed'] + '\n')
        retry_script = Path(record['script'])
        retry_script.write_text(''.join(lines))
        record['script_sha256'] = file_sha256(retry_script)
        record['phase'] = 'ready'
        atomic_json(marker, record)
    if record['phase'] == 'submitting':
        # A controller can die after sbatch succeeds but before it saves the
        # returned ID. The pool's durable accepted receipt closes this gap.
        policy_path = run_policy.find()
        if policy_path is None:
            raise RuntimeError('Timeout submission receipt cannot be checked without its shared policy.')
        matches = []
        for path in (policy_path.parent/'.kinbot_policy_submissions').glob('*.json'):
            receipt = json.loads(path.read_text())
            if (receipt.get('directory') == str(Path.cwd())
                    and receipt.get('command') == ['sbatch', record['script']]):
                matches.append(str(receipt['job_id']))
        if len(matches) != 1:
            raise RuntimeError('Timeout retry submission outcome uncertain; reconcile Slurm before resubmitting.')
        record['replacement_ids'].append(matches[0])
        record['phase'] = 'submitted'
        atomic_json(marker, record)
    if record['phase'] == 'ready':
        if file_sha256(record['script']) != record['script_sha256']:
            raise RuntimeError('Timeout retry Slurm script changed.')
        record['phase'] = 'submitting'
        atomic_json(marker, record)
        try:
            env = dict(os.environ, SBATCH_TIMELIMIT=attempt['walltime'],
                       SBATCH_PARTITION=attempt['partition'])
            response = run_policy.submit(['sbatch', record['script']], cwd=Path.cwd(),
                                         capture_output=True, text=True, timeout=30, env=env)
        except SchedulerUnavailableError:
            # This error precedes submission. All other exceptions, including
            # an OSError while saving an accepted receipt, remain uncertain.
            record['phase'] = 'ready'
            atomic_json(marker, record)
            return 'running'
        except OSError as error:
            if not getattr(error, 'kinbot_submission_not_started', False):
                raise
            record['phase'] = 'ready'
            atomic_json(marker, record)
            return 'running'
        if response is None or response.returncode:
            record['phase'] = 'ready'
            if response is not None:
                record['submission_error'] = response.stderr
            atomic_json(marker, record)
            if response is not None:
                raise RuntimeError(f'Timeout retry sbatch rejected: {response.stderr}')
            return 'running'
        match = re.fullmatch(r'Submitted batch job\s+(\d+)\s*', response.stdout)
        if not match:
            raise RuntimeError('Timeout retry returned an unrecognized Slurm ID; reconcile submission.')
        record['replacement_ids'].append(match.group(1))
        record['phase'] = 'submitted'
        atomic_json(marker, record)
        logger.info('%s: timeout job %s replaced by %s; clean retry with walltime %s; archive %s.',
                    job, old_job, match.group(1), attempt['walltime'], attempt['archive'])
    if record['phase'] == 'submitted':
        qc.job_ids[job] = record['replacement_ids'][-1]
        qc.job_generation_tagged[job] = True
        if record.get('manual_request_file'):
            from kinbot import live_control
            path = Path(record['manual_request_file'])
            request = json.loads(path.read_text())
            if request['id'] != record['manual_request']:
                raise RuntimeError('Manual retry receipt identity changed.')
            live_control.acknowledge(path, request, 'applied', job_id=qc.job_ids[job],
                                     archive=record['attempts'][-1]['archive'])
    return None


def resume(qc, job):
    marker = _marker(job)
    if not marker.exists() or job not in qc.job_ids:
        return None
    with _locked(marker):
        record = json.loads(marker.read_text())
        if not _belongs(record, qc.job_ids[job]):
            return None  # Same logical name, later constrained-path step.
        if (record['phase'] == 'submitted'
                and str(qc.job_ids[job]) == record['replacement_ids'][-1]):
            return None
        return _finish(qc, job, marker, record)


def recover(qc, job, *, manual=None):
    from kinbot.qc import _slurm_job_active
    from kinbot.slurm_retry import _artifacts

    settings = timeout_retry.policy(Path.cwd())
    limit = settings.get('timeout_retries', 0)
    if job not in qc.job_ids:
        return None
    job_id = str(qc.job_ids[job])
    if manual is not None:
        from kinbot.anl.live_control import _ownership
        native = _ownership(Path.cwd(), job_id)
        proof = {'job_id': job_id, 'partition': native['Partition'],
                 'previous_walltime': native['TimeLimit'], 'source': 'operator',
                 'request': manual['id']}
    else:
        proof = (timeout_retry.evidence(Path.cwd(), job_id, [Path('perm')/(job+'.err')])
                 if limit else None)
        if proof is None:
            from kinbot.qc_scf_recovery import evidence
            proof = evidence(qc, job, job_id)
        if proof is None:
            from kinbot.qc_optimizer_recovery import evidence
            proof = evidence(qc, job, job_id)
    if proof is None:
        return None
    if _slurm_job_active(job_id) is not False:
        return 'running'  # Includes unknown scheduler status for an eligible retry.
    waits = getattr(qc, '_timeout_publication_waits', {})
    qc._timeout_publication_waits = waits
    now = time.time()
    if manual is None and now - waits.setdefault(job_id, now) < _GRACE:
        return 'running'
    marker = _marker(job)
    with _locked(marker):
        record = json.loads(marker.read_text()) if marker.exists() else None
        if record and _belongs(record, job_id):
            if manual and record.get('manual_request') == manual['id']:
                return _finish(qc, job, marker, record) or 'running'
            if record['phase'] != 'submitted':
                return _finish(qc, job, marker, record) or 'running'
            history = record['attempts']
        else:
            record, history = None, []
        numerical_sources=('diatomic_symmetry','optimizer_steps','optimizer_ode')
        symmetry_retry = proof['source'] in numerical_sources
        consumed = sum(a.get('source') == proof['source'] for a in history) if symmetry_retry else sum(
            a.get('source') not in ('operator', *numerical_sources) for a in history)
        if manual is None and consumed >= (1 if symmetry_retry else limit):
            qc._timeout_retry_blocked = f'{job}: timeout retry limit reached ({limit}).'
            return None
        from kinbot import live_control
        overrides = live_control.overrides(Path.cwd()/job) if manual else {}
        if manual and manual.get('fallback'):
            overrides.update(manual['fallback'].get('resources',{}))
        walltime = (overrides.get('walltime', proof['previous_walltime']) if manual or symmetry_retry
                    else timeout_retry.next_walltime(settings, proof))
        if walltime is None:
            qc._timeout_retry_blocked = f'{job}: partition maximum prevents a longer allocation.'
            return None
        script, worker = Path(job+'.sbatch').resolve(), Path(job+'.py').resolve()
        if not all(p.is_file() and p.is_relative_to(Path.cwd()) for p in (script, worker)):
            raise RuntimeError('Timeout retry requires the original owned worker and Slurm script.')
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive = marker.parent/stamp
        archive.mkdir()
        if record:
            atomic_json(archive/'previous.json', record)
        previous_script = Path(record['script']) if record and job_id in record['replacement_ids'] else script
        shutil.copy2(previous_script, archive/'batch.original')
        shutil.copy2(worker, archive/worker.name)
        inputs = {str(p): file_sha256(p) for p in (script, worker)}
        effective_worker = worker
        worker_command = job+'.py'
        if record and job_id in record['replacement_ids']:
            previous_worker = Path(record['attempts'][-1]['archive'])/'retry_worker.py'
            if previous_worker.exists():
                effective_worker, worker_command = previous_worker, str(previous_worker)
            elif record.get('effective_worker'):
                effective_worker = Path(record['effective_worker'])
                worker_command = record['worker_command']
            inputs[str(previous_script)] = file_sha256(previous_script)
            inputs[str(effective_worker)] = file_sha256(effective_worker)
        # The existing I/O guard saved the inputs/checkpoints before this
        # allocation started. Restore those verified seeds, not a checkpoint
        # possibly interrupted during a write. Required Guess=Read inputs in
        # IRC/constrained workers therefore survive a clean retry too.
        token = hashlib.sha256(job.encode()).hexdigest()[:20]
        baseline = Path.cwd()/'.kinbot_slurm_io_retries'/token/job_id
        guard_file = baseline/'execution.json'
        if (not guard_file.is_file() and manual and manual.get('pending_only')
                and manual.get('pending_holds', {}).get(job_id, {}).get('phase') == 'cancelled'):
            import secrets
            guard = {'identity': job, 'slurm_job_id': job_id, 'initial_artifacts': {},
                     'random_seed': str(secrets.randbelow(2**31-1))}
            # The pending allocation never started. Preserve its current seeds
            # before archiving, including Guess=Read checkpoints and inputs.
            for source in _artifacts([job, 'perm/'+job], [], marker.parent):
                relative = str(source.relative_to(Path.cwd()))
                saved = baseline/'initial'/relative
                saved.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(source,saved)
                guard['initial_artifacts'][relative] = file_sha256(saved)
            atomic_json(guard_file,guard)
        elif not guard_file.is_file():
            raise RuntimeError(f'{job}: pre-allocation retry seeds are unavailable; manual review required.')
        else:
            guard = json.loads(guard_file.read_text())
        if guard.get('identity') != job or str(guard.get('slurm_job_id')) != job_id:
            raise RuntimeError('Timeout retry seed generation does not match the failed job.')
        seeds = {}
        for relative, digest in guard['initial_artifacts'].items():
            target, saved = (Path.cwd()/relative).resolve(), baseline/'initial'/relative
            if (not target.is_relative_to(Path.cwd()) or not saved.resolve().is_relative_to(baseline)
                    or file_sha256(saved) != digest):
                raise RuntimeError('Timeout retry pre-allocation seed is invalid.')
            if str(target) in inputs:
                if inputs[str(target)] != digest:
                    raise RuntimeError('Timeout worker changed after its failed allocation began.')
                continue
            seeds[str(target)] = {'file': str(saved), 'sha256': digest}
        seed = str(guard.get('random_seed', ''))
        if not seed.isdecimal():
            raise RuntimeError('Timeout retry initial random seed is unavailable.')
        # Archive only the current job's owned artifacts, never sibling
        # calculations, shared database rows or other saved attempts.
        artifacts = _artifacts([job, 'perm/'+job], [], marker.parent)
        moves = {str(p): str(archive/'files'/p.relative_to(Path.cwd()))
                 for p in artifacts if str(p) not in inputs}
        attempt = {**proof, 'walltime': walltime, 'archive': str(archive), 'mode': 'clean_retry'}
        record = {'phase': 'preparing', 'job': job,
                  'original_job_id': record['original_job_id'] if record else job_id,
                  'replacement_ids': record['replacement_ids'] if record else [],
                  'attempts': [*history, attempt], 'inputs': inputs, 'moves': moves,
                  'script': str(archive/'retry.sbatch'), 'seeds': seeds, 'random_seed': seed}
        record.update(effective_worker=str(effective_worker),worker_command=worker_command)
        if manual:
            record.update(manual_request=manual['id'], resource_overrides=overrides,
                          manual_request_file=manual['request_file'],
                          backend=qc.qc, original_worker=str(worker), fallback=manual.get('fallback',{}))
        if proof['source'] in ('optimizer_steps','optimizer_ode'):
            record['numerical_recovery']=proof
            record['effective_worker']=str(effective_worker)
            record['worker_command']=worker_command
        atomic_json(marker, record)
        _finish(qc, job, marker, record)
        return 'running'
