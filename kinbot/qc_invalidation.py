"""Finish an owner-authorized HIR invalidation without reviving old handoffs."""
import os
from pathlib import Path

from ase.db.core import T2000, YEAR


def requested(qc, job):
    """An explicit invalidation is refinement authority, not a QC failure retry."""
    if not hasattr(qc, 'db'):
        return False
    rows = list(qc.db.select(name=job))
    if not rows:
        return False
    marker = rows[-1]
    if marker.data.get('status') != 0 or marker.data.get('reason') != 'HIR restart':
        return False
    return any(row.id == marker.data.get('supersedes_row_id')
               and row.data.get('status') == 'normal' for row in rows[:-1])


def reconcile(qc, job):
    if not hasattr(qc, 'db'):
        return False
    rows = list(qc.db.select(name=job))
    if not rows:
        return False
    marker = rows[-1]
    if marker.data.get('status') != 0 or marker.data.get('reason') != 'HIR restart':
        return False
    boundary = marker.ctime * YEAR + T2000
    # A replacement may already have started or finished. Never retire its
    # files or scheduler identity, even while its database handoff is pending.
    suffixes = ('.log', '.out', '_sella.log', '.chk', '.fchk', '_freq.out',
                '.exitcode', '.pkl', '.result.json', '.optimization_budget.json',
                '.qc_reuse_consumer.json', '.qc_reuse_producer.json',
                '.minimum_reuse.json', '.shared_native.fchk')
    paths = [Path(job + suffix) for suffix in suffixes]
    if any(p.is_file() and p.stat().st_mtime > boundary for p in paths):
        return False
    current = qc.job_ids.get(job)
    retired = marker.data.get('retired_job_id')
    if current:
        if retired is not None and str(current) != str(retired):
            return False
        from kinbot.qc import _slurm_job_active, _pbs_job_active
        active = (_slurm_job_active(current) if qc.queuing == 'slurm' else
                  _pbs_job_active(current) if qc.queuing == 'pbs' else True)
        if active is not False:
            return False
        # Historical markers did not record the retired job ID. Only a
        # pre-invalidation exit receipt can bind that legacy generation.
        if retired is None:
            from kinbot.qc import _worker_exit_status
            exit_record = _worker_exit_status(job + '.exitcode', current)
            if not exit_record or not exit_record[1]:
                return False
    for path in paths:
        if path.is_file():
            archived = Path(f'{path}.restart_{marker.id}')
            if archived.exists():
                # Preserve both observations if interrupted retirement found
                # a second old artifact. Never overwrite historical evidence.
                from kinbot.execution import file_sha256
                archived = Path(str(archived) + '.' + file_sha256(path))
            os.replace(path, archived)
    qc.job_ids.pop(job, None)
    getattr(qc, 'job_generation_tagged', {}).pop(job, None)
    return True
