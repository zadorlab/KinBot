"""Shared evidence and budget rules for allocation-level walltime recovery.

Accounting is optional: retained scontrol records, or a job-specific Slurm
TIME LIMIT diagnostic paired with an observed allocation, provide evidence.
Unknown exits, cancellations and native convergence failures are not timeouts.
"""
import json
from datetime import datetime
import math
from pathlib import Path
import re
import time

from kinbot.execution import atomic_json, slurm_job_record
from kinbot.run_policy import read as read_policy
from kinbot.run_walltime import duration, partition_limits, seconds

_OBSERVED = {}


def policy(directory):
    return read_policy(directory) or {}


def _matching(record, directory, job):
    return (record.get('JobId') == str(job)
            and bool(record.get('WorkDir'))
            and Path(record['WorkDir']).resolve() == Path(directory).resolve())


def observe(directory, job):
    """Save a small allocation record once; no native files or accounting reads."""
    directory = Path(directory).resolve()
    job = str(job)
    if not job.isdecimal() or not policy(directory).get('timeout_retries', 0):
        return
    path = directory / '.kinbot_timeout_allocations' / (job + '.json')
    key = str(path)
    if path.exists() or time.monotonic() - _OBSERVED.get(key, -math.inf) < 60:
        return
    _OBSERVED[key] = time.monotonic()
    record = slurm_job_record(job)
    if (_matching(record, directory, job) and record.get('TimeLimit') and record.get('Partition')
            and record.get('JobState', '').split('+')[0] == 'RUNNING'):
        atomic_json(path, record)


def evidence(directory, job, stderr_files):
    """Caller must independently confirm the exact allocation is no longer live."""
    directory, job = Path(directory).resolve(), str(job)
    if not job.isdecimal():
        return None
    current = slurm_job_record(job)
    if current:
        if not _matching(current, directory, job):
            return None
        if current.get('JobState', '').split('+')[0] != 'TIMEOUT':
            return None  # A positive CANCELLED/FAILED record overrides old text.
        source = 'scontrol'
    else:
        path = directory / '.kinbot_timeout_allocations' / (job + '.json')
        if not path.is_file():
            return None
        current = json.loads(path.read_text())
        if not _matching(current, directory, job):
            return None
        marker = re.compile(r'\bJOB\s+' + re.escape(job)
                            + r'\s+ON\s+\S+\s+CANCELLED AT\s+(\S+)\s+DUE TO TIME LIMIT\b', re.I)
        source = None
        elapsed = None
        for path in stderr_files:
            path = Path(path)
            if not path.is_file():
                continue
            with path.open('rb') as stream:
                stream.seek(max(0, path.stat().st_size - 65536))
                match = marker.search(stream.read().decode(errors='replace'))
                if match:
                    # A manual scontrol increase after the observation may
                    # have changed the real limit. Elapsed allocation time
                    # bounds it conservatively when the terminal record has
                    # already expired. Missing timestamps are not evidence.
                    try:
                        elapsed = (datetime.fromisoformat(match.group(1)) -
                                   datetime.fromisoformat(current['StartTime'])).total_seconds()
                    except (KeyError, TypeError, ValueError):
                        continue
                    if elapsed <= 0:
                        continue
                    source = str(path)
                    break
        if source is None:
            return None
    try:
        previous = seconds(current['TimeLimit'])
    except (KeyError, ValueError):
        return None
    if not 0 < previous < math.inf or not current.get('Partition'):
        return None
    if source != 'scontrol':
        previous = max(previous, elapsed)
    return {'job_id': job, 'previous_walltime': current['TimeLimit'],
            'partition': current['Partition'], 'source': source,
            'previous_allocation_seconds': previous}


def next_walltime(settings, proof):
    """Return a strictly larger limit, bounded by the same partition's MaxTime."""
    previous = max(seconds(proof['previous_walltime']), proof.get('previous_allocation_seconds', 0))
    floor = seconds(settings.get('minimum_walltime', '00:00:00'))
    requested = max(previous * settings.get('timeout_walltime_factor', 2.), floor)
    limits = partition_limits(cached=True)
    try:
        cap = min(limits[name] for name in proof['partition'].split(','))
    except KeyError as error:
        raise RuntimeError('Timeout retry cannot verify the original partition time limit.') from error
    result = min(math.ceil(requested / 60) * 60, cap)
    if result <= previous:
        return None
    return duration(result)
