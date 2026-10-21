"""Finite step-budget extensions; never waive convergence or catch QC errors."""
from pathlib import Path

from kinbot.execution import atomic_json


DEFAULT_STEP_RETRIES = 2


def step_retries(directory=None):
    from kinbot.run_policy import read
    policy = read(directory or Path.cwd()) or {}
    value = policy.get('optimization_step_retries', DEFAULT_STEP_RETRIES)
    if type(value) is not int or not 0 <= value <= 10:
        raise ValueError('optimization_step_retries must be an integer from 0 to 10.')
    return value


def run_optimizer(optimizer, *, fmax, steps, report=None, chunks=None, used=0):
    """Continue the same optimizer/Hessian in bounded chunks from its last point.

    Only a normal False return (step exhaustion) is continued. Calculator,
    constraint and other exceptions propagate without entering another chunk.
    """
    if chunks is None:
        chunks = 1 + step_retries()
    if type(chunks) is not int or chunks < 1 or type(steps) is not int or steps < 1:
        raise ValueError('Optimizer chunks and steps must be positive integers.')
    record = dict(status='running', chunks_used=used, chunk_steps=steps, fmax=fmax)
    for _ in range(chunks):
        if report:
            atomic_json(Path(report), record)
        # Preserve the optimizer object (and its learned Hessian) on extension.
        converged = bool(optimizer.run(fmax=fmax, steps=steps))
        record['chunks_used'] += 1
        record['status'] = 'converged' if converged else 'exhausted'
        if report:
            atomic_json(Path(report), record)
        if converged:
            return True
    return False
