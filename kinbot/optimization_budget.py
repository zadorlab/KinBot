"""Finite step-budget extensions; never waive convergence or catch QC errors."""
from pathlib import Path

from kinbot.execution import atomic_json


DEFAULT_STEP_RETRIES = 2


def step_retries(directory=None):
    import os
    from kinbot.run_policy import read
    policy = read(directory or Path.cwd()) or {}
    value = policy.get('optimization_step_retries', DEFAULT_STEP_RETRIES)
    if 'KINBOT_OPTIMIZATION_STEP_RETRIES' in os.environ:
        value = int(os.environ['KINBOT_OPTIMIZATION_STEP_RETRIES'])
    if type(value) is not int or not 0 <= value <= 10:
        raise ValueError('optimization_step_retries must be an integer from 0 to 10.')
    return value


def run_optimizer(optimizer, *, fmax, steps, report=None, chunks=None, used=0,
                  remaining_steps=None, used_steps=0):
    """Continue the same optimizer/Hessian in bounded chunks from its last point.

    Only a normal False return (step exhaustion) is continued. Calculator,
    constraint and other exceptions propagate without entering another chunk.
    """
    if chunks is None:
        chunks = 1 + step_retries()
    if type(chunks) is not int or chunks < 1 or type(steps) is not int or steps < 1:
        raise ValueError('Optimizer chunks and steps must be positive integers.')
    if remaining_steps is not None and (type(remaining_steps) is not int or remaining_steps < 0):
        raise ValueError('Remaining optimizer steps must be a nonnegative integer.')
    if type(used_steps) is not int or used_steps < 0:
        raise ValueError('Previously used optimizer steps must be a nonnegative integer.')
    record = dict(status='running', chunks_used=used, chunk_steps=steps, fmax=fmax)
    for _ in range(chunks):
        if report:
            atomic_json(Path(report), record)
        # Preserve the optimizer object (and its learned Hessian) on extension.
        allowed=steps if remaining_steps is None else min(steps,remaining_steps)
        before=getattr(optimizer,'nsteps',0)
        converged = bool(optimizer.run(fmax=fmax, steps=allowed))
        if remaining_steps is not None:
            consumed=getattr(optimizer,'nsteps',0)-before
            if not 0 <= consumed <= allowed:
                raise RuntimeError('Optimizer reported an invalid step count.')
            remaining_steps-=consumed
            used_steps+=consumed
            record['remaining_steps']=remaining_steps
        record['chunks_used'] += 1
        if remaining_steps is not None:
            record['steps_used']=used_steps
            record['chunks_used']=max(used,(used_steps+steps-1)//steps)
        record['status'] = 'converged' if converged else 'exhausted'
        if report:
            atomic_json(Path(report), record)
        if converged:
            return True
        if remaining_steps == 0:
            break
    return False
