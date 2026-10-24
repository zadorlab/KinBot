"""Independent, immutable ROTD campaigns with live retry and review receipts.

Ready PES channels may sample before other workers finish. Final aggregation
adopts these exact inputs/results; it never launches a second sampling owner.
"""
import argparse
import ast
from contextlib import contextmanager
import fcntl
import json
import logging
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from kinbot.execution import atomic_json, file_sha256

logger = logging.getLogger('KinBot')


class ReviewRequired(RuntimeError):
    pass


def receipt(inp):
    return Path(inp).with_suffix('.campaign.json')


@contextmanager
def locked(inp):
    with Path(inp).with_suffix('.campaign.lock').open('a+') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        yield


def contract(inp, *, require_ready=True):
    """Read immutable metadata; optionally enforce current scan readiness.

    Ownership survives a scan hold. New sampling and result use still require
    the current complete, reviewed scan and its correction binding.
    """
    inp = Path(inp).resolve()
    tree = ast.parse(inp.read_text())
    meta = next((ast.literal_eval(n.value) for n in tree.body
                 if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name)
                 and t.id == 'validation_metadata' for t in n.targets)), {})
    correction = meta.get('correction', {})
    if correction.get('source'):
        if file_sha256(correction['source']) != correction.get('source_sha256'):
            raise RuntimeError(f'{inp.stem}: correction changed since input generation.')
        if require_ready:
            from kinbot.vrc_scan_recovery import require_ready_for_rates
            require_ready_for_rates(correction['source'])
    return dict(input_sha256=file_sha256(inp), correction=correction,
                review_required=meta.get('name') == 'production')


def validate(inp, record, *, require_ready=True):
    current = contract(inp, require_ready=require_ready)
    if any(record.get(k) != v for k, v in current.items()):
        raise RuntimeError(f'{Path(inp).stem}: campaign scientific contract changed; '
                           'restore it or use a fresh rotdpy_directory.')
    return current


def retained_processes(root, processes):
    """Identify independent campaign owners and their current descendants.

    Refresh may retain only a hash-bound campaign with the exact live process
    identity. A module name or a directory alone is not ownership evidence.
    Cleanup deliberately does not use this filter.
    """
    from kinbot.pes_phase import process_identity, same_process_identity
    root = Path(root).resolve()
    retained = set()
    for process in processes:
        args = shlex.split(process['command'])
        if '-m' not in args:
            continue
        suffix = args[args.index('-m') + 1:]
        if len(suffix) != 3 or suffix[:2] != ['kinbot.rotdpy', 'run']:
            continue
        inp = (Path(process['cwd']) / suffix[2]).resolve()
        if not inp.is_relative_to(root):
            continue
        saved = json.loads(receipt(inp).read_text())
        identity = saved.get('owner_identity')
        current_identity = process_identity(process['pid'])
        if (saved.get('input') != str(inp) or inp.parent != Path(process['cwd']).resolve()
                or not identity or identity.get('pid') != process['pid']
                or identity.get('start') != process['starttime']
                or not same_process_identity(identity, current_identity)):
            raise RuntimeError(f'{inp}: independent campaign process identity changed. '
                               f'Saved={identity!r}; current={current_identity!r}; '
                               f'inventory_start={process["starttime"]!r}; '
                               f'cwd={process["cwd"]!r}; saved_input={saved.get("input")!r}.')
        validate(inp, saved, require_ready=False)
        retained.add(process['pid'])
    while True:
        children = {p['pid'] for p in processes if p.get('ppid') in retained}
        if children <= retained:
            return retained
        retained.update(children)


def replacement_processes(root, processes):
    retained = retained_processes(root, processes)
    return [p for p in processes if p['pid'] not in retained]


def owns_job(root, job):
    """Match a ROTD surface allocation to its campaign and pool receipt."""
    from kinbot import run_policy
    root = Path(root).resolve()
    directory = Path(job['directory']).resolve()
    if not directory.is_relative_to(root):
        return False
    for folder in directory.parents:
        if not folder.is_relative_to(root):
            break
        if not folder.name.startswith('kb_'):
            continue
        inp = folder.with_name(folder.name[3:] + '.py')
        if not receipt(inp).is_file():
            continue
        saved = json.loads(receipt(inp).read_text())
        if saved.get('input') != str(inp):
            raise RuntimeError(f'{inp}: campaign input identity changed.')
        validate(inp, saved, require_ready=False)
        policy = run_policy.find(directory)
        if policy is None:
            return False
        path = policy.parent/'.kinbot_policy_submissions'/f'{job["job_id"]}.json'
        if not path.is_file():
            return False
        claim = json.loads(path.read_text())
        command = claim.get('command', [])
        return (claim.get('job_id') == job['job_id']
                and claim.get('status') == 'submitted'
                and claim.get('directory') == str(directory)
                and len(command) == 2 and command[0] == 'sbatch'
                and (directory/command[-1]).is_file())
    return False


def active(inp, record):
    from kinbot.pes_phase import process_identity, same_process_identity
    identity = record.get('owner_identity')
    if identity and same_process_identity(identity, process_identity(identity['pid'])):
        return True
    # Also detect an explicitly requested `rotdpy run` retry, which uses the
    # same executor lock but has no obligation to rewrite a launch receipt.
    path = Path(inp).with_suffix('.controller.lock')
    if not path.exists():
        return False
    with path.open('r') as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
    return False


def launch(inp, *, retry=False, defer_scan_review=False):
    from kinbot.pes_phase import process_identity
    from kinbot.vrc_scan_recovery import ScanReviewRequired
    inp = Path(inp).resolve()
    with locked(inp):
        path = receipt(inp)
        old = json.loads(path.read_text()) if path.exists() else None
        if old:
            validate(inp, old, require_ready=False)
            if old.get('operator_rejected'):
                return old
            if active(inp, old):
                if old.pop('resume_requested', None):
                    old.pop('resume_hold', None)
                    atomic_json(path, old)
                return old
            if not retry and not old.get('operator_reinclude') and not old.get('resume_requested'):
                return old
        try:
            current = contract(inp)
        except ScanReviewRequired as error:
            if not defer_scan_review or old is None:
                raise
            # The replacement VRC owner calls launch again when this same
            # bound correction is ready. Persist one request, not a retry loop.
            old.update(resume_requested=True, resume_hold=error.channels)
            atomic_json(path, old)
            return old
        record = dict(old or {}, schema=1, input=str(inp), **current)
        record.pop('operator_reinclude', None)
        record.pop('resume_requested', None)
        record.pop('resume_hold', None)
        # Save the contract before spawning, including on a launch failure.
        # An interrupted handoff cannot silently create another generation.
        record.update(status='launching', launches=record.get('launches', 0)+1)
        record.pop('error', None)
        atomic_json(path, record)
        try:
            with inp.with_suffix('.campaign.log').open('ab') as log:
                process = subprocess.Popen([sys.executable, '-u', '-m', 'kinbot.rotdpy',
                                            'run', str(inp)], cwd=inp.parent,
                                           stdout=log, stderr=subprocess.STDOUT,
                                           start_new_session=True)
        except OSError as error:
            record.update(status='failed', error=str(error), owner_identity=None)
            atomic_json(path, record)
            raise
        record.update(status='launched', owner_identity=process_identity(process.pid))
        atomic_json(path, record)
        return record


def resume_retained(root):
    """Adopt each bound campaign independently during a controller handoff.

    Scan holds are reported and never waived. Integrity/ownership errors stay
    fatal; only an absent owner of an unchanged, ready input may be launched.
    """
    from kinbot import rotdpy
    from kinbot.vrc_scan_recovery import ScanReviewRequired
    root = Path(root).resolve()
    outcomes = []
    for path in sorted(root.glob('*/rotdPy/*.campaign.json')):
        saved = json.loads(path.read_text())
        inp = Path(saved['input']).resolve(strict=True)
        if not inp.is_relative_to(root) or receipt(inp) != path:
            raise RuntimeError(f'Campaign ownership mismatch: {path}')
        validate(inp, saved, require_ready=False)
        row = dict(input=str(inp))
        if saved.get('operator_rejected'):
            row['status'] = 'excluded'
        elif active(inp, saved):
            row['status'] = 'owner_retained'
        else:
            execution = rotdpy.execution_path(inp)
            complete = (execution.exists()
                        and json.loads(execution.read_text()).get('status') == 'complete')
            try:
                if complete:
                    validate(inp, saved)
                    rotdpy.read_result(inp)
                    row['status'] = 'complete_retained'
                else:
                    result = launch(inp, retry=True, defer_scan_review=True)
                    if result.get('operator_rejected'):
                        row['status'] = 'excluded'
                    elif result.get('resume_requested'):
                        row.update(status='held', channels=result['resume_hold'], resume_requested=True)
                    else:
                        row.update(status='owner_resumed', owner_identity=result.get('owner_identity'))
            except ScanReviewRequired as error:
                row.update(status='held', channels=error.channels, resume_requested=False)
        outcomes.append(row)
    return outcomes


def require_review(inp):
    """Verify scan readiness and immutable sampling input, without a second vote.

    contract() checks the relaxed-scan evidence and correction binding. The
    legacy acceptance receipt remains audit history, not a new MESS gate.
    """
    path = receipt(inp)
    if not path.exists():
        return
    record = json.loads(path.read_text())
    if record.get('operator_rejected'):
        raise RuntimeError(f'{Path(inp).stem}: channel was rejected by the operator.')
    validate(inp, record)


def stop_for_rejection(inp, archive, request):
    """Stop exactly one bound sampling owner and its verified allocations."""
    import signal
    from kinbot import run_refresh, run_walltime, run_policy
    from kinbot.pes_phase import process_identity, same_process_identity
    from kinbot.anl.live_control import stop_attempt
    inp = Path(inp).resolve()
    with locked(inp):
        path = receipt(inp)
        record = json.loads(path.read_text())
        validate(inp, record, require_ready=False)
        saved = Path(archive)/'campaigns'/path.name
        if not saved.exists():
            atomic_json(saved, record)
        record['operator_rejected'] = dict(request=request['id'], reason=request['reason'])
        atomic_json(path, record)  # Fence restart before signalling the owner.
        if active(inp, record):
            identity = record.get('owner_identity')
            if not identity or not same_process_identity(identity, process_identity(identity['pid'])):
                raise ValueError('Sampling owner identity unavailable; no process signalled.')
            try:
                fd = run_refresh._pidfd_open(identity['pid'])
            except ProcessLookupError:
                return False
            try:
                if not same_process_identity(identity, process_identity(identity['pid'])):
                    return False
                run_refresh._pidfd_send_signal(fd, signal.SIGTERM)
            finally:
                os.close(fd)
            return False  # Re-enter after its forwarding/cleanup handlers exit.
        policy = run_policy.find(inp.parent)
        if policy is None:
            raise ValueError('Sampling cancellation requires the shared allocation ledger.')
        sample_root = inp.with_name('kb_'+inp.stem)
        pending = False
        for job in run_walltime.jobs(policy.parent):
            directory = Path(job['directory']).resolve()
            if not directory.is_relative_to(sample_root):
                continue
            claim = dict(job_id=job['id'], directory=str(directory))
            if not owns_job(policy.parent, claim):
                raise ValueError('Sampling allocation ownership is unverified; retained.')
            cancel_path = saved.parent/(job['id']+'.cancel.json')
            cancel = (json.loads(cancel_path.read_text()) if cancel_path.exists()
                      else dict(request, kind='reaction_decision', cancel_requested=False))
            pending |= not stop_attempt(directory, job['id'], cancel_path, cancel, allow_completed=True)
        return not pending


def accept(inp, note, *, expected=None):
    inp = Path(inp).resolve()
    with locked(inp):
        path = receipt(inp)
        record = json.loads(path.read_text())
        current = validate(inp, record)
        if expected is not None and current != expected:
            raise RuntimeError('Campaign changed after review; inspect the current evidence.')
        record.update(accepted_input_sha256=current['input_sha256'], review_note=note)
        atomic_json(path, record)
    return record


def ensure(inp, *, interval=15):
    """Join an early campaign, keeping a failed campaign open for live retry."""
    from kinbot import live_control, rotdpy
    inp = Path(inp).resolve()
    launch(inp)
    last_error = None
    while True:
        from kinbot.pes_phase import checkpoint
        checkpoint()
        record = json.loads(receipt(inp).read_text())
        validate(inp, record)
        execution = rotdpy.execution_path(inp)
        result = json.loads(execution.read_text()) if execution.exists() else {}
        if active(inp, record):
            time.sleep(interval)
            continue
        try:
            if result.get('status') != 'complete':
                raise RuntimeError(result.get('error') or record.get('error') or
                    f'{inp.stem}: sampling owner exited without a complete result; live retry required.')
            rotdpy.read_result(inp)
            require_review(inp)
            return result
        except RuntimeError as error:
            config = live_control.read(inp.parent) or {}
            if not (isinstance(error, ReviewRequired) or config.get('hold_failures')
                    or config.get('hold_task_failures')):
                raise
            if str(error) != last_error:
                logger.warning('%s; independent campaigns continue.', error)
                last_error = str(error)
            time.sleep(interval)


def managed_worker(parameters):
    root = os.environ.get('KINBOT_PES_PHASE_ROOT')
    return bool(parameters.get('pes') and parameters.get('rotdpy_run') and root
                and Path(root).resolve() == Path.cwd().parent.resolve())


def ready(vts, names):
    """Publish fully validated channels from a managed PES worker immediately."""
    if not managed_worker(vts.par):
        return
    from kinbot.species_routing import routing_name
    par = dict(vts.par, rotdpy_run=0)
    rows = [[routing_name(vts.well), n,
             [routing_name(p) for p in vts.scan_reac[n].products], 0.] for n in names]
    _stage(par, rows)


def ready_leaf(vts, channel, state):
    """An independent physical surface need not wait for its sibling roots."""
    if not managed_worker(vts.par):
        return
    from kinbot.species_routing import routing_name
    correction = Path(state['correction'])
    if file_sha256(correction) != state['correction_sha256']:
        raise RuntimeError('Completed surface correction changed before sampling.')
    leaf = json.loads(correction.read_text())
    par = dict(vts.par)
    par.update(leaf['surface_parameters'])
    reactant = routing_name(vts.well)
    par.update(rotdpy_run=0, vrc_tst_surface_model=None,
               vrc_tst_surface_model_file='', vrc_tst_noscan={},
               vrc_tst_scan={reactant: [state['reaction']]})
    _stage(par, [[reactant, state['reaction'],
                 [routing_name(p) for p in vts.scan_reac[channel].products], 0.]])


def _stage(par, rows):
    root = os.environ.get('KINBOT_PES_PHASE_ROOT')
    if not root or Path(root).resolve() != Path.cwd().parent.resolve():
        raise RuntimeError('Early ROTD sampling requires a managed PES worker.')
    from kinbot.pes import create_rotdpy_inputs
    folder = Path(par.get('rotdpy_directory', 'rotdPy'))
    par['rotdpy_directory'] = str(folder if folder.is_absolute() else Path(root)/folder)
    Path(par['rotdpy_directory']).mkdir(parents=True, exist_ok=True)
    # Global fragment/surface rendering and shared scheduler files have one
    # producer per case at a time. Only the caller thread enters this path.
    with locked(Path(par['rotdpy_directory'])/'staging'):
        files = create_rotdpy_inputs(par, rows, [], correction_root=Path.cwd())
        for inp in files:
            launch(inp)
    logger.info('Ready VRC surface(s) sampling independently: %s', ', '.join(r[1] for r in rows))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['status', 'retry', 'accept'])
    p.add_argument('input', type=Path)
    p.add_argument('--note', default='')
    a = p.parse_args(argv)
    if a.action == 'accept':
        if not a.note.strip():
            p.error('accept requires --note describing the correction review')
        result = accept(a.input, a.note)
    elif a.action == 'retry':
        result = launch(a.input, retry=True)
    else:
        result = json.loads(receipt(a.input).read_text())
        result['active'] = active(a.input, result)
        from kinbot.rotdpy import execution_path
        path = execution_path(a.input)
        result['execution'] = json.loads(path.read_text()) if path.exists() else None
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
