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


def contract(inp):
    """Read generated literal metadata without importing/executing an input."""
    inp = Path(inp).resolve()
    tree = ast.parse(inp.read_text())
    meta = next((ast.literal_eval(n.value) for n in tree.body
                 if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name)
                 and t.id == 'validation_metadata' for t in n.targets)), {})
    correction = meta.get('correction', {})
    if correction.get('source'):
        if file_sha256(correction['source']) != correction.get('source_sha256'):
            raise RuntimeError(f'{inp.stem}: correction changed since input generation.')
        from kinbot.vrc_scan_recovery import require_ready_for_rates
        require_ready_for_rates(correction['source'])
    return dict(input_sha256=file_sha256(inp), correction=correction,
                review_required=meta.get('name') == 'production')


def validate(inp, record):
    current = contract(inp)
    if any(record.get(k) != v for k, v in current.items()):
        raise RuntimeError(f'{Path(inp).stem}: campaign scientific contract changed; '
                           'restore it or use a fresh rotdpy_directory.')
    return current


def active(inp, record):
    from kinbot.pes_phase import process_identity
    identity = record.get('owner_identity')
    if identity and process_identity(identity['pid']) == identity:
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


def launch(inp, *, retry=False):
    from kinbot.pes_phase import process_identity
    inp = Path(inp).resolve()
    with locked(inp):
        path = receipt(inp)
        old = json.loads(path.read_text()) if path.exists() else None
        if old:
            validate(inp, old)
            if active(inp, old) or not retry:
                return old
        current = contract(inp)
        record = dict(old or {}, schema=1, input=str(inp), **current)
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


def require_review(inp):
    path = receipt(inp)
    if not path.exists():
        return
    record = json.loads(path.read_text())
    current = validate(inp, record)
    if current['review_required'] and record.get('accepted_input_sha256') != current['input_sha256']:
        raise ReviewRequired(f'{Path(inp).stem}: completed correction requires review before MESS. '
                           'Use kinbot.rotdpy_campaign accept after inspecting its native evidence.')


def accept(inp, note):
    inp = Path(inp).resolve()
    with locked(inp):
        path = receipt(inp)
        record = json.loads(path.read_text())
        current = validate(inp, record)
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
