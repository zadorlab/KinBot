"""Adopt valid orphaned QC and archive incompatible geometry attempts.

The old pathname remains a compatibility alias. It must remain until this
run is retired: submitted scripts and accepted provenance contain that name.
"""
import argparse
from contextlib import ExitStack
import ctypes
import errno
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time

from kinbot.execution import atomic_json, file_sha256


def _exchange(old, new):
    """Linux atomic directory/symlink exchange; never leave a missing pathname."""
    library = ctypes.CDLL(None, use_errno=True)
    function = getattr(library, 'renameat2', None)
    if function is None:
        raise OSError(errno.ENOSYS, 'renameat2(RENAME_EXCHANGE) is unavailable.')
    function.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
                         ctypes.c_char_p, ctypes.c_uint]
    function.restype = ctypes.c_int
    if function(-100, os.fsencode(old), -100, os.fsencode(new), 2):
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), str(old))


def rehome(old, new):
    """Keep existing directory inodes and all old absolute paths usable."""
    old, new = Path(old).absolute(), Path(new).absolute()
    if old.parent != new.parent or old == new:
        raise ValueError('Source and destination must be distinct sibling directories.')
    if old.is_symlink():
        if old.resolve() != new or not new.is_dir():
            raise ValueError('Source is an unrelated symlink.')
        return new
    if not old.is_dir():
        raise ValueError('Source run does not exist.')
    if new.is_symlink() and new.resolve() == old:
        return old  # Safe compatibility view on filesystems without exchange.
    placeholder = (new.is_symlink() and os.readlink(new) == str(new)
                   and (old / 'salvage.json').exists())
    if (new.exists() or new.is_symlink()) and not placeholder:
        raise ValueError('Destination already exists; refusing to merge runs.')
    # Before exchange this is an intentionally self-referential placeholder;
    # afterwards the link occupies old and resolves to the moved directory.
    if not placeholder:
        new.symlink_to(new, target_is_directory=True)
    try:
        _exchange(old, new)
    except BaseException as error:
        if new.is_symlink() and os.readlink(new) == str(new):
            new.unlink()
        if isinstance(error, OSError) and error.errno in (errno.ENOSYS, errno.EINVAL, errno.EOPNOTSUPP):
            # Some shared filesystems cannot atomically exchange directory and
            # symlink entries. Leave every live path/inode intact; NEW is an
            # alternate entry point to the same physical run, never a copy.
            new.symlink_to(old, target_is_directory=True)
            return old
        raise
    return new


def _cases(root, updates=None):
    record = json.loads((root / 'prepared.json').read_text())
    # Resolve only after rehome: old paths remain valid aliases.
    paths = [Path(p).resolve() for p in record['inputs']]
    if not paths or any(not p.is_relative_to(root) for p in paths):
        raise ValueError('Prepared inputs are outside the requested run.')
    for name, digest in record['sha256'].items():
        allowed = {digest}
        if updates and name in updates:
            allowed.update(updates[name].values())
        if file_sha256(root / name) not in allowed:
            raise ValueError(f'Prepared input changed: {name}')
    return record, paths


def _input_update(path, values, root, archive, history):
    """Journal the exact allowed bytes before changing a prepared input."""
    _save(path, root, archive)
    name = str(path.relative_to(root))
    digest = hashlib.sha256((json.dumps(values, indent=2, sort_keys=True)+'\n').encode()).hexdigest()
    changes = history.setdefault('input_updates', {})
    if name in changes and changes[name]['after'] != digest:
        raise ValueError(f'Interrupted migration settings changed: {name}')
    changes.setdefault(name, {'before': file_sha256(path), 'after': digest})
    atomic_json(root / 'salvage.json', history)
    atomic_json(path, values)


def _inventory(root):
    from kinbot.run_cleanup import _jobs, _processes
    jobs = _jobs(root)  # scheduler outage is a hard stop, not an empty queue
    processes = _processes(root)
    if processes:
        raise RuntimeError('Controllers/processes are still alive; no files or jobs changed: '
                           + ', '.join(str(p['pid']) for p in processes))
    for marker_name in ('*/launch.json', '*/pes_status.json', '*/[0-9]*/kinbot_status.json'):
        for marker in root.glob(marker_name):
            host = json.loads(marker.read_text()).get('host')
            if host and host.split('.')[0] != socket.gethostname().split('.')[0]:
                raise RuntimeError(f'{marker}: run salvage on original controller host {host}.')
    owned = {}
    for path in root.glob('*/[0-9]*/anl_automatic/states/*/composite/state.json'):
        state = json.loads(path.read_text())
        for ident, entry in state['tasks'].items():
            if entry['status'] == 'submitting':
                raise RuntimeError(f'{path}/{ident}: ambiguous submission.')
            if entry['status'] == 'submitted':
                owned[str(entry['job_id'])] = (path.parent / 'tasks' / ident).resolve()
    for job in jobs:
        if owned.get(job['job_id']) != Path(job['directory']).resolve():
            raise RuntimeError(f'Live job {job["job_id"]} lacks matching dispatcher ownership.')
    return {'jobs': jobs, 'processes': processes}


def _cancel_incompatible(root, history):
    """Cancel only owned, live tasks whose geometry provenance is obsolete."""
    from kinbot.anl.salvage import geometry_invalidations
    from kinbot.anl import dispatch
    from kinbot.run_cleanup import _cancel
    targets = {}
    for path in root.glob('*/[0-9]*/anl_automatic/states/*/composite/workflow.json'):
        graph, spec, state = dispatch._load(path.parent)
        by_id = {task['id']: task for task in spec['tasks']}
        for ident in geometry_invalidations(spec, state):
            entry = state['tasks'].get(ident, {})
            if entry.get('status') == 'submitted':
                dispatch._verify_stage_files(graph, by_id[ident], entry)
                targets[str(entry['job_id'])] = str((path.parent / 'tasks' / ident).resolve())
    inventory = _inventory(root)
    live = {j['job_id']: j for j in inventory['jobs']}
    selected = [live[j] for j in targets if j in live]
    if any(str(Path(j['directory']).resolve()) != targets[j['job_id']] for j in selected):
        raise RuntimeError('Cancellation ownership changed; stopping recovery.')
    history.setdefault('incompatible_jobs', {}).update(targets)
    atomic_json(root / 'salvage.json', history)
    if selected:
        print('Stopping incompatible geometry attempts: ' + ', '.join(j['job_id'] for j in selected), flush=True)
        _cancel(root, selected)
    deadline = time.monotonic() + 120
    while True:
        remaining = [j for j in _inventory(root)['jobs'] if j['job_id'] in targets]
        if not remaining:
            return
        if time.monotonic() >= deadline:
            raise RuntimeError('Cancellation is pending; rerun the same command after these jobs leave Slurm. No task files were moved.')
        time.sleep(2)


def _save(path, root, archive, *, move=False):
    if not path.exists():
        return
    if path.is_symlink() or not path.resolve().is_relative_to(root):
        raise ValueError(f'Refusing to archive an external artifact: {path}')
    target = archive / path.relative_to(root)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        # A completed move is idempotent; never overwrite preserved evidence.
        if move:
            raise ValueError(f'Archive destination already exists: {target}')
        return
    if move:
        path.rename(target)
    elif path.is_dir():
        shutil.copytree(path, target)
    else:
        shutil.copy2(path, target)


def _controller_history(root):
    """Read a completed launch receipt without modifying a fresh campaign.

    The first controller refresh of a fresh run uses the normal handoff
    transaction. Its salvage.json is written after handoff preflight succeeds;
    fresh_run.json remains the immutable record of the original launch.
    """
    marker = root / 'salvage.json'
    if marker.exists():
        return json.loads(marker.read_text())
    fresh = json.loads((root / 'fresh_run.json').read_text())
    if fresh.get('schema') != 1 or fresh.get('phase') != 'launched':
        raise RuntimeError('Controller refresh requires a completed fresh-run launch.')
    _, inputs = _cases(root)
    pin = Path(fresh['code'])
    if (pin.parent != root or not pin.name.startswith('.fresh_code.')
            or pin.is_symlink() or not pin.is_dir()):
        raise ValueError('Fresh launch code must be a local pinned directory.')
    for path in inputs:
        receipt = json.loads((path.parent / 'launch.json').read_text())
        if (receipt != fresh.get('launches', {}).get(path.parent.name)
                or receipt.get('code') != str(pin)
                or receipt.get('input') != str(path)
                or type(receipt.get('pid')) is not int or receipt['pid'] <= 0):
            raise RuntimeError(f'{path.parent}: fresh launch receipt does not match; '
                               'no controllers or allocations changed.')
    return dict(status='launched', inputs=[str(p) for p in inputs], code=str(pin),
                fresh_launch_record=str(root / 'fresh_run.json'))


def _pin_code(root, *, refresh=False):
    import kinbot
    source = Path(kinbot.__file__).parent
    record = _controller_history(root)
    pin = Path(record.get('code', root / '.salvage_code'))
    if refresh:
        hashes = {str(p.relative_to(source)): file_sha256(p)
                  for p in source.rglob('*') if p.is_file()
                  and '__pycache__' not in p.parts and p.suffix != '.pyc'}
        digest = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
        pin = root / ('.salvage_code.' + digest[:16])
    if pin.parent != root or not pin.name.startswith('.salvage_code') or pin.is_symlink():
        raise ValueError('Pinned source must be a local salvage code directory.')
    if not pin.exists():
        temporary = pin.with_name(pin.name + '.preparing')
        if temporary.exists():
            shutil.rmtree(temporary)
        temporary.mkdir()
        shutil.copytree(source, temporary / 'kinbot',
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        hashes = {str(p.relative_to(temporary)): file_sha256(p)
                  for p in (temporary / 'kinbot').rglob('*') if p.is_file()}
        atomic_json(temporary / 'manifest.json', hashes)
        temporary.rename(pin)
    for name, digest in json.loads((pin / 'manifest.json').read_text()).items():
        if file_sha256(pin / name) != digest:
            raise ValueError('Pinned KinBot source changed; refusing a mixed-version restart.')
    return pin


def _worker_route_plans(root):
    """Check prospective worker profiles before changing any route ledger."""
    from kinbot.parameters import Parameters
    from kinbot.profiled_qc import relaxed_scf_routes
    _, paths = _cases(root, json.loads((root / 'salvage.json').read_text()).get('input_updates'))
    plans = []
    for path in paths:
        parameters = Parameters(str(path)).par
        for worker in sorted(path.parent.glob('[0-9]*')):
            if not worker.is_dir():
                continue
            ledger = worker / parameters['profiled_route_file']
            if ledger.is_symlink() or not ledger.resolve().is_relative_to(worker.resolve()):
                raise ValueError(f'{ledger}: theory ledger must remain inside its worker.')
            if ledger.exists():
                payload = json.loads(ledger.read_text())
                revised = relaxed_scf_routes(payload, parameters)
                plans.append((ledger, payload, revised))
    return plans


def _migrate_worker_routes(root, archive):
    # _inventory rejects every live allocation except dispatcher-owned ANL
    # tasks. Those use their own task/spec ledgers and are not edited here.
    _inventory(root)
    plans = _worker_route_plans(root)
    for ledger, previous, revised in plans:
        if previous == revised:
            continue
        _save(ledger, root, archive)
        atomic_json(ledger, revised)


def _preflight_worker_routes(root):
    """Exercise the real QC router used at worker startup, without submitting."""
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    _, paths = _cases(root, _controller_history(root).get('input_updates'))
    previous = Path.cwd()
    try:
        for path in paths:
            parameters = Parameters(str(path)).par
            for worker in sorted(path.parent.glob('[0-9]*')):
                if not worker.is_dir():
                    continue
                os.chdir(worker)
                QuantumChemistry(parameters)
    finally:
        os.chdir(previous)


def resume_failed_controllers(root):
    """Repair the observed startup profile mismatch, leaving QC jobs intact."""
    root = Path(root).resolve(strict=True)
    marker = root / 'salvage.json'
    history = json.loads(marker.read_text())
    _inventory(root)
    restart = history.get('controller_restart')
    if restart and restart['status'] == 'prepared':
        return history
    if restart and restart['status'] == 'launched':
        raise RuntimeError('The reviewed controller restart was already launched; inspect its new status before retrying.')
    if not restart:
        if history['status'] != 'launched':
            raise ValueError('Failed-controller recovery requires a previously launched salvage.')
        _, paths = _cases(root)
        for path in paths:
            status = path.parent / 'pes_status.json'
            record = json.loads(status.read_text()) if status.exists() else {}
            if (record.get('status') != 'failed' or
                    'Theory profiles changed since jobs were submitted' not in record.get('error', '')):
                raise RuntimeError(f'{status}: expected the reviewed theory-profile startup failure.')
        # Verify the scientific changes before archiving any controller state.
        _worker_route_plans(root)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive = Path(history['archive']) / ('controllers.' + stamp)
        archive.mkdir()
        restart = {'status': 'preparing', 'archive': str(archive),
                   'previous_code': history['code']}
        history['controller_restart'] = restart
        atomic_json(marker, history)
    archive = Path(restart['archive'])
    _migrate_worker_routes(root, archive)
    _preflight_worker_routes(root)
    _, paths = _cases(root)
    for path in paths:
        for directory in [path.parent, *sorted(path.parent.glob('[0-9]*'))]:
            if not directory.is_dir():
                continue
            for name in ('driver.log', 'pes.log', 'pes_summary.txt', 'launch.json',
                         'pes_status.json', 'kinbot.log', 'kinbot.err', 'kinbot.out',
                         'kinbot_status.json'):
                _save(directory / name, root, archive, move=True)
    # Keep the previous snapshot in place; a new snapshot serves only the
    # replacement controllers and their subsequently created QC jobs.
    pin = _pin_code(root, refresh=True)
    restart['status'] = 'prepared'
    history.update(status='prepared', code=str(pin))
    atomic_json(marker, history)
    return history


def prepare(old, new, *, pivot_backend='molpro', pivot_method='UKS/B3LYP3',
            pivot_basis='cc-pVDZ', scf=1.e-8, cc=1.e-7, l3_max_retries=1):
    """Repair an idle controller tree while its owned QC allocations continue."""
    from kinbot.anl.salvage import repair_graph
    from kinbot.parameters import Parameters
    old, new = Path(old).absolute(), Path(new).absolute()
    root = new if old.is_symlink() and old.resolve() == new else old
    marker = root / 'salvage.json'
    history = json.loads(marker.read_text()) if marker.exists() else None
    if history and any(history[k] != value for k, value in (
            ('scf', scf), ('cc', cc), ('l3_max_retries', l3_max_retries), ('pivot_backend', pivot_backend),
            ('pivot_method', pivot_method), ('pivot_basis', pivot_basis))):
        raise ValueError('Interrupted salvage options changed.')
    if history and history.get('status') in ('prepared', 'launched'):
        return history
    if type(l3_max_retries) is not int or l3_max_retries < 0:
        raise ValueError('l3_max_retries must be a nonnegative integer.')
    inventory = _inventory(root)
    record, paths = _cases(root, (history or {}).get('input_updates'))
    # The journal travels with the directory; an interrupted exchange remains
    # discoverable from either pathname. Only incompatible owned jobs are stopped.
    if history is None:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive = old.with_name(old.name + '.salvage_archive.' + stamp)
        archive.mkdir()
        history = {'schema': 1, 'status': 'preparing', 'source': str(old),
                   'destination': str(new), 'archive': str(archive),
                   'scf': scf, 'cc': cc, 'l3_max_retries': l3_max_retries, 'pivot_backend': pivot_backend,
                   'pivot_method': pivot_method, 'pivot_basis': pivot_basis,
                   'inventory': inventory, 'graphs': []}
        atomic_json(marker, history)
    archive = Path(history['archive'])
    _cancel_incompatible(root, history)
    root = rehome(old, new)
    marker = root / 'salvage.json'
    # Old accepted QC and scripts retain their exact contents and hashes.
    # Only launch inputs and the future numerical policy are migrated.
    history['physical_root'] = str(root)
    history['layout'] = 'exchanged' if root == new else 'destination_alias'
    atomic_json(marker, history)
    print(f'Recovery view: {new}; physical directory: {root}. Both paths remain usable.', flush=True)
    record, paths = _cases(root, history.get('input_updates'))
    _save(root / 'prepared.json', root, archive)
    for path in paths:
        _save(path, root, archive)
        values = json.loads(path.read_text())
        values.update(anl_scf_convergence=scf, anl_cc_convergence=cc,
                      anl_require_l3_geometry=1, anl_l3_max_retries=l3_max_retries,
                      vrc_tst_scf_convergence=scf,
                      vrc_tst_correlation_convergence=cc,
                      vrc_tst_pivot_backend=pivot_backend)
        if pivot_backend == 'molpro':
            values.update(vrc_tst_pivot_method=pivot_method, vrc_tst_pivot_basis=pivot_basis)
        values['input_file'] = str(path)
        _input_update(path, values, root, archive, history)
        resolved = Parameters(str(path)).par
        from kinbot.qc import _export_slurm_excluded_nodes
        _export_slurm_excluded_nodes(resolved.get('slurm_exclude_nodes', []))
        print(f'Reconciling {path.parent.name} calculations...', flush=True)
        _save(path.parent / 'resolved_parameters.json', root, archive)
        _input_update(path.parent / 'resolved_parameters.json', resolved, root, archive, history)
        # Checkpoint-free fragment orbitals are recomputed by the new Molpro
        # backend. Existing Gaussian DB entries cannot authorize that path.
        for worker in sorted(path.parent.glob('[0-9]*')):
            if not worker.is_dir():
                continue
            plan_file = worker / 'anl_automatic/network_request.json'
            if plan_file.exists():
                _save(plan_file, root, archive)
                plan = json.loads(plan_file.read_text())
                plan['resources'].update(scf_convergence=scf, cc_convergence=cc, require_l3_geometry=True, l3_max_retries=l3_max_retries)
                atomic_json(plan_file, plan)
            # Reference runs that never produced an accepted model are
            # restarted locally; successful references are not shared across
            # PES runs or recomputed.
            for reference in worker.glob('anl_automatic/reference_l2/*'):
                if (reference.is_dir() and not (reference / 'thermochemistry_model.json').exists()):
                    _save(reference, root, archive, move=True)
            for graph_file in sorted(worker.glob('anl_automatic/states/*/composite/workflow.json')):
                graph = graph_file.parent
                outcome = repair_graph(graph, archive / 'graphs' / graph.relative_to(root),
                                       scf=scf, cc=cc, l3_max_retries=l3_max_retries)
                history['graphs'] = [r for r in history['graphs'] if r['graph'] != str(graph)]
                history['graphs'].append({'graph': str(graph), **outcome})
                atomic_json(marker, history)
            # Failed VRC scans cannot be reused as accepted points. Preserve
            # their evidence separately and retain the fragment/native files.
            for target in ('vrctst/molpro_scan', '.kinbot_post_reaction_branches.json'):
                item = worker / target
                if item.exists():
                    if item.is_dir():
                        for state in item.glob('*/state.json'):
                            entries = json.loads(state.read_text())['tasks'].values()
                            if any(e['status'] in ('complete', 'submitted', 'submitting') for e in entries):
                                raise RuntimeError(f'{state}: VRC result/ownership needs explicit review.')
                    _save(item, root, archive, move=True)
        for directory, dirs, files in os.walk(path.parent):
            dirs[:] = [d for d in dirs if d not in ('tasks', 'attempts')
                       and not Path(directory, d).is_symlink()]
            for name in files:
                if name in ('driver.log', 'pes.log', 'kinbot.log', 'kinbot.err',
                            'kinbot.out', 'launch.json', 'kinbot_status.json', 'pes_status.json'):
                    _save(Path(directory, name), root, archive, move=True)
    record['inputs'] = [str(p) for p in paths]
    record['sha256'] = {name: file_sha256(root / name) for name in record['sha256']}
    atomic_json(root / 'prepared.json', record)
    _migrate_worker_routes(root, archive)
    _preflight_worker_routes(root)
    pin = _pin_code(root)
    history.update(status='prepared', code=str(pin), inputs=record['inputs'])
    atomic_json(marker, history)
    return history


def launch(root):
    """Launch only absent controllers; repeated invocations never duplicate them."""
    from kinbot.run_cleanup import _processes
    root = Path(root).resolve(strict=True)
    marker = root / 'salvage.json'
    record = json.loads(marker.read_text())
    if record['status'] == 'launched':
        return record
    if record['status'] != 'prepared':
        raise ValueError('Salvage has not finished preparing.')
    _pin_code(root)
    _preflight_worker_routes(root)
    from kinbot.parameters import Parameters
    from kinbot.anl.automatic import validate_automatic_environment
    for filename in record['inputs']:
        print(f'Checking launch environment for {Path(filename).parent.name}...', flush=True)
        validate_automatic_environment(Parameters(filename).par)
    for filename in record['inputs']:
        path = Path(filename)
        case = path.parent
        launch_file = case / 'launch.json'
        from kinbot.rotdpy_campaign import replacement_processes
        processes = replacement_processes(case, _processes(case))
        if processes:
            if not launch_file.exists():
                raise RuntimeError(f'{case}: unknown live controller; no duplicate launched.')
            continue
        if launch_file.exists():
            # It might have started and failed quickly. Do not repeatedly
            # restart it simply because the overall launch was interrupted.
            continue
        env = dict(os.environ)
        env['PYTHONPATH'] = record['code'] + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
        env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', '')
        # Match fresh launch: only native batch jobs use the configured QC
        # core count. Controller imports must not create node-sized BLAS pools.
        for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            env[key] = '1'
        with (case / 'driver.log').open('ab') as log:
            process = subprocess.Popen([sys.executable, '-u', '-m', 'kinbot.pes', str(path)],
                cwd=case, env=env, stdin=subprocess.DEVNULL, stdout=log,
                stderr=subprocess.STDOUT, start_new_session=True)
        atomic_json(launch_file, {'pid': process.pid, 'host': socket.gethostname(),
                    'input': str(path), 'code': record['code']})
    record['status'] = 'launched'
    if record.get('controller_restart'):
        record['controller_restart']['status'] = 'launched'
    atomic_json(marker, record)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--launch', action='store_true')
    parser.add_argument('--resume-failed', action='store_true',
                        help='Resume controllers stopped by the reviewed theory-profile startup mismatch.')
    parser.add_argument('--pivot-backend', choices=('molpro', 'gaussian'), default='molpro')
    parser.add_argument('--pivot-method', default='UKS/B3LYP3')
    parser.add_argument('--pivot-basis', default='cc-pVDZ')
    parser.add_argument('--scf', type=float, default=1.e-8)
    parser.add_argument('--cc', type=float, default=1.e-7)
    parser.add_argument('--l3-max-retries', type=int, default=1)
    args = parser.parse_args(argv)
    lock_path = args.destination.absolute().with_suffix('.salvage.lock')
    with lock_path.open('a') as lock, ExitStack() as locks:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        source = args.source.resolve(strict=True)
        for path in source.glob('*/[0-9]*/.kinbot_run.json'):
            handle = locks.enter_context(path.open())
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.resume_failed:
            result = resume_failed_controllers(args.destination)
        else:
            result = prepare(args.source, args.destination, pivot_backend=args.pivot_backend,
                             pivot_method=args.pivot_method, pivot_basis=args.pivot_basis,
                             scf=args.scf, cc=args.cc, l3_max_retries=args.l3_max_retries)
        # Child drivers acquire their own worker locks after this process
        # returns; release migration locks before starting them.
        locks.close()
        if args.launch:
            result = launch(args.destination)
        print(json.dumps({k: result[k] for k in ('status', 'destination', 'archive', 'code')}, indent=2))
        for graph in result['graphs']:
            print(Path(graph['graph']).parent.name, 'complete=', len(graph['complete']),
                  'adopted=', graph['adopted_jobs'], 'retry=', graph['retries'])


if __name__ == '__main__':
    main()
