"""Stage, preflight and launch the two reviewed examples; never use SSH.

This example helper uses the general KinBot drivers and scoped cleanup API.
Preparation is read-only with respect to old runs. Retirement is explicit.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from kinbot.execution import atomic_json, file_sha256
from kinbot.parameters import Parameters
from kinbot.run_format import ensure_current_run
from install_orbitals import install

PACKAGE = Path(__file__).resolve().parent
REPO = PACKAGE.parents[2]
CASES = ('ch3s', 'ch3ss')
OLD_NAMES = ('ch3s_ch3ss_reference_checks_20261015',) + tuple(
    'minimum_active_space_molpro_confirmation_20261016' + suffix for suffix in (
        '', '_r2', '_r2_recovery', '_r2_cabs3_gates', '_r2_multistate',
        '_r2_multistate_solver_retry', '_r2_multistate_molpro23_probe', '_shift01'))


def immutable_copy(source, target):
    if target.exists():
        if file_sha256(source) != file_sha256(target):
            raise ValueError(f'Existing prepared artifact differs: {target}; use a fresh destination.')
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def immutable_json(path, value):
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError(f'Existing prepared settings differ: {path}; use a fresh destination.')
    else:
        atomic_json(path, value)


def prepare(root):
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    assets = root/'assets'
    assets.mkdir(exist_ok=True)
    for name in ('reviewed_orbitals.tar.gz', 'orbital_inventory.json',
                 'ch3s_cs_state_model.json', 'ch3ss_cs_state_model.json', 'ch3ss_ss_state_model.json'):
        immutable_copy(PACKAGE/name, assets/name)
    install(assets)
    inputs = []
    for case in CASES:
        directory = root/case
        ensure_current_run(directory, create=True)
        source = REPO/'examples/anl'/f'{case}_anl1_f12'
        surfaces = json.loads((source/f'{case}_surface_model.json').read_text())
        for channel in surfaces['channels']:
            for stage in channel['bottlenecks']:
                p = stage['parameters']
                p['vrc_tst_state_model_file'] = str(assets/Path(p['vrc_tst_state_model_file']).name)
        immutable_json(directory/'surfaces.json', surfaces)
        settings = json.loads((source/f'{case}.json').read_text())
        settings['vrc_tst_surface_model_file'] = str(directory/'surfaces.json')
        if settings.get('vrc_tst_state_model_file'):
            settings['vrc_tst_state_model_file'] = str(assets/Path(settings['vrc_tst_state_model_file']).name)
        settings['input_file'] = str(directory/'input.json')
        immutable_json(directory/'input.json', settings)
        resolved = Parameters(str(directory/'input.json')).par
        immutable_json(directory/'resolved_parameters.json', resolved)
        from kinbot.vrc_multistate import load_surface_model
        load_surface_model(resolved)
        if (resolved['queue_job_limit'] != 10 or resolved['simultaneous_kinbot'] != 1
                or resolved['vrc_tst_max_nodes'] != 10 or resolved['rotdpy_max_jobs'] != 10):
            raise ValueError('The reviewed per-PES admission budget has changed.')
        inputs.append(str(directory/'input.json'))
    files = [p for p in root.rglob('*') if p.is_file() and
             (p.parent == assets or p.parent == assets/'orbitals'
              or p.name in ('input.json', 'surfaces.json', 'resolved_parameters.json'))]
    immutable_json(root/'prepared.json', {'schema': 1, 'inputs': inputs,
        'sha256': {str(p.relative_to(root)): file_sha256(p) for p in sorted(files)}})
    return root


def verify(root):
    record = json.loads((root/'prepared.json').read_text())
    for name, expected in record['sha256'].items():
        if file_sha256(root/name) != expected:
            raise ValueError(f'Prepared artifact changed: {name}')
    return record


def preflight(root):
    from kinbot.anl.automatic import validate_automatic_environment
    from kinbot.anl.site import render_site_setup
    record = verify(root)
    report = {'time_utc': datetime.now(timezone.utc).isoformat(), 'cases': {}}
    # No calculation is submitted by these checks. The actual numerical
    # backends will still be checked inside their allocated jobs.
    for path in record['inputs']:
        par = Parameters(path).par
        report['cases'][Path(path).parent.name] = validate_automatic_environment(par)
    required = {'molpro': {'molpro'}, 'cfour': {'xcfour'}, 'mrcc': {'dmrcc'}}
    setup = root/'site_preflight.sh'
    setup.write_text(render_site_setup(required))
    import shlex
    for backend, programs in required.items():
        script = f'set -e\nexport KINBOT_BACKEND={backend}\nsource {shlex.quote(str(setup))}\n'
        script += '\n'.join('command -v '+program for program in sorted(programs))
        result = subprocess.run(['bash', '-c', script], capture_output=True, text=True, timeout=60)
        if result.returncode:
            raise RuntimeError(f'{backend} preflight failed: {result.stderr or result.stdout}')
        report[backend] = result.stdout.strip()
    report['commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    report['inputs_sha256'] = {path: file_sha256(path) for path in record['inputs']}
    native = root/('backend_preflight.'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    env = dict(os.environ)
    first = report['cases'][CASES[0]]
    env['KINBOT_TEST_MESS'] = first['mess_command']
    env['KINBOT_TEST_MESSPF'] = first['messpf_command']
    result = subprocess.run([sys.executable, '-m', 'pytest', '-q',
        str(REPO/'tests/test_vrc_backend_native.py'), '--basetemp='+str(native)],
        cwd=REPO, env=env, capture_output=True, text=True, timeout=120)
    (root/'backend_preflight.log').write_text(result.stdout+'\n'+result.stderr)
    if result.returncode:
        raise RuntimeError('Native RotdPy/MESS/MESSPF preflight failed; see '+str(root/'backend_preflight.log'))
    report['native_backend_check'] = {'directory': str(native), 'output': result.stdout.strip()}
    atomic_json(root/'preflight.json', report)
    return report


def retire(root, old_parent):
    from kinbot import run_cleanup
    from kinbot.orbital_campaign import stop_controller
    old_parent = Path(old_parent).expanduser().resolve()
    candidates = [old_parent/name for name in OLD_NAMES if (old_parent/name).exists()]
    inventories = {}
    for old in candidates:
        if old == root or old in root.parents or root in old.parents:
            raise ValueError('Production and retired runs must be separate directories.')
        if not ((old/'campaign.json').is_file() or (old/'confirmation.json').is_file()):
            raise ValueError(f'Unrecognized old campaign: {old}')
        inventories[str(old)] = run_cleanup.inspect(old)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    atomic_json(root/f'retirement_inventory.{stamp}.json', inventories)
    reports = []
    for old in candidates:
        if (old/'campaign.json').is_file():
            stop_controller(old)
        report = run_cleanup.stop(old, archive=True, timeout=300.)
        reports.append(report)
        atomic_json(root/f'retirement_results.{stamp}.json', reports)
    return reports


def assert_not_started(root):
    record = verify(root)
    from kinbot.run_cleanup import inspect
    for path in record['inputs']:
        directory = Path(path).parent
        live = inspect(directory)
        if live['jobs'] or live['processes'] or (directory/'launch.json').exists():
            raise RuntimeError(f'{directory} already started; inspect/resume it, do not duplicate the launch.')
    return record


def start(root):
    record = assert_not_started(root)
    launched = []
    for path in record['inputs']:
        directory = Path(path).parent
        env = dict(os.environ)
        env['PATH'] = str(Path(sys.executable).parent)+os.pathsep+env.get('PATH', '')
        env['PYTHONPATH'] = str(REPO)+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
        with (directory/'driver.log').open('ab') as log:
            process = subprocess.Popen([sys.executable, '-u', '-m', 'kinbot.pes', path],
                cwd=directory, env=env, stdin=subprocess.DEVNULL, stdout=log,
                stderr=subprocess.STDOUT, start_new_session=True)
        item = {'pid': process.pid, 'input': path, 'log': str(directory/'driver.log'),
                'time_utc': datetime.now(timezone.utc).isoformat()}
        atomic_json(directory/'launch.json', item)
        launched.append(item)
    return launched


def monitor(root):
    from kinbot.run_cleanup import inspect
    print(datetime.now().astimezone().isoformat(timespec='seconds'))
    for case in CASES:
        directory = root/case
        report = inspect(directory)
        print(f'\n{case}: {len(report["jobs"])} Slurm allocations; {len(report["processes"])} local processes')
        for job in report['jobs']:
            print(job['job_id'], job['state'], job['reason'], job['directory'])
        log = directory/'driver.log'
        if log.exists():
            with log.open('rb') as handle:
                handle.seek(max(0, log.stat().st_size-8192))
                print('\n'.join(handle.read().decode(errors='replace').splitlines()[-12:]))
    print('\nNo active jobs is not a success verdict. Review driver/worker logs and final ANL, VRC and MESS reports.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'preflight', 'launch', 'monitor'])
    parser.add_argument('root', type=Path)
    parser.add_argument('--retire-parent', type=Path,
                        help='Explicitly stop/archive only the listed superseded campaigns under this directory.')
    args = parser.parse_args()
    if args.action == 'monitor':
        monitor(args.root.expanduser().resolve(strict=True))
        sys.exit(0)
    root = prepare(args.root)
    with (root/'.handoff.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = {'prepared': str(root)}
        if args.action in ('preflight', 'launch'):
            result['preflight'] = preflight(root)
        if args.action == 'launch':
            assert_not_started(root)
            if args.retire_parent:
                result['retired'] = retire(root, args.retire_parent)
            result['launched'] = start(root)
        print(json.dumps(result, indent=2))
