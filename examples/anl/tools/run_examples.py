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
    inputs = []
    files = set()
    for case in CASES:
        directory = root/case
        ensure_current_run(directory, create=True)
        source = REPO/'examples/anl'/f'{case}_anl1_f12'
        settings = json.loads((source/f'{case}.json').read_text())
        surfaces = settings['vrc_tst_surface_model']
        for channel in surfaces['channels']:
            for stage in channel['bottlenecks']:
                p = stage['parameters']
                for model in p['vrc_tst_state_model']['basis_models'].values():
                    selection = model['orbital_selection']
                    seed = (source/selection['path']).resolve()
                    if not seed.is_relative_to(source.resolve()):
                        raise ValueError('Example orbital data must be inside its own directory.')
                    if file_sha256(seed) != selection['sha256']:
                        raise ValueError(f'Checkpoint differs from the single input: {seed}')
                    relative = Path('orbital_reviews/orbitals')/seed.name
                    immutable_copy(seed, directory/relative)
                    files.add(directory/relative)
                    selection['path'] = relative.as_posix()
        settings['input_file'] = str(directory/'input.json')
        immutable_json(directory/'input.json', settings)
        resolved = Parameters(str(directory/'input.json')).par
        immutable_json(directory/'resolved_parameters.json', resolved)
        files.update((directory/'input.json', directory/'resolved_parameters.json'))
        from kinbot.vrc_multistate import load_surface_model
        load_surface_model(resolved)
        if (resolved['queue_job_limit'] != 10 or resolved['simultaneous_kinbot'] != 1
                or resolved['vrc_tst_max_nodes'] != 10 or resolved['rotdpy_max_jobs'] != 10):
            raise ValueError('The reviewed per-PES admission budget has changed.')
        inputs.append(str(directory/'input.json'))
    # Only declared launch inputs belong in this immutable manifest. On a
    # repeated handoff, runtime reference inputs/WFUs must not change it.
    immutable_json(root/'prepared.json', {'schema': 1, 'inputs': inputs,
        'sha256': {str(p.relative_to(root)): file_sha256(p) for p in sorted(files)}})
    return root


def verify(root):
    record = json.loads((root/'prepared.json').read_text())
    for name, expected in record['sha256'].items():
        if file_sha256(root/name) != expected:
            raise ValueError(f'Prepared artifact changed: {name}')
    return record


def archive_legacy_example_layout(root):
    """Preserve obsolete untracked files that git pull cannot relocate.

    These are example/review artifacts, not jobs or production configuration.
    Check every directory against the current index before moving any of them.
    """
    names = ('caspt2_production', 'orbital_reviews', 'acetaldehyde_cbh_anl0_f12')
    candidates = [REPO/'examples/anl'/name for name in names
                  if (REPO/'examples/anl'/name).exists()]
    if not candidates:
        return []
    for path in candidates:
        if path.is_symlink() or not path.is_dir():
            raise ValueError(f'Unexpected obsolete example path: {path}')
        tracked = subprocess.check_output(
            ['git', 'ls-files', '-z', '--', str(path.relative_to(REPO))], cwd=REPO)
        if tracked:
            raise ValueError(f'Cannot archive tracked example files: {path}; update the repository first.')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive = root.parent/('anl_example_layout.stopped.'+stamp)
    archive.mkdir()
    report = []
    for path in candidates:
        destination = archive/path.name
        path.rename(destination)
        report.append({'source': str(path), 'archive': str(destination)})
        atomic_json(root/'repository_layout_archive.json', report)
    return report


def preflight(root):
    from kinbot.anl.automatic import validate_automatic_environment
    from kinbot.anl.site import render_site_setup
    record = verify(root)
    report = {'time_utc': datetime.now(timezone.utc).isoformat(), 'cases': {}}
    report['rendered_gaussian_inputs'] = render_gaussian_inputs(root, record)
    report['rendered_vrc_inputs'] = render_vrc_inputs(root, record)
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


def render_vrc_inputs(root, record):
    """Resolve every declared channel/root/basis through production renderers."""
    from ase.data import atomic_numbers
    from copy import deepcopy
    from kinbot.vrc_multistate import load_surface_model
    from kinbot.vrc_correction import correction_components, spin_splitting_settings
    from kinbot.molpro import stage_vrc_state_model, _vrc_multireference_method, vrc_multireference_energy_contract
    paths = []
    for filename in record['inputs']:
        par = Parameters(filename).par
        surface, _ = load_surface_model(par)
        for ci, channel in enumerate(surface['channels']):
            for stage in channel['bottlenecks']:
                envelope = stage['parameters']['vrc_tst_state_model']
                molecule = envelope['molecule']
                symbols, charge = molecule['symbols'], molecule['charge']
                mults = [e['multiplicity'] for e in next(iter(envelope['basis_models'].values()))['state_average']]
                # Render the declared target/reference blocks. Actual fragment
                # coupling is checked against accepted products by the driver.
                for state in stage['states']:
                    leaf = deepcopy(par)
                    leaf.update(deepcopy(stage['parameters']))
                    leaf['vrc_tst_selected_roots'] = {str(state['multiplicity']): state['root']}
                    leaf['vrc_tst_spin_splitting'].update(target_multiplicity=state['multiplicity'],
                                                       reference_multiplicity=max(mults))
                    settings = spin_splitting_settings(leaf)
                    if settings['reference_multiplicity'] not in mults:
                        raise ValueError('Reference multiplicity missing from the reviewed model.')
                    for i, component in enumerate(correction_components(leaf, symbols=symbols, charge=charge, settings=settings)):
                        if vrc_multireference_energy_contract(component['method']) is None:
                            continue
                        directory = root/'input_preflight'/Path(filename).parent.name/f'channel{ci}_{stage["name"]}_root{state["root"]}'
                        model = stage_vrc_state_model(component['state_model'], directory, f'vrc{i}')
                        deck, _ = _vrc_multireference_method(component['method'],
                            sum(atomic_numbers[s] for s in symbols)-charge, component['multiplicity']-1,
                            state_model=model, shift=component.get('shift'), ipea=component.get('ipea'))
                        path = directory/f'vrc{i}.inp'
                        path.write_text(deck)
                        paths.append({'path': str(path.relative_to(root)), 'sha256': file_sha256(path),
                                      'method': component['method'], 'basis': component['basis'],
                                      'active_indices': model['orbital_selection']['active_indices']})
    return {'execution': 'render only; Molpro not executed', 'files': paths}


def render_gaussian_inputs(root, record):
    """Exercise the actual profile merge/writer, without running Gaussian."""
    from copy import deepcopy
    from ase import Atoms
    from ase.data import atomic_numbers
    from kinbot.qc import QuantumChemistry
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.ase_modules.io.gaussian import write_gaussian_in
    from kinbot.ase_modules.calculators.gaussian import fixed_frequency_parameters
    variants = {'well': {}, 'ts': {'ts': 1, 'step': 1, 'max_step': 1},
                'irc': {'irc': 'Forward'}, 'rotor': {'hir': 1},
                'rigid_rotor': {'hir': 1, 'rigid': 1}, 'scan': {'scan': 1}}
    paths = []
    previous = Path.cwd()
    previous_exclusions = os.environ.get('SBATCH_EXCLUDE')
    try:
        for filename in record['inputs']:
            directory = root/'input_preflight'/Path(filename).parent.name
            directory.mkdir(parents=True, exist_ok=True)
            os.chdir(directory)
            par = Parameters(filename).par
            species = StationaryPoint('preflight', par['charge'], par['mult'],
                                      smiles=par['smiles'], structure=par['structure'])
            atoms = Atoms(species.atom, positions=species.geom)
            electrons = sum(atomic_numbers[a] for a in species.atom)-par['charge']
            router = QuantumChemistry(par)
            for level in ('l1', 'l2'):
                qc = getattr(router, level)
                if qc.qc != 'gauss':
                    raise ValueError('Reviewed production L1 and L2 must use Gaussian.')
                for name, flags in variants.items():
                    options = qc.merge_kwargs(qc.get_qc_arguments(
                        'preflight', par['mult'], par['charge'], electrons,
                        high_level=int(level == 'l2'), **flags))
                    options.pop('label', None)
                    if qc.use_sella:
                        # Native force/Hessian calls use the external optimizer;
                        # do not leave Gaussian optimization/IRC directives in previews.
                        for key in ('opt', 'freq'):
                            options.pop(key, None)
                        if flags.get('irc'):
                            for key in ('irc', 'geom', 'guess'):
                                options.pop(key, None)
                    for mode, props, kwargs in (
                            ('energy', ['energy'], options),
                            ('forces', ['energy', 'forces'], options),
                            ('hessian', ['energy'], fixed_frequency_parameters(options))):
                        path = directory/f'{level}_{name}_{mode}.com'
                        with path.open('w') as handle:
                            write_gaussian_in(handle, atoms, properties=props, **deepcopy(kwargs))
                        paths.append({'path': str(path.relative_to(root)), 'sha256': file_sha256(path)})
    finally:
        os.chdir(previous)
        # QC construction propagates exclusions to child submitters. Preview
        # construction must not leak one case's settings into later preflights.
        if previous_exclusions is None:
            os.environ.pop('SBATCH_EXCLUDE', None)
        else:
            os.environ['SBATCH_EXCLUDE'] = previous_exclusions
    return {'execution': 'render only; Gaussian not executed', 'files': paths}


def retire_production(root, previous, *, known=(), allow_partial=False):
    """Stop each marked PES, then archive the complete old handoff tree."""
    import socket
    from kinbot import run_cleanup
    root = Path(root).resolve(strict=True)
    previous = Path(previous).expanduser().resolve(strict=True)
    if previous == root or previous in root.parents or root in previous.parents:
        raise ValueError('Replacement and previous production must be separate directories.')
    expected = {previous/case/'input.json' for case in CASES}
    manifest = previous/'prepared.json'
    if allow_partial and not manifest.exists():
        # An interrupted prepare may have created only the first case.
        actual = {p for p in expected if p.parent.exists()}
        if not actual:
            raise ValueError(f'No marked PES directories in {previous}; not archived.')
        for path in actual:
            ensure_current_run(path.parent)
    else:
        record = json.loads(manifest.read_text())
        actual = {Path(p) for p in record.get('inputs', [])}
        if actual != expected:
            raise ValueError('Unrecognized previous production input directories.')
    if any(path.is_symlink() or path.parent.is_symlink() for path in actual):
        raise ValueError('Unrecognized previous production input directories.')
    for path in actual:
        for name in ('pes_status.json', 'launch.json'):
            marker = path.parent/name
            if marker.exists():
                host = json.loads(marker.read_text()).get('host')
                if host and host.split('.')[0] != socket.gethostname().split('.')[0]:
                    raise RuntimeError(f'{path.parent}: retire from controller host {host}; '
                                       'local /proc cannot exclude a live remote driver.')
    # Check every tree before the first stop. An unavailable scheduler aborts.
    inventories = {str(p.parent): run_cleanup.inspect(p.parent) for p in sorted(actual)}
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    audit = root/f'production_retirement.{stamp}.json'
    report = {'previous': str(previous), 'inventories': inventories, 'stopped': []}
    atomic_json(audit, report)
    for index, path in enumerate(sorted(actual)):
        options = {'known': known} if index == 0 and known else {}
        report['stopped'].append(run_cleanup.stop(path.parent, archive=False, timeout=300., **options))
        atomic_json(audit, report)
    # Recheck both before moving their shared parent.
    for path in sorted(actual):
        live = run_cleanup.inspect(path.parent)
        if live['jobs'] or live['processes']:
            raise RuntimeError(f'Work appeared during retirement: {path.parent}; not archived.')
    destination = previous.with_name(previous.name+'.stopped.'+stamp)
    report['archive_intended'] = str(destination)
    atomic_json(audit, report)
    previous.rename(destination)
    report['archive'] = str(destination)
    atomic_json(audit, report)
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
    return start_inputs(record['inputs'])


def start_inputs(inputs):
    """Caller must hold the handoff lock and check these inputs are unstarted."""
    import socket
    launched = []
    for path in inputs:
        directory = Path(path).parent
        env = dict(os.environ)
        env['PATH'] = str(Path(sys.executable).parent)+os.pathsep+env.get('PATH', '')
        env['PYTHONPATH'] = str(REPO)+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
        with (directory/'driver.log').open('ab') as log:
            process = subprocess.Popen([sys.executable, '-u', '-m', 'kinbot.pes', path],
                cwd=directory, env=env, stdin=subprocess.DEVNULL, stdout=log,
                stderr=subprocess.STDOUT, start_new_session=True)
        item = {'pid': process.pid, 'input': path, 'log': str(directory/'driver.log'),
                'host': socket.gethostname(),
                'time_utc': datetime.now(timezone.utc).isoformat()}
        atomic_json(directory/'launch.json', item)
        launched.append(item)
    return launched


def monitor(root, args=None):
    from kinbot import run_monitor
    if args is None:
        print(run_monitor.render(run_monitor.snapshot(root)))
    else:
        run_monitor.run(root, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'preflight', 'launch', 'resume', 'monitor'])
    parser.add_argument('root', type=Path)
    parser.add_argument('--retire-parent', type=Path,
                        help='Explicitly stop/archive only the listed superseded campaigns under this directory.')
    parser.add_argument('--replace-production', type=Path,
                        help='On launch, stop both PES trees in this previous handoff and archive its parent.')
    from kinbot.run_monitor import add_arguments
    add_arguments(parser)
    args = parser.parse_args()
    if (args.retire_parent or args.replace_production) and args.action != 'launch':
        parser.error('Retirement options require launch.')
    if args.retire_parent and args.replace_production:
        parser.error('Choose old campaigns or a previous production handoff, not both.')
    if args.action == 'monitor':
        monitor(args.root.expanduser().resolve(strict=True), args)
        sys.exit(0)
    root = (args.root.expanduser().resolve(strict=True)
            if args.action == 'resume' else prepare(args.root))
    with (root/'.handoff.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = {'prepared': str(root)}
        if args.action in ('preflight', 'launch', 'resume'):
            result['preflight'] = preflight(root)
        if args.action == 'resume':
            from resume_downstream import repair_before_resume
            result['repair'] = repair_before_resume(root, verify(root))
            result['launched'] = start(root)
        if args.action == 'launch':
            assert_not_started(root)
            result['archived_repository_layout'] = archive_legacy_example_layout(root)
            if args.retire_parent:
                result['retired'] = retire(root, args.retire_parent)
            if args.replace_production:
                result['retired_production'] = retire_production(root, args.replace_production)
            result['launched'] = start(root)
        # Keep the terminal handoff short; all paths and check results remain
        # in the machine-readable report for later review.
        report_path = root / (args.action + '_report.json')
        atomic_json(report_path, result)
        summary = {'run': str(root), 'report': str(report_path)}
        if 'preflight' in result:
            summary['preflight'] = 'passed; details in preflight.json'
        if 'repair' in result:
            repair = result['repair']
            summary['archive'] = repair['archive']
            summary['accepted_composite_tasks_preserved'] = sum(
                g['complete_tasks_preserved'] for g in repair['graphs'])
            summary['invalid_vrc_trees_archived'] = len(repair['vrc_archived'])
        if 'launched' in result:
            summary['launched'] = result['launched']
        print(json.dumps(summary, indent=2))
