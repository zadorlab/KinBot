"""General dispatch behavior, with CH4 as the sole chemistry fixture."""

import json
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write
import numpy as np
import pytest

from examples.anl.ch4_dispatch import ch4_spec
from kinbot.ase_modules.calculators.factory import build_calculator
from kinbot.anl.dispatch import (
    _geometry_hash, _molpro_stack_mw, _molpro_total_mw, _runtime_profile, advance, prepare,
    preflight, retry_failed, run_task, validate_spec,
)


def _write_spec(directory, spec):
    path = Path(directory) / 'spec.json'
    path.write_text(json.dumps(spec))
    return path


def _mock_execution(run_dir, ident, *, geometry=None):
    task_dir = Path(run_dir) / 'tasks' / ident
    task = json.loads((task_dir / 'task.json').read_text())
    spec = task['task']
    if geometry:
        write(task_dir / geometry, Atoms(
            symbols=task['molecule']['symbols'],
            positions=read(task_dir / 'geometry.xyz').positions))
        details = {'final_geometry_sha256': _geometry_hash(read(task_dir / geometry))}
    else:
        details = {}
    names = {'geometry.xyz', 'task.json'}
    if spec['kind'] == 'external':
        names.update({spec['input_name'], spec.get('stdout', 'stdout.txt'),
                      spec.get('stderr', 'stderr.txt')})
        names.update(spec['required_outputs'])
        names.update(spec.get('files_from_env', {}))
    else:
        names.update({'final.xyz', 'optimization.log'})
        if spec['profile']['calculator'] == 'molpro':
            generated = [f'{ident}_step_0001{suffix}'
                         for suffix in ('.inp', '.out', '.log', '.xyz')]
            names.update(generated)
            details['generated_files'] = generated
    for name in names:
        path = task_dir / name
        if not path.exists():
            path.write_text('mock output\n')
    artifacts = {name: hashlib.sha256((task_dir / name).read_bytes()).hexdigest()
                 for name in names}
    (task_dir / 'execution.json').write_text(json.dumps({
        'schema': 1,
        'task_id': ident, 'geometry_sha256': task['geometry_sha256'],
        'status': 'executed', 'details': details, 'artifacts': artifacts,
    }))


def test_ch4_general_graph_stages_geometry_then_all_independent_jobs():
    spec = validate_spec(ch4_spec())
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        assert {path.name for path in (run_dir / 'tasks').iterdir()} == {'l2_geometry'}
        l2_script = (run_dir / 'tasks' / 'l2_geometry' / 'job.slurm').read_text()
        assert l2_script.index('#SBATCH --exclusive') < l2_script.index('set -euo pipefail')

        _mock_execution(run_dir, 'l2_geometry', geometry='final.xyz')
        state = advance(run_dir)
        assert state['tasks']['l2_geometry']['status'] == 'complete'
        assert state['tasks']['l3_geometry']['status'] == 'staged'
        assert 'harmonic' not in state['tasks']
        molpro_record = json.loads(
            (run_dir / 'tasks' / 'l3_geometry' / 'task.json').read_text())
        assert molpro_record['task']['kind'] == 'ase_optimize'
        assert molpro_record['task']['profile']['calculator'] == 'molpro'
        assert molpro_record['task']['profile']['optimizer'] == 'sella'
        assert not list((run_dir / 'tasks' / 'l3_geometry').glob('*.inp'))
        molpro_script = (run_dir / 'tasks' / 'l3_geometry' / 'job.slurm').read_text()
        assert 'export OMP_NUM_THREADS=1' in molpro_script
        assert 'export MKL_NUM_THREADS=1' in molpro_script

        _mock_execution(run_dir, 'l3_geometry', geometry='final.xyz')
        state = advance(run_dir)
        assert all(state['tasks'][task['id']]['status'] == 'staged'
                   for task in spec['tasks'][2:])
        for task in spec['tasks']:
            script = (run_dir / 'tasks' / task['id'] / 'job.slurm').read_text()
            assert script.count('#SBATCH --exclusive') == 1
            assert script.index('#SBATCH --exclusive') < script.index('set -euo pipefail')
        f12 = (run_dir / 'tasks' / 'f12_tz' / 'f12_tz.inp').read_text()
        harmonic = (run_dir / 'tasks' / 'harmonic' / 'harmonic.inp').read_text()
        assert 'frequencies,numerical' in harmonic
        assert 'set,charge=0\nset,spin=0' in harmonic
        assert 'ccsd(t)-f12,scale_trip=1' in f12
        assert 'kb_f12b=energy(2)' in f12
        assert '-m' in spec['tasks'][3]['command']
        assert '{molpro_stack_mw}' in spec['tasks'][3]['command']
        cfour = (run_dir / 'tasks' / 'cfour_dboc' / 'ZMAT').read_text()
        assert 'COORD=CARTESIAN' in cfour
        assert 'DBOC=ON' in cfour
        assert 'MEM_UNIT=MB,MEMORY_SIZE=11200' in cfour
        gaussian = (run_dir / 'tasks' / 'gaussian_vpt2' / 'vpt2.com').read_text()
        assert 'Freq=Anharmonic' in gaussian
        assert 'Opt=(Tight,CalcFC)' in gaussian
        assert '\n0 1\n' in gaussian


def test_ch4_scheduler_respects_node_cap_after_geometry():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        _mock_execution(run_dir, 'l2_geometry', geometry='final.xyz')
        advance(run_dir)
        _mock_execution(run_dir, 'l3_geometry', geometry='final.xyz')
        advance(run_dir)
        submitted = []

        def fake_sbatch(argv, **kwargs):
            ident = Path(kwargs['cwd']).name
            submitted.append(ident)
            return type('Response', (), {'returncode': 0,
                                         'stdout': f'{100 + len(submitted)}\n',
                                         'stderr': ''})()

        with patch('kinbot.anl.dispatch.subprocess.run', side_effect=fake_sbatch):
            state = advance(run_dir, submit=True)
        assert len(submitted) == 3
        assert sum(item['status'] == 'submitted' for item in state['tasks'].values()) == 3
        assert set(submitted) <= {task['id'] for task in ch4_spec()['tasks'][2:]}

        with patch('kinbot.anl.dispatch._job_active', return_value=True), \
                patch('kinbot.anl.dispatch.subprocess.run', side_effect=fake_sbatch):
            advance(run_dir, submit=True)
        assert len(submitted) == 3


def test_generic_external_task_runs_with_isolated_io_and_accepted_geometry():
    spec = {
        'schema': 1, 'name': 'generic-test',
        'molecule': {'symbols': ['H', 'H'],
                     'positions': [[0, 0, 0], [0, 0, 0.74]],
                     'charge': 0, 'multiplicity': 1},
        'limits': {'max_nodes': 1, 'max_cores_per_node': 2,
                   'max_memory_mb_per_node': 1024},
        'tasks': [{
            'id': 'geometry', 'kind': 'external', 'backend': 'fake',
            'geometry_from': 'initial', 'geometry_output': 'final.xyz',
            'resources': {'cores': 1, 'memory_mb': 512, 'walltime': '00:10:00'},
            'input_name': 'job.inp', 'input_template': '{{XYZ}}\n',
            'command': [sys.executable, '-c',
                        "from pathlib import Path; "
                        "Path('final.xyz').write_text(Path('geometry.xyz').read_text()); "
                        "Path('job.out').write_text('DONE')"],
            'stdout': 'launcher.stdout', 'stderr': 'launcher.stderr',
            'required_outputs': ['job.out', 'final.xyz'],
            'success_marker': {'file': 'job.out', 'contains': 'DONE'},
        }],
    }
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        result = run_task(run_dir / 'tasks' / 'geometry' / 'task.json')
        assert result['status'] == 'executed'
        assert 'job.out' in result['artifacts']
        state = advance(run_dir)
        assert state['tasks']['geometry']['status'] == 'complete'
        with pytest.raises(RuntimeError, match='already has execution'):
            run_task(run_dir / 'tasks' / 'geometry' / 'task.json')


def test_program_failure_text_is_rejected_even_with_zero_exit_status():
    spec = ch4_spec()
    spec['tasks'] = [{
        'id': 'check', 'kind': 'external', 'backend': 'fake',
        'resources': {'cores': 1, 'memory_mb': 512,
                      'walltime': '00:10:00'},
        'input_name': 'check.inp', 'input_template': '{{XYZ}}\n',
        'command': [sys.executable, '-c',
                    "from pathlib import Path; "
                    "Path('check.out').write_text('DONE\\nFatal error\\n')"],
        'required_outputs': ['check.out'],
        'success_marker': {'file': 'check.out', 'contains': 'DONE'},
        'failure_markers': [{'file': 'check.out', 'contains': 'Fatal error'}],
    }]
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        with pytest.raises(RuntimeError, match='Failure marker'):
            run_task(run_dir / 'tasks' / 'check' / 'task.json')
        assert advance(run_dir)['tasks']['check']['status'] == 'failed'


def test_ch4_ase_sella_runner_with_analytic_toy_calculator():
    class Quadratic(Calculator):
        implemented_properties = ['energy', 'forces']

        def __init__(self, target):
            super().__init__()
            self.target = target

        def calculate(self, atoms=None, properties=('energy',),
                      system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            assert atoms.info['charge'] == 0
            assert atoms.info['spin'] == 1
            delta = atoms.positions - self.target
            self.results = {'energy': float(np.sum(delta * delta)),
                            'forces': -2 * delta}

    spec = ch4_spec()
    spec['tasks'] = [spec['tasks'][0]]
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        target = read(run_dir / 'tasks' / 'l2_geometry' / 'geometry.xyz').positions.copy()
        target[1, 0] += 0.04
        with patch('kinbot.ase_modules.calculators.factory.build_calculator',
                   return_value=Quadratic(target)):
            result = run_task(run_dir / 'tasks' / 'l2_geometry' / 'task.json')
        assert result['status'] == 'executed'
        assert result['details']['optimizer'] == 'sella'
        assert (run_dir / 'tasks' / 'l2_geometry' / 'final.xyz').is_file()
        assert advance(run_dir)['tasks']['l2_geometry']['status'] == 'complete'


def test_ch4_gaussian_ase_input_rendering_without_gaussian_executable():
    spec = ch4_spec()
    task = spec['tasks'][0]
    atoms = Atoms(symbols=spec['molecule']['symbols'],
                  positions=spec['molecule']['positions'])
    profile = _runtime_profile(task)
    with TemporaryDirectory() as directory:
        calc = build_calculator(profile, directory,
                                {'name': 'ch4', 'charge': 0, 'mult': 1})
        calc.write_input(atoms, ['energy', 'forces'])
        gaussian = (Path(directory) / 'ch4.com').read_text()
    assert '%nprocshared=4' in gaussian
    assert '%mem=11200MB' in gaussian
    assert 'B2PLYP/cc-pVTZ' in gaussian
    assert 'EmpiricalDispersion(GD3BJ)' in gaussian
    assert '\n0 1\n' in gaussian
    assert 'force' in gaussian.lower()


def test_validation_rejects_cycles_and_oversized_exclusive_jobs():
    spec = ch4_spec()
    spec['tasks'][0]['geometry_from'] = 'l3_geometry'
    with pytest.raises(ValueError, match='cycle'):
        validate_spec(spec)
    spec = ch4_spec()
    spec['tasks'][1]['resources']['cores'] = 17
    with pytest.raises(ValueError, match='exceeds'):
        validate_spec(spec)
    spec = ch4_spec()
    spec['tasks'][1]['profile']['calculator_kwargs'] = {'nproc': 4}
    with pytest.raises(ValueError, match='Molpro nproc'):
        validate_spec(spec)
    spec = ch4_spec()
    spec['tasks'][6]['input_name'] = 'dboc.inp'
    with pytest.raises(ValueError, match='ZMAT'):
        validate_spec(spec)
    spec = ch4_spec()
    spec['tasks'][3]['stdout'] = 'execution.json'
    with pytest.raises(ValueError, match='reserved'):
        validate_spec(spec)
    spec = ch4_spec()
    spec['tasks'][0]['profile']['calculator_kwargs']['nprocshared'] = 8
    with pytest.raises(ValueError, match='must match Slurm cores'):
        validate_spec(spec)
    assert _molpro_total_mw({'cores': 8, 'memory_mb': 32000}) == 1800
    assert _molpro_stack_mw({'cores': 8, 'memory_mb': 32000}) == 225
    spec = ch4_spec()
    spec['tasks'][1]['resources']['memory_mb'] = 2000
    with pytest.raises(ValueError, match='per-process overhead'):
        validate_spec(spec)


def test_missing_sbatch_does_not_leave_uncertain_submission_state():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        with patch('kinbot.anl.dispatch.subprocess.run',
                   side_effect=FileNotFoundError('sbatch')):
            with pytest.raises(FileNotFoundError):
                advance(run_dir, submit=True)
        state = json.loads((run_dir / 'state.json').read_text())
        assert state['tasks']['l2_geometry']['status'] == 'staged'


def test_completed_geometry_or_artifact_mutation_blocks_children():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        _mock_execution(run_dir, 'l2_geometry', geometry='final.xyz')
        advance(run_dir)
        final = run_dir / 'tasks' / 'l2_geometry' / 'final.xyz'
        final.write_text(final.read_text().replace('0.630', '0.635'))
        with pytest.raises(RuntimeError, match='artifact|geometry'):
            advance(run_dir)


def test_staged_input_mutation_blocks_submission():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        job = run_dir / 'tasks' / 'l2_geometry' / 'job.slurm'
        job.write_text(job.read_text().replace('--exclusive', '--oversubscribe'))
        with pytest.raises(RuntimeError, match='staged job.slurm changed'):
            advance(run_dir, submit=True)


def test_failed_task_can_be_archived_and_retried():
    spec = ch4_spec()
    spec['tasks'] = [{
        'id': 'sample', 'kind': 'external', 'backend': 'fake',
        'resources': {'cores': 1, 'memory_mb': 512,
                      'walltime': '00:10:00'},
        'input_name': 'sample.inp',
        'input_template': '{{XYZ}}\n',
        'command': [sys.executable, '-c', 'raise SystemExit(2)'],
        'required_outputs': ['sample.out'],
    }]
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        with pytest.raises(RuntimeError, match='Program exited'):
            run_task(run_dir / 'tasks' / 'sample' / 'task.json')
        assert advance(run_dir)['tasks']['sample']['status'] == 'failed'
        archive = retry_failed(run_dir, 'sample')
        assert (archive / 'execution.json').is_file()
        assert not (run_dir / 'tasks' / 'sample' / 'execution.json').exists()
        assert advance(run_dir)['tasks']['sample']['status'] == 'staged'


def test_preflight_reports_missing_program_before_submission():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        spec = ch4_spec()
        spec['tasks'][0]['profile']['command'] = str(root / 'missing_g16')
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        with patch('kinbot.anl.dispatch.shutil.which', return_value='/bin/true'):
            with pytest.raises(RuntimeError, match='Missing executable'):
                preflight(run_dir)


def test_preflight_sources_site_setup_and_validates_slurm_without_submitting():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        bindir = root / 'bin'
        bindir.mkdir()
        for program in ('sbatch', 'squeue'):
            executable = bindir / program
            executable.write_text('#!/usr/bin/env bash\nexit 0\n')
            executable.chmod(0o755)
        backend_programs = {'gaussian': 'g16', 'molpro': 'molpro',
                            'cfour': 'xcfour'}
        for backend, program in backend_programs.items():
            backend_dir = root / backend
            backend_dir.mkdir()
            executable = backend_dir / program
            executable.write_text('#!/usr/bin/env bash\nexit 0\n')
            executable.chmod(0o755)
        genbas = root / 'GENBAS'
        genbas.write_text('basis fixture\n')
        (run_dir / 'site_setup.sh').write_text(
            'case "$KINBOT_BACKEND" in\n'
            + ''.join(f'  {backend}) export PATH={root / backend}:$PATH ;;\n'
                      for backend in backend_programs)
            + '  *) exit 7 ;;\n'
            + 'esac\n'
            + 'if [ "$KINBOT_BACKEND" = cfour ]; then '
            + f'export CFOUR_GENBAS={genbas}; fi\n'
            f'export PYTHONPATH={Path(__file__).resolve().parents[1]}:${{PYTHONPATH:-}}\n'
        )
        with patch.dict(os.environ, {'PATH': f'{bindir}:{os.environ["PATH"]}'}):
            result = preflight(run_dir)
        assert result['slurm_scripts_tested'] == 1
        assert result['exclusive_jobs_checked'] == 1
        assert result['programs'] == ['g16', 'molpro', 'xcfour']


def test_prepare_discovers_vendor_setup_and_cfour_genbas():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        gaussian = root / 'gaussian' / 'g16'
        molpro = root / 'molpro' / 'bin'
        cfour = root / 'cfour' / 'bin'
        for directory, program in ((gaussian, 'g16'), (molpro, 'molpro'),
                                   (cfour, 'xcfour')):
            directory.mkdir(parents=True)
            executable = directory / program
            executable.write_text('#!/usr/bin/env bash\nexit 0\n')
            executable.chmod(0o755)
        (gaussian / 'bsd').mkdir()
        (gaussian / 'bsd' / 'g16.profile').write_text('export PROFILE_SOURCED=yes\n')
        (root / 'cfour' / 'basis').mkdir()
        genbas = root / 'cfour' / 'basis' / 'GENBAS'
        genbas.write_text('basis fixture\n')
        (root / 'cfour' / 'payload').mkdir()
        (cfour / 'xcfour').replace(root / 'cfour' / 'payload' / 'xcfour')
        (cfour / 'xcfour').symlink_to(root / 'cfour' / 'payload' / 'xcfour')
        spec = ch4_spec()
        for task in spec['tasks']:
            task['resources']['partition'] = 'chosen_by_user'
        base_path = os.environ['PATH']
        environment = {'PATH': f'{gaussian}:{molpro}:{cfour}:{base_path}',
                       'LOADEDMODULES': 'molpro/molpro24:cfour/2.1',
                       'CFOUR_GENBAS': ''}
        with patch.dict(os.environ, environment):
            run_dir = prepare(_write_spec(root, spec), root / 'run')
            setup = (run_dir / 'site_setup.sh').read_text()
            assert 'module load molpro/molpro24' in setup
            assert 'module load cfour/2.1' in setup
            for backend in ('gaussian', 'molpro', 'cfour'):
                program = {'gaussian': 'g16', 'molpro': 'molpro',
                           'cfour': 'xcfour'}[backend]
                command = (
                    'set -euo pipefail\n'
                    'module() { test "$1" = load; export LOADED_MODULE="$2"; }\n'
                    f'export KINBOT_BACKEND={backend}\n'
                    'source site_setup.sh\n'
                    f'command -v {program}\n'
                    'printf "|%s|%s|%s|%s" "${PROFILE_SOURCED:-}" '
                    '"${g16root:-}" "${CFOUR_GENBAS:-}" "${LOADED_MODULE:-}"\n'
                )
                child_env = os.environ.copy()
                child_env['PATH'] = base_path
                result = subprocess.run(['bash', '-c', command], cwd=run_dir,
                                        env=child_env, capture_output=True,
                                        text=True, check=True)
                if backend == 'gaussian':
                    assert result.stdout.startswith(f'{gaussian / "g16"}\n')
                    assert f'|yes|{root / "gaussian"}|' in result.stdout
                elif backend == 'cfour':
                    assert result.stdout.startswith(f'{cfour / "xcfour"}\n')
                    assert f'|{genbas}|cfour/2.1' in result.stdout
                    cfour_command = command
                else:
                    assert result.stdout.startswith(f'{molpro / "molpro"}\n')
                    assert result.stdout.endswith('|molpro/molpro24')
            override = root / 'different-GENBAS'
            override.write_text('override fixture\n')
            child_env['CFOUR_GENBAS'] = str(override)
            result = subprocess.run(['bash', '-c', cfour_command],
                                    cwd=run_dir, env=child_env, capture_output=True,
                                    text=True, check=True)
            assert f'|{override}|cfour/2.1' in result.stdout


def test_prepare_selects_fitting_slurm_partition_and_keeps_explicit_choice():
    display = ('short-cpu*|96|126000|30:00|up\n'
               'day-long-cpu|96|126000|1-00:00:00|up\n'
               'week-long-cpu|96|126000|7-00:00:00|up\n'
               'drained|96|126000|31-00:00:00|down\n')
    spec = ch4_spec()
    spec['tasks'][0]['resources']['partition'] = 'week-long-cpu'
    actual_which = shutil.which

    def which(program):
        return '/usr/bin/sinfo' if program == 'sinfo' else actual_which(program)

    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        with patch('kinbot.anl.site.shutil.which', side_effect=which), \
                patch('kinbot.anl.site.subprocess.run', return_value=
                      subprocess.CompletedProcess([], 0, display, '')):
            run_dir = prepare(_write_spec(root, spec), root / 'run')
        prepared = json.loads((run_dir / 'workflow.json').read_text())
        partitions = {task['id']: task['resources']['partition']
                      for task in prepared['tasks']}
        assert partitions['l2_geometry'] == 'week-long-cpu'
        assert set(partitions.values()) == {'week-long-cpu', 'day-long-cpu'}
        script = (run_dir / 'tasks' / 'l2_geometry' / 'job.slurm').read_text()
        assert '#SBATCH --partition=week-long-cpu\n' in script
        with patch('kinbot.anl.site.shutil.which', side_effect=which), \
                patch('kinbot.anl.site.subprocess.run', return_value=
                      subprocess.CompletedProcess([], 0, display.splitlines()[0] + '\n', '')):
            with pytest.raises(RuntimeError, match='no available Slurm partition'):
                prepare(_write_spec(root, spec), root / 'impossible')
        assert not (root / 'impossible').exists()
