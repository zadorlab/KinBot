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
    _geometry_hash, _molpro_stack_mw, _molpro_total_mw, _runtime_profile, advance, main, prepare,
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
        assert '#SBATCH --ntasks=1\n' in l2_script
        assert '#SBATCH --cpus-per-task=4\n' in l2_script

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
        assert '#SBATCH --ntasks=8\n' in molpro_script
        assert '#SBATCH --cpus-per-task=1\n' in molpro_script
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
            backend = (task['backend'] if task['kind'] == 'external'
                       else task['profile']['calculator'])
            if backend == 'molpro':
                assert f"#SBATCH --ntasks={task['resources']['cores']}\n" in script
                assert '#SBATCH --cpus-per-task=1\n' in script
            else:
                assert '#SBATCH --ntasks=1\n' in script
                assert f"#SBATCH --cpus-per-task={task['resources']['cores']}\n" in script
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
        assert 'MEM_UNIT=MB\nMEMORY_SIZE=11200)' in cfour
        assert max(map(len, cfour.splitlines())) <= 72
        gaussian = (run_dir / 'tasks' / 'gaussian_vpt2' / 'vpt2.com').read_text()
        vpt2_task = next(task for task in spec['tasks']
                         if task['id'] == 'gaussian_vpt2')
        assert vpt2_task['geometry_from'] == 'l2_geometry'
        assert vpt2_task['depends_on'] == ['l3_geometry']
        assert state['tasks']['gaussian_vpt2']['geometry_from'] == 'l2_geometry'
        assert 'B2PLYP/cc-pVTZ' in gaussian
        assert 'EmpiricalDispersion=GD3BJ' in gaussian
        assert 'Freq=Anharmonic' in gaussian
        assert 'Opt=' not in gaussian
        assert '\n0 1\n' in gaussian


def test_frequency_only_vpt2_requires_matching_l2_geometry_surface():
    spec = ch4_spec()
    vpt2 = next(task for task in spec['tasks'] if task['id'] == 'gaussian_vpt2')
    vpt2['result_parser']['method'] = 'B3LYP'
    vpt2['input_template'] = vpt2['input_template'].replace('B2PLYP/', 'B3LYP/')
    with pytest.raises(ValueError, match='matching Gaussian geometry source'):
        validate_spec(spec)

    lower = ch4_spec()
    l2 = next(task for task in lower['tasks'] if task['id'] == 'l2_geometry')
    l2['profile']['method'] = 'B3LYP'
    l2['profile']['calculator_kwargs'].pop('EmpiricalDispersion')
    lower_vpt2 = next(task for task in lower['tasks']
                      if task['id'] == 'gaussian_vpt2')
    lower_vpt2['input_template'] = lower_vpt2['input_template'].replace(
        'B2PLYP/', 'B3LYP/').replace('EmpiricalDispersion=GD3BJ ', '')
    lower_vpt2['result_parser']['method'] = 'B3LYP'
    lower_vpt2['result_parser'].pop('dispersion')
    validate_spec(lower)


def test_retry_rewrites_old_cfour_keyword_line_without_changing_workflow():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        original = ch4_spec()
        cfour = next(task for task in original['tasks'] if task['id'] == 'cfour_dboc')
        cfour['geometry_from'] = 'initial'
        cfour['input_template'] = (
            'CH4\n{{CARTESIAN}}\n\n'
            '*CFOUR(CALC=SCF,BASIS=cc-pVTZ,DBOC=ON,COORD=CARTESIAN,'
            'UNITS=ANGSTROM,CHARGE={{CHARGE}},MULTIPLICITY={{MULT}},'
            'MEM_UNIT=MB,MEMORY_SIZE={{WORK_MEMORY_MB}})\n')
        spec = {'schema': 1, 'name': 'cfour-retry',
                'molecule': original['molecule'], 'limits': {'max_nodes': 1},
                'tasks': [cfour]}
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        workflow_hash = hashlib.sha256((run_dir / 'workflow.json').read_bytes()).hexdigest()
        task_dir = run_dir / 'tasks' / 'cfour_dboc'
        head, keyword_block = (task_dir / 'ZMAT').read_text().split('*CFOUR(', 1)
        old_input = (head + '*CFOUR('
                     + ','.join(keyword_block.rstrip().removesuffix(')').splitlines())
                     + ')\n')
        (task_dir / 'ZMAT').write_text(old_input)
        state_path = run_dir / 'state.json'
        state = json.loads(state_path.read_text())
        state['tasks']['cfour_dboc'].update(
            status='failed', input_sha256=hashlib.sha256(old_input.encode()).hexdigest())
        state_path.write_text(json.dumps(state))
        archive = retry_failed(run_dir, 'cfour_dboc')
        assert (archive / 'ZMAT').read_text() == old_input
        new_input = (task_dir / 'ZMAT').read_text()
        assert 'CHARGE=0\nMULTIPLICITY=1\nMEM_UNIT=MB' in new_input
        assert max(map(len, new_input.splitlines())) <= 72
        assert hashlib.sha256((run_dir / 'workflow.json').read_bytes()).hexdigest() == workflow_hash


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


def test_selected_submission_keeps_other_ready_jobs_staged():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        _mock_execution(run_dir, 'l2_geometry', geometry='final.xyz')
        advance(run_dir)
        _mock_execution(run_dir, 'l3_geometry', geometry='final.xyz')
        advance(run_dir)
        with patch('kinbot.anl.dispatch.subprocess.run', return_value=type(
                'Response', (), {'returncode': 0, 'stdout': '123\n', 'stderr': ''})()):
            state = advance(run_dir, submit=True, submit_only={'molpro_dz_sp'})
        assert state['tasks']['molpro_dz_sp']['status'] == 'submitted'
        assert all(state['tasks'][task['id']]['status'] == 'staged'
                   for task in ch4_spec()['tasks'][2:] if task['id'] != 'molpro_dz_sp')
        with pytest.raises(ValueError, match='Unknown submission'):
            advance(run_dir, submit_only={'missing_task'})


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
    assert 'integral(UltraFine)' in gaussian
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


def test_auto_cores_use_safe_node_memory_and_efficient_rank_counts():
    spec = ch4_spec()
    spec['limits'] = {'max_nodes': 2}
    spec['tasks'] = [dict(spec['tasks'][1], geometry_from='initial'),
                     spec['tasks'][0]]
    spec['tasks'][0]['resources'] = {'walltime': '24:00:00'}
    spec['tasks'][1]['resources'] = {'walltime': '08:00:00'}
    node_groups = [
        {'name': 'day', 'default': True, 'cores': 96,
         'memory_mb': 126000, 'seconds': 86400},
        {'name': 'day', 'default': True, 'cores': 64,
         'memory_mb': 64000, 'seconds': 86400},
    ]
    with TemporaryDirectory() as temporary, \
            patch('kinbot.anl.site._partitions', return_value=node_groups):
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        resolved = json.loads((run_dir / 'workflow.json').read_text())
        molpro = resolved['tasks'][0]['resources']
        gaussian = resolved['tasks'][1]['resources']
        assert molpro['partition'] == gaussian['partition'] == 'day'
        assert molpro['memory_mb'] == gaussian['memory_mb'] == 64000
        assert molpro['use_all_node_memory'] is True
        assert gaussian['use_all_node_memory'] is True
        assert molpro['cores'] == 4
        assert gaussian['cores'] == 12
        assert _molpro_stack_mw(molpro) >= molpro['min_stack_mw']
        script = (run_dir / 'tasks' / 'l3_geometry' / 'job.slurm').read_text()
        assert '#SBATCH --ntasks=4\n' in script
        assert '#SBATCH --cpus-per-task=1\n' in script
        assert '#SBATCH --mem=0\n' in script


def test_auto_molpro_cores_fail_when_node_cannot_meet_stack_minimum():
    spec = ch4_spec()
    spec['limits'] = {'max_nodes': 1}
    spec['tasks'] = [dict(spec['tasks'][1], geometry_from='initial')]
    spec['tasks'][0]['resources'] = {'walltime': '24:00:00',
                                     'min_stack_mw': 512}
    node_groups = [{'name': 'day', 'default': True, 'cores': 16,
                    'memory_mb': 3000, 'seconds': 86400}]
    with TemporaryDirectory() as temporary, \
            patch('kinbot.anl.site._partitions', return_value=node_groups):
        root = Path(temporary)
        with pytest.raises(RuntimeError, match='cannot support one core'):
            prepare(_write_spec(root, spec), root / 'run')


def test_auto_molpro_skips_faster_partition_without_enough_rank_memory():
    spec = ch4_spec(auto_resources=True)
    spec['tasks'] = [dict(spec['tasks'][1], geometry_from='initial')]
    node_groups = [
        {'name': 'short', 'default': True, 'cores': 32,
         'memory_mb': 3000, 'seconds': 86400},
        {'name': 'day', 'default': False, 'cores': 96,
         'memory_mb': 64000, 'seconds': 172800},
    ]
    with TemporaryDirectory() as temporary, \
            patch('kinbot.anl.site._partitions', return_value=node_groups):
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        resolved = json.loads((run_dir / 'workflow.json').read_text())
        assert resolved['tasks'][0]['resources']['partition'] == 'day'
        assert resolved['tasks'][0]['resources']['cores'] == 4


def test_auto_resource_ceiling_allows_a_lower_method_specific_rank_cap():
    spec = ch4_spec(auto_resources=True)
    assert spec['limits'] == {'max_nodes': 3}
    spec['tasks'][1]['resources']['max_cores'] = 8
    node_groups = [{'name': 'day', 'default': True, 'cores': 96,
                    'memory_mb': 512000, 'seconds': 86400}]
    with TemporaryDirectory() as temporary, \
            patch('kinbot.anl.site._partitions', return_value=node_groups):
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        resolved = json.loads((run_dir / 'workflow.json').read_text())
        by_id = {task['id']: task for task in resolved['tasks']}
        assert by_id['l3_geometry']['resources']['cores'] == 8
        assert by_id['harmonic']['resources']['cores'] == 16
        assert by_id['harmonic']['resources']['use_all_node_memory'] is True
        script = (run_dir / 'tasks' / 'l2_geometry' / 'job.slurm').read_text()
        assert '#SBATCH --mem=0\n' in script


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


def test_failed_job_missing_from_squeue_can_be_reconciled_and_retried():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        state_file = run_dir / 'state.json'
        state = json.loads(state_file.read_text())
        entry = state['tasks']['l2_geometry']
        entry.update(status='submitted', job_id='52301030')
        state_file.write_text(json.dumps(state))
        (run_dir / 'tasks' / 'l2_geometry' / 'execution.json').write_text(json.dumps({
            'schema': 1, 'task_id': 'l2_geometry',
            'geometry_sha256': entry['geometry_sha256'], 'status': 'failed',
            'error': 'MPI initialization failed',
        }))
        with patch('kinbot.anl.dispatch.subprocess.run', return_value=
                   subprocess.CompletedProcess([], 1, '',
                       'slurm_load_jobs error: Unable to contact slurm controller')):
            with pytest.raises(RuntimeError, match='squeue failed'):
                advance(run_dir)
        assert json.loads(state_file.read_text())['tasks']['l2_geometry']['status'] == 'submitted'
        invalid_id = subprocess.CompletedProcess(
            [], 1, '', 'slurm_load_jobs error: Invalid job id specified')
        with patch('kinbot.anl.dispatch.subprocess.run', return_value=invalid_id):
            assert advance(run_dir)['tasks']['l2_geometry']['status'] == 'failed'
            archive = retry_failed(run_dir, 'l2_geometry')
        assert (archive / 'execution.json').is_file()
        assert (run_dir / 'tasks' / 'l2_geometry' / 'job.slurm').is_file()
        assert json.loads(state_file.read_text())['tasks']['l2_geometry']['status'] == 'staged'


def test_status_refreshes_finished_job_without_submitting(capsys):
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        state_file = run_dir / 'state.json'
        state = json.loads(state_file.read_text())
        state['tasks']['l2_geometry'].update(status='submitted', job_id='1234')
        state_file.write_text(json.dumps(state))
        _mock_execution(run_dir, 'l2_geometry', geometry='final.xyz')
        with patch('kinbot.anl.dispatch._job_active', return_value=False), \
                patch('kinbot.anl.dispatch.subprocess.run',
                      side_effect=AssertionError('status must not submit')):
            main(['status', str(run_dir)])
        summary = json.loads(capsys.readouterr().out)
        assert summary['l2_geometry'] == 'complete'
        assert summary['l3_geometry'] == 'staged'
        assert json.loads(state_file.read_text())['tasks']['l2_geometry']['status'] == 'complete'


def test_preflight_reports_missing_program_before_submission():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        spec = ch4_spec()
        spec['tasks'][0]['profile']['command'] = str(root / 'missing_g16')
        run_dir = prepare(_write_spec(root, spec), root / 'run')
        with patch('kinbot.anl.dispatch.shutil.which', return_value='/bin/true'):
            with pytest.raises(RuntimeError, match='Missing executable'):
                preflight(run_dir)


def test_preflight_reports_silent_site_setup_exit_status():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        run_dir = prepare(_write_spec(root, ch4_spec()), root / 'run')
        (run_dir / 'site_setup.sh').write_text('exit 7\n')
        with patch('kinbot.anl.dispatch.shutil.which', return_value='/bin/true'):
            with pytest.raises(RuntimeError, match='exit status 7 without diagnostics'):
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
        (gaussian / 'bsd' / 'g16.profile').write_text(
            'false  # vendor profiles may use an unsuccessful probe\n'
            'export PROFILE_SOURCED=yes\n'
            'export PROFILE_SCRATCH="$GAUSS_SCRDIR"\n'
            'export GAUSS_SCRDIR=/profile-default\n')
        (root / 'cfour' / 'basis').mkdir()
        genbas = root / 'cfour' / 'basis' / 'GENBAS'
        genbas.write_text('basis fixture\n')
        (root / 'cfour' / 'payload').mkdir()
        (cfour / 'xcfour').replace(root / 'cfour' / 'payload' / 'xcfour')
        (cfour / 'xcfour').symlink_to(root / 'cfour' / 'payload' / 'xcfour')
        spec = ch4_spec()
        for task in spec['tasks']:
            task['resources']['partition'] = 'chosen_by_user'
        blocked_scratch = root / 'not-a-directory'
        blocked_scratch.write_text('blocks mkdir\n')
        node_scratch = root / 'node-scratch'
        site_scratch = root / 'site-scratch'
        base_path = os.environ['PATH']
        environment = {'PATH': f'{gaussian}:{molpro}:{cfour}:{base_path}',
                       'LOADEDMODULES': 'molpro/molpro24:cfour/2.1',
                       'CFOUR_GENBAS': '',
                       'GAUSS_SCRDIR': str(blocked_scratch),
                       'SLURM_TMPDIR': str(node_scratch),
                       'SCRATCH': str(site_scratch)}
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
                    'printf "|%s|%s|%s|%s|%s|%s" "${PROFILE_SOURCED:-}" '
                    '"${g16root:-}" "${CFOUR_GENBAS:-}" "${LOADED_MODULE:-}" '
                    '"${GAUSS_SCRDIR:-}" "${PROFILE_SCRATCH:-}"\n'
                )
                child_env = os.environ.copy()
                child_env['PATH'] = base_path
                result = subprocess.run(['bash', '-c', command], cwd=run_dir,
                                        env=child_env, capture_output=True,
                                        text=True, check=True)
                if backend == 'gaussian':
                    assert result.stdout.startswith(f'{gaussian / "g16"}\n')
                    assert f'|yes|{root / "gaussian"}|' in result.stdout
                    assert result.stdout.endswith(f'|{node_scratch}|{node_scratch}')
                    gaussian_command = command
                elif backend == 'cfour':
                    assert result.stdout.startswith(f'{cfour / "xcfour"}\n')
                    assert f'|{genbas}|cfour/2.1' in result.stdout
                    cfour_command = command
                else:
                    assert result.stdout.startswith(f'{molpro / "molpro"}\n')
                    assert '|molpro/molpro24|' in result.stdout
            override = root / 'different-GENBAS'
            override.write_text('override fixture\n')
            child_env['CFOUR_GENBAS'] = str(override)
            result = subprocess.run(['bash', '-c', cfour_command],
                                    cwd=run_dir, env=child_env, capture_output=True,
                                    text=True, check=True)
            assert f'|{override}|cfour/2.1' in result.stdout
            child_env['SLURM_TMPDIR'] = str(blocked_scratch)
            result = subprocess.run(['bash', '-c', gaussian_command],
                                    cwd=run_dir, env=child_env, capture_output=True,
                                    text=True, check=True)
            assert result.stdout.endswith(f'|{site_scratch}|{site_scratch}')
            explicit_scratch = root / 'explicit-scratch'
            child_env['GAUSS_SCRDIR'] = str(explicit_scratch)
            result = subprocess.run(['bash', '-c', gaussian_command],
                                    cwd=run_dir, env=child_env, capture_output=True,
                                    text=True, check=True)
            assert result.stdout.endswith(f'|{explicit_scratch}|{explicit_scratch}')
            (gaussian / 'bsd' / 'g16.profile').write_text('false\n')
            result = subprocess.run(['bash', '-c', gaussian_command],
                                    cwd=run_dir, env=child_env, capture_output=True,
                                    text=True, check=False)
            assert result.returncode == 1
            assert 'Gaussian profile failed with exit status 1' in result.stderr


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
