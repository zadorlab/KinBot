"""General dispatch behavior, with CH4 as the sole chemistry fixture."""

import json
from pathlib import Path
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
    _geometry_hash, _molpro_total_mw, _runtime_profile, advance, prepare, run_task,
    validate_spec,
)


def _write_spec(directory, spec):
    path = Path(directory) / 'spec.json'
    path.write_text(json.dumps(spec))
    return path


def _mock_execution(run_dir, ident, *, geometry=None):
    task_dir = Path(run_dir) / 'tasks' / ident
    task = json.loads((task_dir / 'task.json').read_text())
    if geometry:
        write(task_dir / geometry, Atoms(
            symbols=task['molecule']['symbols'],
            positions=read(task_dir / 'geometry.xyz').positions))
        details = {'final_geometry_sha256': _geometry_hash(read(task_dir / geometry))}
    else:
        details = {}
    (task_dir / 'execution.json').write_text(json.dumps({
        'task_id': ident, 'geometry_sha256': task['geometry_sha256'],
        'status': 'executed', 'details': details,
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
        molpro_geom = (run_dir / 'tasks' / 'l3_geometry' /
                       'l3_geometry.inp').read_text()
        assert 'optg,savexyz=l3_geometry.xyz' in molpro_geom
        assert 'C' in molpro_geom
        assert 'memory,' not in molpro_geom.lower()
        assert 'orient,noorient' in molpro_geom

        _mock_execution(run_dir, 'l3_geometry', geometry='l3_geometry.xyz')
        state = advance(run_dir)
        assert all(state['tasks'][task['id']]['status'] == 'staged'
                   for task in spec['tasks'][2:])
        for task in spec['tasks']:
            script = (run_dir / 'tasks' / task['id'] / 'job.slurm').read_text()
            assert script.count('#SBATCH --exclusive') == 1
            assert script.index('#SBATCH --exclusive') < script.index('set -euo pipefail')
        f12 = (run_dir / 'tasks' / 'f12_tz' / 'f12_tz.inp').read_text()
        assert 'ccsd(t)-f12,scale_trip=1' in f12
        assert 'kb_f12b=energy(2)' in f12
        cfour = (run_dir / 'tasks' / 'cfour_dboc' / 'ZMAT').read_text()
        assert 'COORDINATES=CARTESIAN' in cfour
        assert 'DBOC=ON' in cfour
        assert 'MEM_UNIT=MB,MEMORY_SIZE=11200' in cfour
        mrcc = (run_dir / 'tasks' / 'mrcc_ccsdtq' / 'MINP').read_text()
        assert 'calc=CCSDT(Q)' in mrcc
        assert 'mem=22400MB' in mrcc
        assert 'geom=xyz\n5\n' in mrcc
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
        _mock_execution(run_dir, 'l3_geometry', geometry='l3_geometry.xyz')
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
    spec['tasks'][0]['optimizer']['sella_kwargs'] = {'internal': False}
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
    spec['tasks'][1]['stdout'] = 'l3_geometry.out'
    with pytest.raises(ValueError, match='collide|Molpro'):
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
