"""Native-output components require complete, hash-checked task records."""

import hashlib
import json
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
from ase.io import read, write
from ase.units import Hartree, invcm, kJ, mol

from kinbot.anl.dispatch import _geometry_hash, advance, prepare
from kinbot.anl.model import IncompleteRecipeError
from kinbot.anl.recipes import recipe
from kinbot.anl.results import parse_result
from kinbot.anl.workflow import cbs_task_component, task_component
from tests.anl_fixture import dispatch_spec, molpro_task


_F12 = """basis=cc-pVTZ-F12
ccsd(t)-f12,scale_trip=1
 !CCSD(T)-F12b total energy -40.454906199189
 CCSD(T)-F12/cc-pVTZ-F12 energy= -40.454906199189
 Molpro calculation terminated
"""

_F12_QZ = """basis=cc-pVQZ-F12
ccsd(t)-f12,scale_trip=1
 !CCSD(T)-F12b total energy -40.456608306474
 CCSD(T)-F12/cc-pVQZ-F12 energy= -40.456608306474
 Molpro calculation terminated
"""


def _complete_task(run_dir, task, output):
    task_id = task['id']
    directory = run_dir / 'tasks' / task_id
    output_name = task['result_parser']['file']
    (directory / output_name).write_text(output)
    stdout_name = task.get('stdout', 'stdout.txt')
    stderr_name = task.get('stderr', 'stderr.txt')
    for name in (stdout_name, stderr_name):
        if name != output_name:
            (directory / name).write_text('')
    state_path = run_dir / 'state.json'
    state = json.loads(state_path.read_text())
    state['tasks'][task_id]['status'] = 'complete'
    state_path.write_text(json.dumps(state))
    artifacts = {
        name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
        for name in {'geometry.xyz', 'task.json', task['input_name'],
                     output_name, stdout_name, stderr_name}
    }
    (directory / 'execution.json').write_text(json.dumps({
        'schema': 1, 'task_id': task_id, 'status': 'executed',
        'geometry_sha256': state['tasks'][task_id]['geometry_sha256'],
        'artifacts': artifacts,
        'details': {'parsed_result': parse_result(output, task['result_parser'])},
    }))
    return directory


def _completed_f12(root):
    spec = dispatch_spec()
    geometry = next(task for task in spec['tasks']
                    if task['id'] == 'l3_geometry')
    geometry['geometry_from'] = 'initial'
    task = next(task for task in spec['tasks'] if task['id'] == 'f12_tz')
    spec['tasks'] = [geometry, task]
    spec_file = root / 'spec.json'
    spec_file.write_text(json.dumps(spec))
    run_dir = prepare(spec_file, root / 'run')
    _complete_optimization(run_dir, geometry)
    advance(run_dir, submit=False)
    directory = _complete_task(run_dir, task, _F12)
    return run_dir, directory


def _completed_f12_pair(root, *, split_geometry=False):
    spec = dispatch_spec()
    geometry = next(task for task in spec['tasks']
                    if task['id'] == 'l3_geometry')
    geometry['geometry_from'] = 'initial'
    tasks = [task for task in spec['tasks']
             if task['id'] in ('f12_tz', 'f12_qz')]
    optimizations = [geometry]
    if split_geometry:
        other = deepcopy(geometry)
        other['id'] = 'l3_other'
        tasks[1]['geometry_from'] = other['id']
        optimizations.append(other)
    spec['tasks'] = [*optimizations, *tasks]
    spec_file = root / 'spec.json'
    spec_file.write_text(json.dumps(spec))
    run_dir = prepare(spec_file, root / 'run')
    for optimization in optimizations:
        _complete_optimization(run_dir, optimization)
    advance(run_dir, submit=False)
    for task, output in zip(tasks, (_F12, _F12_QZ)):
        _complete_task(run_dir, task, output)
    return run_dir


def _complete_optimization(run_dir, task, *, move=0.0):
    ident = task['id']
    directory = run_dir / 'tasks' / ident
    atoms = read(directory / 'geometry.xyz')
    atoms.positions[1, 0] += move
    write(directory / 'final.xyz', atoms)
    (directory / 'optimization.log').write_text('Synthetic accepted optimization\n')
    generated = []
    if task['profile']['calculator'] == 'molpro':
        generated = [f'{ident}_step_0001{suffix}'
                     for suffix in ('.inp', '.out', '.log', '.xyz')]
        for name in generated:
            (directory / name).write_text('Synthetic native artifact\n')
    artifact_names = {'geometry.xyz', 'task.json', 'final.xyz',
                      'optimization.log', *generated}
    artifacts = {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                 for name in artifact_names}
    geometry_hash = _geometry_hash(atoms)
    (directory / 'execution.json').write_text(json.dumps({
        'schema': 1, 'task_id': ident, 'status': 'executed',
        'geometry_sha256': _geometry_hash(read(directory / 'geometry.xyz')),
        'artifacts': artifacts,
        'details': {'final_geometry_sha256': geometry_hash,
                    'generated_files': generated},
    }))
    state_path = run_dir / 'state.json'
    state = json.loads(state_path.read_text())
    state['tasks'][ident]['status'] = 'complete'
    state['tasks'][ident]['final_geometry_sha256'] = geometry_hash
    state_path.write_text(json.dumps(state))
    return geometry_hash


def _harmonic_output(basis, shift):
    modes = [1343.29, 1343.90, 1344.76, 1570.56, 1570.88,
             3033.49, 3151.55, 3152.19, 3153.44]
    modes = [value + shift for value in modes]
    zpe_cm = sum(modes) / 2
    zpe_hartree = zpe_cm * invcm / Hartree
    zpe_kj = zpe_hartree * Hartree * mol / kJ
    lines = '\n'.join(f' {i} {value:.2f}' for i, value in enumerate(modes, 1))
    return (f'basis={basis}\nccsd(t)\nfrequencies,numerical\n'
            'PROGRAM * FREQUENCIES (Calculation of harmonic vibrational '
            'spectra for CCSD(T))\nVibration Wavenumber\n'
            f'{lines}\n\nZero point energy: {zpe_hartree:.8f} [H] '
            f'{zpe_cm:.2f} [1/CM] {zpe_kj:.2f} [KJ/MOL]\n'
            'Molpro calculation terminated\n')


def _vpt2_output(basis, correction):
    harmonic = 9800.0
    return (f'#p B2PLYP/{basis} Freq=Anharmonic '
            'EmpiricalDispersion=GD3BJ\n'
            'Anharmonic Zero Point Energy\n'
            f'Harmonic       : cm-1 = {harmonic:.5f} ;\n'
            f'Anharmonic Pot.: cm-1 = {correction:.5f} ;\n'
            'Watson+Coriolis: cm-1 = 0.00000 ;\n'
            f'Total Anharm   : cm-1 = {harmonic + correction:.5f} ;\n'
            'Normal termination of Gaussian 16\n')


def _completed_zpe_pair(root, *, kind):
    spec = dispatch_spec()
    is_harmonic = kind == 'harmonic'
    source = next(task for task in spec['tasks'] if task['id'] == (
        'l3_geometry' if is_harmonic else 'l2_geometry'))
    optimizations, evaluations = [], []
    for suffix, basis in (('tz', 'cc-pVTZ'), ('qz', 'cc-pVQZ')):
        opt = deepcopy(source)
        opt['id'] = f'geom_{suffix}'
        opt['geometry_from'] = 'initial'
        opt['profile']['basis'] = basis
        optimizations.append(opt)
        if is_harmonic:
            task = molpro_task(
                f'freq_{suffix}',
                f'basis={basis}\nhf\nccsd(t)\nfrequencies,numerical\n',
                geometry_from=opt['id'],
                result_parser={'kind': 'molpro_harmonic', 'basis': basis})
        else:
            task = deepcopy(next(item for item in spec['tasks']
                                 if item['id'] == 'gaussian_vpt2'))
            task['id'] = f'vpt2_{suffix}'
            task['geometry_from'] = opt['id']
            task['depends_on'] = []
            task['input_name'] = f"{task['id']}.com"
            task['stdin'] = task['input_name']
            task['stdout'] = f"{task['id']}.log"
            task['stderr'] = f"{task['id']}.err"
            task['required_outputs'] = [task['stdout']]
            task['success_marker']['file'] = task['stdout']
            task['result_parser']['file'] = task['stdout']
            task['result_parser']['basis'] = basis
            task['input_template'] = task['input_template'].replace(
                'B2PLYP/cc-pVTZ', f'B2PLYP/{basis}')
        evaluations.append(task)
    spec['tasks'] = optimizations + evaluations
    spec_file = root / 'spec.json'
    spec_file.write_text(json.dumps(spec))
    run_dir = prepare(spec_file, root / 'run')
    lower_geometry = _complete_optimization(run_dir, optimizations[0])
    upper_geometry = _complete_optimization(run_dir, optimizations[1], move=0.01)
    advance(run_dir, submit=False)
    for index, task in enumerate(evaluations):
        basis = ('cc-pVTZ', 'cc-pVQZ')[index]
        output = (_harmonic_output(basis, index)
                  if is_harmonic else _vpt2_output(basis, -137.0 + index))
        _complete_task(run_dir, task, output)
    return run_dir, lower_geometry, upper_geometry


def test_read_completed_f12_component_reparses_and_checks_artifacts():
    with TemporaryDirectory() as temporary:
        run_dir, directory = _completed_f12(Path(temporary))
        result = task_component(run_dir, 'f12_tz', key='f12_tz',
                                state_id='state-A')
        assert result.value_hartree == pytest.approx(-40.454906199189)
        assert result.method == 'CCSD(T)-F12b'
        assert result.settings['scale_trip'] == 1
        assert result.source_sha256 == hashlib.sha256(
            (directory / 'f12_tz.out').read_bytes()).hexdigest()

        (directory / 'f12_tz.out').write_text(_F12.replace('-40.454906199189',
                                                          '-40.454000000000'))
        with pytest.raises(RuntimeError, match='artifact'):
            task_component(run_dir, 'f12_tz', key='f12_tz',
                           state_id='state-A')


def test_old_unparsed_dispatch_result_cannot_enter_recipe():
    with TemporaryDirectory() as temporary:
        run_dir, directory = _completed_f12(Path(temporary))
        execution_path = directory / 'execution.json'
        execution = json.loads(execution_path.read_text())
        del execution['details']['parsed_result']
        execution_path.write_text(json.dumps(execution))
        with pytest.raises(ValueError, match='differs from native'):
            task_component(run_dir, 'f12_tz', key='reference_cbs',
                           state_id='state-A')
        state_path = run_dir / 'state.json'
        state = json.loads(state_path.read_text())
        state['tasks']['f12_tz']['status'] = 'submitted'
        state_path.write_text(json.dumps(state))
        with pytest.raises(IncompleteRecipeError, match='not complete'):
            task_component(run_dir, 'f12_tz', key='reference_cbs',
                           state_id='state-A')


def test_verified_f12_task_pair_builds_recipe_cbs_component():
    with TemporaryDirectory() as temporary:
        run_dir = _completed_f12_pair(Path(temporary))
        requirement = next(item for item in recipe('ANL0-F12').requirements
                           if item.key == 'reference_cbs')
        kwargs = {'requirement': requirement, 'state_id': 'state-A',
                  'lower_basis': 'cc-pVTZ-F12',
                  'upper_basis': 'cc-pVQZ-F12'}
        result = cbs_task_component(run_dir, 'f12_tz', 'f12_qz', **kwargs)
        assert result.value_hartree == pytest.approx(-40.457504545051066)
        assert result.backend == 'composite'
        assert result.basis == 'CBS(cc-pVTZ-F12,cc-pVQZ-F12)'
        assert result.settings['scale_trip'] == 1
        assert result.settings['upper_cardinal'] == 4
        assert len(result.source_sha256) == 64

        with pytest.raises(ValueError, match='basis order'):
            cbs_task_component(run_dir, 'f12_qz', 'f12_tz', **kwargs)
        with pytest.raises(ValueError, match='distinct'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_tz', **kwargs)
        with pytest.raises(ValueError, match='basis pair'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz',
                               **{**kwargs, 'upper_basis': 'cc-pV5Z-F12'})
        harmonic = next(item for item in recipe('ANL1').requirements
                        if item.key == 'harmonic_zpe')
        with pytest.raises(ValueError, match='geometry level'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz',
                               **{**kwargs, 'requirement': harmonic})


def test_cbs_pair_rejects_incomplete_or_changed_native_output():
    with TemporaryDirectory() as temporary:
        run_dir = _completed_f12_pair(Path(temporary))
        requirement = next(item for item in recipe('ANL0-F12').requirements
                           if item.key == 'reference_cbs')
        kwargs = {'requirement': requirement, 'state_id': 'state-A',
                  'lower_basis': 'cc-pVTZ-F12',
                  'upper_basis': 'cc-pVQZ-F12'}
        state_path = run_dir / 'state.json'
        state = json.loads(state_path.read_text())
        state['tasks']['f12_qz']['status'] = 'submitted'
        state_path.write_text(json.dumps(state))
        with pytest.raises(IncompleteRecipeError, match='not complete'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz', **kwargs)
        state['tasks']['f12_qz']['status'] = 'complete'
        state_path.write_text(json.dumps(state))
        output = run_dir / 'tasks' / 'f12_qz' / 'f12_qz.out'
        output.write_text(output.read_text().replace('-40.456608306474',
                                                    '-40.456000000000'))
        with pytest.raises(RuntimeError, match='artifact'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz', **kwargs)


def test_electronic_cbs_rejects_distinct_optimizers_even_at_same_coordinates():
    with TemporaryDirectory() as temporary:
        run_dir = _completed_f12_pair(Path(temporary), split_geometry=True)
        requirement = next(item for item in recipe('ANL0-F12').requirements
                           if item.key == 'reference_cbs')
        with pytest.raises(ValueError, match='different geometry sources'):
            cbs_task_component(
                run_dir, 'f12_tz', 'f12_qz', requirement=requirement,
                state_id='state-A', lower_basis='cc-pVTZ-F12',
                upper_basis='cc-pVQZ-F12')


def test_harmonic_cbs_accepts_two_verified_basis_optimized_geometries():
    with TemporaryDirectory() as temporary:
        run_dir, lower_geometry, upper_geometry = _completed_zpe_pair(
            Path(temporary), kind='harmonic')
        assert lower_geometry != upper_geometry
        requirement = next(item for item in recipe('ANL1').requirements
                           if item.key == 'harmonic_zpe')
        result = cbs_task_component(
            run_dir, 'freq_tz', 'freq_qz', requirement=requirement,
            state_id='state-A', lower_basis='cc-pVTZ',
            upper_basis='cc-pVQZ')
        lower = task_component(run_dir, 'freq_tz', key='freq_tz',
                               state_id='state-A')
        upper = task_component(run_dir, 'freq_qz', key='freq_qz',
                               state_id='state-A')
        assert result.value_hartree == pytest.approx(
            upper.value_hartree + 0.5265464668174269 *
            (upper.value_hartree - lower.value_hartree))
        assert result.quantity == 'zpe'
        assert result.geometry_sha256 == upper_geometry
        assert result.settings['geometry_mode'] == 'basis_optimized'
        assert result.value_hartree > 0

        final = run_dir / 'tasks' / 'geom_tz' / 'final.xyz'
        final.write_text(final.read_text() + '\n')
        with pytest.raises(RuntimeError, match='artifact'):
            cbs_task_component(run_dir, 'freq_tz', 'freq_qz',
                               requirement=requirement, state_id='state-A',
                               lower_basis='cc-pVTZ', upper_basis='cc-pVQZ')


def test_profiled_vpt2_cbs_accepts_two_matching_gaussian_geometries():
    with TemporaryDirectory() as temporary:
        run_dir, lower_geometry, upper_geometry = _completed_zpe_pair(
            Path(temporary), kind='vpt2')
        assert lower_geometry != upper_geometry
        equation = recipe('ANL1', vpt2_method='B2PLYP-D3BJ', vpt2_cbs=True)
        assert equation.name == 'profiled:ANL1:B2PLYP-D3BJ:VPT2-CBS'
        requirement = next(item for item in equation.requirements
                           if item.key == 'vpt2_correction')
        result = cbs_task_component(
            run_dir, 'vpt2_tz', 'vpt2_qz', requirement=requirement,
            state_id='state-A', lower_basis='cc-pVTZ',
            upper_basis='cc-pVQZ')
        expected_cm = -136.0 + 0.5265464668174269
        assert result.value_hartree == pytest.approx(expected_cm * invcm / Hartree)
        assert result.quantity == 'correction'
        assert result.geometry_sha256 == upper_geometry
        assert result.settings['dispersion'] == 'GD3BJ'
