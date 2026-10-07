"""The external-site validation graph is general and cannot overclaim ANL."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
from ase import Atoms
from ase.db import connect
from ase.io import read, write

from kinbot.anl.dispatch import (_geometry_hash, prepare,
                                 validate_spec)
from kinbot.anl.model import ComponentResult
from kinbot.anl.recipes import recipe
from kinbot.anl import validation as validation_module
from kinbot.anl.validation import (common_corrections_validation_spec,
                                   assemble_profiled_anl0_f12,
                                   audit_anl0_post_geometry_run,
                                   audit_composite_run,
                                   audit_current_base_run,
                                   audit_higher_order_run,
                                   audit_interface_run,
                                   audit_post_geometry_run,
                                   current_base_validation_spec,
                                   composite_validation_spec,
                                   interface_validation_spec,
                                   higher_order_validation_spec,
                                   audit_kinbot_run,
                                   main as validation_main,
                                   molecule_from_completed_run,
                                   molecule_from_database,
                                   molecule_from_smiles,
                                   post_geometry_validation_spec)
from kinbot.anl import site
from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.stationary_pt import StationaryPoint


def _molecule():
    return {'symbols': ['C', 'H', 'H', 'H', 'H'],
            'positions': [[0., 0., 0.], [0.63, 0.63, 0.63],
                          [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63],
                          [0.63, -0.63, -0.63]],
            'charge': 0, 'multiplicity': 1}


def test_non_mrcc_graph_has_geometry_barrier_and_parallel_fanout():
    spec = interface_validation_spec(
        _molecule(), max_nodes=4, partition='day-long-cpu')
    resolved = deepcopy(spec)
    for task in resolved['tasks']:
        task['resources'].update(cores=4, memory_mb=64000,
                                 partition='test')
    validate_spec(resolved)
    tasks = {task['id']: task for task in spec['tasks']}
    assert tasks['l3_geometry']['geometry_from'] == 'l2_geometry'
    for ident in ('harmonic', 'f12_tz', 'f12_qz', 'ccsdt_dz',
                  'cfour_dboc'):
        assert tasks[ident]['geometry_from'] == 'l3_geometry'
    assert tasks['gaussian_vpt2']['geometry_from'] == 'l2_geometry'
    assert tasks['gaussian_vpt2']['depends_on'] == ['l3_geometry']
    assert spec['intent']['mrcc_enabled'] is False
    assert spec['intent']['claim'] == 'interface-validation-only'
    assert all(task['resources']['cores'] == 'auto' for task in spec['tasks'])
    assert all(task['resources']['partition'] == 'day-long-cpu'
               for task in spec['tasks'])


def test_current_base_reuses_l2_but_reruns_unrestricted_l3_and_base_terms():
    molecule = deepcopy(_molecule())
    molecule['source'] = {
        'run_dir': '/accepted/interface', 'task_id': 'l2_geometry',
        'geometry_sha256': '2' * 64, 'artifact_sha256': 'a' * 64,
        'profile': {'calculator': 'gaussian', 'method': 'B2PLYP',
                    'basis': 'cc-pVTZ'},
    }
    spec = current_base_validation_spec(
        molecule, max_nodes=3, partition='day-long-cpu')
    tasks = {task['id']: task for task in spec['tasks']}
    assert set(tasks) == {
        'l3_geometry', 'harmonic', 'f12_tz', 'f12_qz', 'cfour_dboc'}
    assert tasks['l3_geometry']['geometry_from'] == 'initial'
    assert tasks['l3_geometry']['profile']['method'] == 'CCSD(T)'
    for ident in ('harmonic', 'f12_tz', 'f12_qz'):
        assert tasks[ident]['geometry_from'] == 'l3_geometry'
        assert tasks[ident]['result_parser']['reference'] == 'RHF'
        assert 'uccsd(t)' in tasks[ident]['input_template'].lower()
        assert 'uhf_uccsd' not in tasks[ident]['input_template'].lower()
    resolved = deepcopy(spec)
    for task in resolved['tasks']:
        task['resources'].update(cores=4, memory_mb=64000,
                                 partition='test')
    validate_spec(resolved)


def test_legacy_interface_audit_is_readable_but_recipe_incompatible(monkeypatch):
    spec = interface_validation_spec(_molecule())
    for task in spec['tasks']:
        parser = task.get('result_parser', {})
        if parser.get('kind') in ('molpro_energy', 'molpro_harmonic'):
            parser.pop('reference', None)
            task['input_template'] = (task['input_template']
                                      .replace('rhf\n', 'hf\n')
                                      .replace('uccsd(t)', 'ccsd(t)')
                                      .replace('uccsd(t)-f12b',
                                               'ccsd(t)-f12'))
    state = {'tasks': {task['id']: {'status': 'complete'}
                       for task in spec['tasks']}}
    parsed = {
        'harmonic': {'kind': 'molpro_harmonic'},
        'f12_tz': {'kind': 'molpro_energy', 'energy_hartree': -10.0},
        'f12_qz': {'kind': 'molpro_energy', 'energy_hartree': -10.1},
        'ccsdt_dz': {'kind': 'molpro_energy', 'energy_hartree': -9.9},
        'cfour_dboc': {'kind': 'cfour_dboc'},
        'gaussian_vpt2': {'kind': 'gaussian_vpt2'},
    }
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (Path(run_dir), spec, state))
    monkeypatch.setattr(
        validation_module, '_verified_task_result',
        lambda run_dir, task_id: (None, None, None, None, None,
                                  parsed[task_id]))
    result = audit_interface_run('/legacy/interface')
    assert result['status'] == 'interface_complete_legacy_recipe_incompatible'
    assert set(result['legacy_molpro_tasks']) == {
        'harmonic', 'f12_tz', 'f12_qz', 'ccsdt_dz'}
    assert 'verified_f12_cbs_hartree' not in result


def test_higher_order_probe_routes_closed_and_open_shell_without_pair_locking():
    closed = higher_order_validation_spec(
        _molecule(), max_nodes=3, partition='day-long-cpu')
    closed_tasks = {task['id']: task for task in closed['tasks']}
    assert closed_tasks['ccsdtq_tz']['backend'] == 'mrcc'
    assert closed_tasks['ccsdtq_dz']['backend'] == 'mrcc'
    assert closed_tasks['ccsdtqp_dz']['backend'] == 'mrcc'
    assert closed_tasks['ccsdtqp_dz']['result_parser']['reference'] == 'RHF'
    assert closed_tasks['ccsdtqp_dz']['result_parser']['correlation'] == \
        'unrestricted'
    assert closed_tasks['ccsdtqp_dz']['resources']['walltime'] == '7-00:00:00'
    assert closed_tasks['ccsdtqp_dz']['resources']['max_cores'] == 8
    assert closed_tasks['ccsdt_tz']['backend'] == 'molpro'
    for ident in ('ccsdtq_tz', 'ccsdtq_dz'):
        assert closed_tasks[ident]['result_parser']['reference'] == 'RHF'
        assert closed_tasks[ident]['result_parser']['correlation'] == 'unrestricted'
        assert closed_tasks[ident]['result_parser']['program'] == 'mrcc'
        assert 'scftype=RHF' in closed_tasks[ident]['input_template']
        assert 'ccprog=mrcc' in closed_tasks[ident]['input_template']
    assert closed['intent']['equation'] == (
        'CCSDT(Q)/TZ - CCSD(T)/TZ + CCSDTQ(P)/DZ - CCSDT(Q)/DZ')

    radical = deepcopy(_molecule())
    radical['symbols'] = radical['symbols'][:-1]
    radical['positions'] = radical['positions'][:-1]
    radical['multiplicity'] = 2
    opened = higher_order_validation_spec(radical)
    open_tasks = {task['id']: task for task in opened['tasks']}
    assert '\nuccsd(t)\n' in \
        open_tasks['ccsdt_tz']['input_template'].lower()
    assert 'uhf_uccsd' not in open_tasks['ccsdt_tz']['input_template'].lower()
    for ident in ('ccsdtq_tz', 'ccsdtq_dz', 'ccsdtqp_dz'):
        assert open_tasks[ident]['backend'] == 'mrcc'
        assert open_tasks[ident]['result_parser']['reference'] == 'ROHF'
        assert open_tasks[ident]['result_parser']['correlation'] == 'unrestricted'
        assert open_tasks[ident]['result_parser']['program'] == 'mrcc'
        assert 'scftype=ROHF' in open_tasks[ident]['input_template']
        assert 'ccprog=mrcc' in open_tasks[ident]['input_template']
        assert 'rohftype=semicanonical' in open_tasks[ident]['input_template']

    for spec in (closed, opened):
        resolved = deepcopy(spec)
        for task in resolved['tasks']:
            task['resources'].update(cores=4, memory_mb=64000,
                                     partition='test')
        validate_spec(resolved)


def test_higher_order_probe_can_select_a_minimal_literature_pair():
    spec = higher_order_validation_spec(
        _molecule(), task_ids=('ccsdtq_dz', 'ccsdtqp_dz'))
    assert [task['id'] for task in spec['tasks']] == [
        'ccsdtq_dz', 'ccsdtqp_dz']
    assert spec['intent']['selected_tasks'] == [
        'ccsdtq_dz', 'ccsdtqp_dz']
    with pytest.raises(ValueError, match='Invalid higher-order task'):
        higher_order_validation_spec(_molecule(), task_ids=())


def test_common_corrections_graph_pins_core_and_relativistic_differences():
    spec = common_corrections_validation_spec(
        _molecule(), max_nodes=4, partition='day-long-cpu')
    tasks = {task['id']: task for task in spec['tasks']}
    assert set(tasks) == {'cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz',
                          'rel_dkh', 'rel_nonrel'}
    for ident in ('cv_ae_tz', 'cv_ae_qz', 'rel_dkh', 'rel_nonrel'):
        assert tasks[ident]['result_parser']['core'] == 'all-electron'
        assert ';core}' in tasks[ident]['input_template'].lower()
    for ident in ('cv_fc_tz', 'cv_fc_qz'):
        assert tasks[ident]['result_parser']['core'] == 'frozen'
        assert ';core}' not in tasks[ident]['input_template'].lower()
    assert tasks['rel_dkh']['result_parser']['relativistic'] == 'DKH2'
    assert 'set,dkho=2' in tasks['rel_dkh']['input_template'].lower()
    assert tasks['rel_nonrel']['result_parser']['relativistic'] == 'none'
    resolved = deepcopy(spec)
    for task in resolved['tasks']:
        task['resources'].update(cores=4, memory_mb=64000,
                                 partition='test')
    validate_spec(resolved)


def test_post_geometry_graph_has_one_global_parallel_limit():
    spec = post_geometry_validation_spec(
        _molecule(), max_nodes=3, partition='day-long-cpu')
    assert spec['name'] == 'anl-post-geometry-validation'
    assert spec['limits'] == {'max_nodes': 3}
    tasks = {task['id']: task for task in spec['tasks']}
    assert set(tasks) == {
        'ccsdt_tz', 'ccsdt_dz', 'ccsdtq_tz', 'ccsdtq_dz',
        'ccsdtqp_dz', 'cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz',
        'rel_dkh', 'rel_nonrel'}
    assert all(task['geometry_from'] == 'initial' for task in tasks.values())
    resolved = deepcopy(spec)
    for task in resolved['tasks']:
        task['resources'].update(cores=4, memory_mb=64000,
                                 partition='test')
    validate_spec(resolved)


def test_post_geometry_graph_can_limit_execution_to_anl0():
    spec = post_geometry_validation_spec(
        _molecule(), max_nodes=2, partition='day-long-cpu', anl0_only=True)
    tasks = {task['id']: task for task in spec['tasks']}
    assert set(tasks) == {
        'ccsdt_dz', 'ccsdtq_dz',
        'cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz',
        'rel_dkh', 'rel_nonrel'}
    assert spec['intent']['selected_higher_order_tasks'] == [
        'ccsdt_dz', 'ccsdtq_dz']
    assert 'ccsdtq_tz' not in tasks
    assert 'ccsdtqp_dz' not in tasks


def test_composite_graph_releases_every_property_after_one_l3_barrier():
    spec = composite_validation_spec(
        _molecule(), max_nodes=4, partition='day-long-cpu', anl0_only=True)
    assert spec['name'] == 'anl-composite-validation'
    assert spec['limits'] == {'max_nodes': 4}
    tasks = {task['id']: task for task in spec['tasks']}
    assert set(tasks) == {
        'l2_geometry', 'l3_geometry', 'harmonic', 'f12_tz', 'f12_qz',
        'cfour_dboc', 'gaussian_vpt2', 'ccsdt_dz', 'ccsdtq_dz',
        'cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz',
        'rel_dkh', 'rel_nonrel'}
    assert tasks['l2_geometry']['geometry_from'] == 'initial'
    assert tasks['l3_geometry']['geometry_from'] == 'l2_geometry'
    assert tasks['gaussian_vpt2']['geometry_from'] == 'l2_geometry'
    assert tasks['gaussian_vpt2']['depends_on'] == ['l3_geometry']
    property_ids = set(tasks) - {'l2_geometry', 'l3_geometry',
                                 'gaussian_vpt2'}
    assert all(tasks[ident]['geometry_from'] == 'l3_geometry'
               for ident in property_ids)
    assert spec['intent']['anl0_only'] is True
    resolved = deepcopy(spec)
    for task in resolved['tasks']:
        task['resources'].update(cores=4, memory_mb=64000,
                                 partition='test')
    validate_spec(resolved)


def test_atomic_composite_omits_fictitious_vibrational_jobs():
    atom = molecule_from_smiles('[H]', charge=0, multiplicity=2)
    spec = composite_validation_spec(
        atom, max_nodes=8, partition='day-long-cpu', anl0_only=True)
    tasks = {task['id']: task for task in spec['tasks']}
    assert 'harmonic' not in tasks
    assert 'gaussian_vpt2' not in tasks
    assert {'l2_geometry', 'l3_geometry', 'f12_tz', 'f12_qz',
            'cfour_dboc', 'ccsdt_dz', 'cv_ae_tz', 'cv_ae_qz',
            'cv_fc_tz', 'cv_fc_qz', 'rel_dkh', 'rel_nonrel'} <= set(tasks)
    assert 'ccsdtq_dz' not in tasks
    for ident in ('f12_tz', 'f12_qz'):
        assert '\nuccsd-f12b\n' in tasks[ident]['input_template'].lower()
        assert 'scale_trip' not in tasks[ident]['input_template'].lower()
        assert tasks[ident]['result_parser']['rank_exact_electrons'] == 1
    resolved = deepcopy(spec)
    for task in resolved['tasks']:
        task['resources'].update(
            cores=1, memory_mb=64000, partition='test')
    validate_spec(resolved)


def test_two_electron_post_graph_skips_impossible_higher_rank_job():
    hydrogen = {
        'symbols': ['H', 'H'],
        'positions': [[0., 0., 0.], [0., 0., 0.74]],
        'charge': 0, 'multiplicity': 1,
    }
    spec = post_geometry_validation_spec(
        hydrogen, max_nodes=2, partition='day-long-cpu', anl0_only=True)
    tasks = {task['id']: task for task in spec['tasks']}
    assert 'ccsdt_dz' in tasks
    assert 'ccsdtq_dz' not in tasks
    assert spec['intent']['rank_exact_higher_order']['electron_count'] == 2
    assert spec['intent']['rank_exact_higher_order']['correction_hartree'] == 0.


def test_two_electron_post_audit_replaces_failed_high_job_by_rank_identity(
        monkeypatch):
    hydrogen = {
        'symbols': ['H', 'H'],
        'positions': [[0., 0., 0.], [0., 0., 0.74]],
        'charge': 0, 'multiplicity': 1,
        'source': {'geometry_sha256': '4' * 64},
    }
    spec = post_geometry_validation_spec(
        hydrogen, anl0_only=True)
    # Reproduce an archive prepared before rank-aware staging: the now
    # unnecessary MRCC task is present and failed, while every required task
    # completed.
    obsolete_high = higher_order_validation_spec(
        hydrogen, task_ids=('ccsdtq_dz',))['tasks'][0]
    spec['tasks'].append(obsolete_high)
    state = {'tasks': {task['id']: {'status': 'complete'}
                       for task in spec['tasks']}}
    state['tasks']['ccsdtq_dz']['status'] = 'failed'
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (Path(run_dir), spec, state))

    low = ComponentResult(
        key='hoe_low', value_hartree=-1.151234,
        quantity='electronic', method='CCSD(T)', basis='cc-pVDZ',
        backend='molpro', state_id='anl0-post-geometry-validation-state',
        charge=0, multiplicity=1, geometry_sha256='4' * 64,
        source_sha256='a' * 64, source='ccsdt_dz.out',
        settings={'reference': 'RHF', 'correlation': 'unrestricted'})

    def component(run_dir, task_id, *, key, state_id):
        assert task_id == 'ccsdt_dz'
        return low

    monkeypatch.setattr(validation_module, 'task_component', component)
    monkeypatch.setattr(
        validation_module, 'audit_common_corrections_run',
        lambda run_dir: {'status': 'common_corrections_interface_complete'})
    result = audit_anl0_post_geometry_run('/synthetic/hydrogen-post')
    assert result['higher_order_dz_hartree'] == 0.
    assert result['rank_exact_higher_order'] is True
    assert result['higher_order_components']['ccsdtq_dz']['backend'] == \
        'known_zero'
    assert result['other_task_statuses']['ccsdtq_dz'] == 'failed'


def test_post_geometry_audit_requires_and_combines_both_groups(monkeypatch):
    spec = post_geometry_validation_spec(_molecule())
    spec['molecule']['source'] = {'geometry_sha256': '4' * 64}
    state = {'tasks': {task['id']: {'status': 'complete'}
                       for task in spec['tasks']}}
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (Path(run_dir), spec, state))
    monkeypatch.setattr(
        validation_module, 'audit_higher_order_run',
        lambda run_dir: {'status': 'higher_order_interface_complete'})
    monkeypatch.setattr(
        validation_module, 'audit_common_corrections_run',
        lambda run_dir: {'status': 'common_corrections_interface_complete'})
    result = audit_post_geometry_run('/synthetic/post')
    assert result['status'] == 'post_geometry_interface_complete'
    assert result['geometry_sha256'] == '4' * 64
    assert len(result['task_statuses']) == 11


def test_composite_audit_combines_base_vpt2_and_anl0_fanout(monkeypatch):
    spec = composite_validation_spec(_molecule(), anl0_only=True)
    state = {'tasks': {task['id']: {'status': 'complete'}
                       for task in spec['tasks']}}
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (Path(run_dir), spec, state))
    monkeypatch.setattr(
        validation_module, '_verified_task_result',
        lambda run_dir, task_id: (
            None, None, None, None, None,
            {'review_required': True, 'warnings': ['native warning'],
             'anharmonic_correction_hartree': -0.001}))
    monkeypatch.setattr(
        validation_module, 'audit_current_base_run',
        lambda run_dir: {'status': 'current_base_interface_complete'})
    monkeypatch.setattr(
        validation_module, 'audit_anl0_post_geometry_run',
        lambda run_dir: {'status': 'anl0_post_geometry_interface_complete'})
    result = audit_composite_run('/synthetic/composite')
    assert result['status'] == 'composite_anl_interface_complete'
    assert result['requested_ladder_head'] == 'ANL0-F12'
    assert result['base']['status'] == 'current_base_interface_complete'
    assert result['post_geometry']['status'] == \
        'anl0_post_geometry_interface_complete'
    assert result['vpt2']['review_required'] is True


def test_anl0_post_geometry_audit_ignores_anl1_only_failures(monkeypatch):
    spec = post_geometry_validation_spec(_molecule())
    spec['molecule']['source'] = {'geometry_sha256': '4' * 64}
    state = {'tasks': {task['id']: {'status': 'complete'}
                       for task in spec['tasks']}}
    state['tasks']['ccsdtq_tz']['status'] = 'failed'
    state['tasks']['ccsdtqp_dz']['status'] = 'failed'
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (Path(run_dir), spec, state))

    def component(run_dir, task_id, *, key, state_id):
        value = -9.91 if task_id == 'ccsdtq_dz' else -9.90
        return ComponentResult(
            key=key, value_hartree=value, quantity='electronic',
            method='CCSDT(Q)' if task_id == 'ccsdtq_dz' else 'CCSD(T)',
            basis='cc-pVDZ', backend='mrcc' if task_id == 'ccsdtq_dz'
            else 'molpro', state_id=state_id, charge=0, multiplicity=1,
            geometry_sha256='4' * 64,
            source_sha256=hashlib.sha256(task_id.encode()).hexdigest(),
            source=f'{task_id}.out')

    monkeypatch.setattr(validation_module, 'task_component', component)
    monkeypatch.setattr(
        validation_module, 'audit_common_corrections_run',
        lambda run_dir: {'status': 'common_corrections_interface_complete'})
    result = audit_anl0_post_geometry_run('/synthetic/post')
    assert result['status'] == 'anl0_post_geometry_interface_complete'
    assert abs(result['higher_order_dz_hartree'] + 0.01) < 1e-12
    assert result['other_task_statuses']['ccsdtq_tz'] == 'failed'
    assert result['other_task_statuses']['ccsdtqp_dz'] == 'failed'


def test_higher_order_audit_returns_cross_program_correction(monkeypatch):
    spec = higher_order_validation_spec(_molecule())
    state = {'tasks': {task['id']: {'status': 'complete'}
                       for task in spec['tasks']}}
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (Path(run_dir), spec, state))
    values = {
        'ccsdt_tz': -10.000, 'ccsdt_dz': -9.900,
        'ccsdtq_tz': -10.020, 'ccsdtq_dz': -9.915,
        'ccsdtqp_dz': -9.920,
    }

    def component(run_dir, task_id, *, key, state_id):
        task = next(item for item in spec['tasks'] if item['id'] == task_id)
        parser = task['result_parser']
        return ComponentResult(
            key=key, value_hartree=values[task_id], quantity='electronic',
            method=parser['method'], basis=parser['basis'],
            backend=task['backend'], state_id=state_id, charge=0,
            multiplicity=1, geometry_sha256='3' * 64,
            source_sha256=hashlib.sha256(task_id.encode()).hexdigest(),
            source=f'{task_id}.out',
            settings={name: parser[name] for name in
                      ('reference', 'correlation', 'driver', 'program')
                      if name in parser})

    monkeypatch.setattr(validation_module, 'task_component', component)
    result = audit_higher_order_run('/synthetic/higher')
    assert result['status'] == 'higher_order_interface_complete'
    assert abs(result['correction_hartree'] + 0.025) < 1e-12
    assert result['components']['ccsdtq_tz']['backend'] == 'mrcc'
    assert result['components']['ccsdtqp_dz']['backend'] == 'mrcc'


def test_targeted_higher_order_audit_returns_only_available_difference(
        monkeypatch):
    spec = higher_order_validation_spec(
        _molecule(), task_ids=('ccsdtq_dz', 'ccsdtqp_dz'))
    state = {'tasks': {task['id']: {'status': 'complete'}
                       for task in spec['tasks']}}
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (Path(run_dir), spec, state))
    values = {'ccsdtq_dz': -40.3875184,
              'ccsdtqp_dz': -40.387527542311}

    def component(run_dir, task_id, *, key, state_id):
        task = next(item for item in spec['tasks'] if item['id'] == task_id)
        parser = task['result_parser']
        return ComponentResult(
            key=key, value_hartree=values[task_id], quantity='electronic',
            method=parser['method'], basis=parser['basis'],
            backend=task['backend'], state_id=state_id, charge=0,
            multiplicity=1, geometry_sha256='3' * 64,
            source_sha256=hashlib.sha256(task_id.encode()).hexdigest(),
            source=f'{task_id}.out',
            settings={name: parser[name] for name in
                      ('reference', 'correlation', 'driver', 'program')
                      if name in parser})

    monkeypatch.setattr(validation_module, 'task_component', component)
    result = audit_higher_order_run('/synthetic/targeted')
    assert result['status'] == 'targeted_higher_order_interface_complete'
    assert set(result['corrections_hartree']) == {'delta_p_dz'}
    assert abs(result['corrections_hartree']['delta_p_dz']
               + 0.000009142311) < 1e-12
    assert 'correction_hartree' not in result


def test_completed_run_geometry_export_is_hash_verified():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        spec = interface_validation_spec(_molecule())
        task = next(item for item in spec['tasks']
                    if item['id'] == 'l3_geometry')
        task['geometry_from'] = 'initial'
        task['resources'].update(cores=2, memory_mb=32000)
        spec['tasks'] = [task]
        run_dir = prepare(spec, root / 'run')
        directory = run_dir / 'tasks' / task['id']
        atoms = read(directory / 'geometry.xyz')
        write(directory / 'final.xyz', atoms)
        (directory / 'optimization.log').write_text('accepted\n')
        generated = [f"{task['id']}_step_0001{suffix}"
                     for suffix in ('.inp', '.out', '.log', '.xyz')]
        for name in generated:
            (directory / name).write_text('native\n')
        names = {'geometry.xyz', 'task.json', 'final.xyz',
                 'optimization.log', *generated}
        artifacts = {name: hashlib.sha256(
            (directory / name).read_bytes()).hexdigest() for name in names}
        final_hash = _geometry_hash(atoms)
        (directory / 'execution.json').write_text(json.dumps({
            'schema': 1, 'task_id': task['id'], 'status': 'executed',
            'geometry_sha256': final_hash, 'artifacts': artifacts,
            'details': {'final_geometry_sha256': final_hash,
                        'generated_files': generated}}))
        state_file = run_dir / 'state.json'
        state = json.loads(state_file.read_text())
        state['tasks'][task['id']].update(
            status='complete', final_geometry_sha256=final_hash)
        state_file.write_text(json.dumps(state))
        molecule = molecule_from_completed_run(run_dir)
        assert molecule['source']['geometry_sha256'] == final_hash
        assert molecule['multiplicity'] == 1
        (directory / 'final.xyz').write_text(
            (directory / 'final.xyz').read_text() + '\n')
        try:
            molecule_from_completed_run(run_dir)
        except RuntimeError as error:
            assert 'artifact' in str(error)
        else:
            raise AssertionError('Changed optimized geometry was accepted.')


def test_profiled_anl0_f12_assembler_uses_all_required_components(monkeypatch):
    interface = Path('/synthetic/interface').resolve()
    base = Path('/synthetic/current-base').resolve()
    interface_spec = {
        'name': 'anl1-f12-non-mrcc-interface-validation',
        'molecule': _molecule(), 'tasks': []}
    base_spec = {
        'name': 'anl-current-base-validation',
        'molecule': {**_molecule(), 'source': {
            'run_dir': str(interface), 'task_id': 'l2_geometry',
            'geometry_sha256': '2' * 64}}, 'tasks': []}
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: ((base, base_spec, {'tasks': {}})
                         if Path(run_dir).resolve() == base else
                         (interface, interface_spec, {'tasks': {}})))
    l2_hash, l3_hash = '2' * 64, '3' * 64
    monkeypatch.setattr(
        validation_module, 'molecule_from_completed_run',
        lambda run_dir, task: {'source': {'geometry_sha256':
                                           l2_hash if task == 'l2_geometry'
                                           else l3_hash}})
    monkeypatch.setattr(
        validation_module, '_same_imported_l3_geometry',
        lambda *args: None)
    equation = recipe('ANL0-F12', vpt2_method='B2PLYP-D3BJ')
    values = {item.key: (index + 1) / 1000
              for index, item in enumerate(equation.requirements)}
    components = {}
    for item in equation.requirements:
        if item.key == 'spin_orbit':
            continue
        geometry = l2_hash if item.geometry_role == 'l2' else l3_hash
        components[item.key] = ComponentResult(
            key=item.key, value_hartree=values[item.key],
            quantity=item.quantity, method=item.method, basis=item.basis,
            backend=item.backends[0], state_id='ethane-singlet', charge=0,
            multiplicity=1, geometry_sha256=geometry,
            source_sha256=hashlib.sha256(item.key.encode()).hexdigest(),
            source=f'{item.key}.out', settings=dict(item.settings))
    monkeypatch.setattr(
        validation_module, 'cbs_task_component',
        lambda *args, **kwargs: components['reference_cbs'])
    monkeypatch.setattr(
        validation_module, 'task_component',
        lambda run_dir, task_id, *, key, state_id: components[key])
    monkeypatch.setattr(
        validation_module, 'core_valence_task_component',
        lambda *args, **kwargs: components['core_valence_cbs'])
    monkeypatch.setattr(
        validation_module, 'scalar_relativistic_task_component',
        lambda *args, **kwargs: components['scalar_relativistic'])
    result = assemble_profiled_anl0_f12(
        interface, 'higher', 'corrections', state_id='ethane-singlet',
        spin_orbit_hartree=0.002, spin_orbit_source='test datum',
        base_run=base)
    assert result['status'] == 'complete'
    assert result['recipe'].startswith('profiled:ANL0-F12:')
    assert set(result['components']) == {
        item.key for item in equation.requirements}
    expected_electronic = sum(
        (0.002 if term.component == 'spin_orbit'
         else values[term.component]) * term.coefficient
        for term in equation.electronic_terms)
    expected_zpe = sum(values[term.component] * term.coefficient
                       for term in equation.zero_point_terms)
    assert result['electronic_hartree'] == expected_electronic
    assert result['zero_point_hartree'] == expected_zpe
    assert result['zero_k_hartree'] == expected_electronic + expected_zpe


def test_profiled_anl0_f12_assembler_accepts_one_composite_run(monkeypatch):
    run = Path('/synthetic/composite').resolve()
    spec = {'name': 'anl-composite-validation',
            'molecule': _molecule(), 'tasks': []}
    monkeypatch.setattr(
        validation_module, '_load',
        lambda run_dir: (run, spec, {'tasks': {}}))
    l2_hash, l3_hash = '2' * 64, '3' * 64
    monkeypatch.setattr(
        validation_module, 'molecule_from_completed_run',
        lambda run_dir, task: {'source': {'geometry_sha256':
                                           l2_hash if task == 'l2_geometry'
                                           else l3_hash}})

    def unexpected_import(*args):
        raise AssertionError('A single composite run must not be treated as '
                             'an imported child graph.')

    monkeypatch.setattr(
        validation_module, '_same_imported_geometry', unexpected_import)
    monkeypatch.setattr(
        validation_module, '_same_imported_l3_geometry', unexpected_import)
    equation = recipe('ANL0-F12', vpt2_method='B2PLYP-D3BJ')
    components = {}
    for index, item in enumerate(equation.requirements, 1):
        if item.key == 'spin_orbit':
            continue
        components[item.key] = ComponentResult(
            key=item.key, value_hartree=index / 1000,
            quantity=item.quantity, method=item.method, basis=item.basis,
            backend=item.backends[0], state_id='methane-singlet', charge=0,
            multiplicity=1,
            geometry_sha256=(l2_hash if item.geometry_role == 'l2'
                             else l3_hash),
            source_sha256=hashlib.sha256(item.key.encode()).hexdigest(),
            source=f'{item.key}.out', settings=dict(item.settings))
    monkeypatch.setattr(
        validation_module, 'cbs_task_component',
        lambda *args, **kwargs: components['reference_cbs'])
    monkeypatch.setattr(
        validation_module, 'task_component',
        lambda run_dir, task_id, *, key, state_id: components[key])
    monkeypatch.setattr(
        validation_module, 'core_valence_task_component',
        lambda *args, **kwargs: components['core_valence_cbs'])
    monkeypatch.setattr(
        validation_module, 'scalar_relativistic_task_component',
        lambda *args, **kwargs: components['scalar_relativistic'])
    result = assemble_profiled_anl0_f12(
        run, run, run, state_id='methane-singlet',
        spin_orbit_hartree=0., spin_orbit_source='closed-shell test',
        spin_orbit_backend='known_zero')
    assert result['status'] == 'complete'
    assert result['geometry_sha256'] == {'l2': l2_hash, 'l3': l3_hash}


def test_quality_review_is_bound_to_exact_flagged_native_output(tmp_path):
    component = ComponentResult(
        key='vpt2_correction', value_hartree=-0.001,
        quantity='correction', method='B2PLYP-D3BJ', basis='cc-pVTZ',
        backend='gaussian', state_id='ethane-singlet', charge=0,
        multiplicity=1, geometry_sha256='2' * 64,
        source_sha256='a' * 64, source='gaussian_vpt2.log',
        settings={'dispersion': 'GD3BJ'}, review_required=True)
    review = tmp_path / 'vpt2_review.json'
    review.write_text(json.dumps({
        'schema': 1, 'task_id': 'gaussian_vpt2',
        'native_output_sha256': 'a' * 64, 'decision': 'accept',
        'reviewer': 'test reviewer',
        'rationale': 'Reviewed the named Gaussian warnings and mode table.',
    }))
    accepted = validation_module._apply_quality_review(
        component, 'gaussian_vpt2', review)
    assert accepted.review_required is False
    assert accepted.settings['quality_review']['native_output_sha256'] == \
        component.source_sha256
    assert accepted.source_sha256 != component.source_sha256
    review.write_text(review.read_text().replace('a' * 64, 'b' * 64))
    try:
        validation_module._apply_quality_review(
            component, 'gaussian_vpt2', review)
    except ValueError as error:
        assert 'exact native output' in str(error)
    else:
        raise AssertionError('Mismatched VPT2 review hash was accepted.')


def test_database_export_requires_one_complete_accepted_l2_record():
    with TemporaryDirectory() as temporary:
        database = Path(temporary) / 'kinbot.db'
        db = connect(database)
        atoms = Atoms('CH4', positions=_molecule()['positions'])
        db.write(atoms, name='methane_well_high', data={
            'status': 'normal', 'energy': -10., 'zpe': .1,
            'frequencies': [100., 200.]})
        molecule = molecule_from_database(
            database, 'methane_well_high', charge=0, multiplicity=1)
        assert molecule['symbols'] == ['C', 'H', 'H', 'H', 'H']
        assert molecule['charge'] == 0
        db.write(atoms, name='bad_well_high', data={'status': 'error'})
        try:
            molecule_from_database(database, 'bad_well_high', charge=0,
                                   multiplicity=1)
        except ValueError as error:
            assert 'not normal' in str(error)
        else:
            raise AssertionError('An incomplete L2 result was exported.')


def test_ethane_external_example_finds_only_requested_scission():
    root = Path(__file__).resolve().parents[1]
    parameters = Parameters(
        root / 'examples/anl/ethane_profiled_hpc/ethane.json',
        show_warnings=False).par
    species = StationaryPoint(
        'well0', parameters['charge'], parameters['mult'],
        smiles=parameters['smiles'], structure=parameters['structure'])
    species.characterize()
    assert str(species.chemid) == '301020900180000000001'
    species.name = str(species.chemid)
    ReactionFinder(species, parameters, None).find_reactions()
    assert species.reac_type == ['hom_sci']
    assert species.reac_inst == [[0, 1]]


def test_prepare_from_database_cli_stages_general_graph(monkeypatch):
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        database = root / 'kinbot.db'
        connect(database).write(
            Atoms('CH4', positions=_molecule()['positions']),
            name='methane_well_high', data={
                'status': 'normal', 'energy': -10., 'zpe': .1,
                'frequencies': [100., 200.]})
        monkeypatch.setattr(site, '_partitions', lambda: [{
            'name': 'test-cpu', 'default': True, 'cores': 32,
            'memory_mb': 128000, 'seconds': 7 * 24 * 3600}])
        run_dir = root / 'run'
        assert validation_main([
            'prepare-from-db', str(database), 'methane_well_high',
            str(run_dir), '--max-nodes', '2', '--partition', 'test-cpu']) == 0
        workflow = json.loads((run_dir / 'workflow.json').read_text())
        assert workflow['limits']['max_nodes'] == 2
        assert all(task['resources']['partition'] == 'test-cpu'
                   for task in workflow['tasks'])
        state = json.loads((run_dir / 'state.json').read_text())
        assert set(state['tasks']) == {'l2_geometry'}


def test_smiles_initial_geometry_is_deterministic_and_state_checked():
    first = molecule_from_smiles('C', charge=0, multiplicity=1)
    second = molecule_from_smiles('C', charge=0, multiplicity=1)
    assert first == second
    assert first['symbols'].count('C') == 1
    assert first['symbols'].count('H') == 4
    assert first['smiles'] == 'C'
    hydrogen = molecule_from_smiles(
        '[H][H]', charge=0, multiplicity=1)
    assert hydrogen['symbols'] == ['H', 'H']
    assert hydrogen['positions'][0] != hydrogen['positions'][1]
    hydrogen_tasks = {
        task['id']: task for task in interface_validation_spec(hydrogen)['tasks']}
    assert hydrogen_tasks['l2_geometry']['optimizer']['sella_kwargs'] == {
        'internal': False}
    assert hydrogen_tasks['l3_geometry']['optimizer']['sella_kwargs'] == {
        'internal': False}
    with pytest.raises(ValueError, match='formal charge'):
        molecule_from_smiles('[NH4+]', charge=0, multiplicity=1)
    with pytest.raises(ValueError, match='radical electrons'):
        molecule_from_smiles('[CH3]', charge=0, multiplicity=1)


def test_prepare_from_smiles_cli_stages_general_graph(monkeypatch):
    captured = {}

    def fake_prepare(spec, run_dir):
        captured['spec'] = spec
        captured['run_dir'] = Path(run_dir)
        return Path(run_dir)

    monkeypatch.setattr(validation_module, 'prepare', fake_prepare)
    assert validation_main([
        'prepare-from-smiles', 'C', '/synthetic/methane',
        '--max-nodes', '2', '--partition', 'day-long-cpu']) == 0
    assert captured['run_dir'] == Path('/synthetic/methane')
    assert captured['spec']['molecule']['smiles'] == 'C'
    assert captured['spec']['limits']['max_nodes'] == 2
    assert all(task['resources']['partition'] == 'day-long-cpu'
               for task in captured['spec']['tasks'])


def test_prepare_composite_from_smiles_cli_stages_one_graph(monkeypatch):
    captured = {}

    def fake_prepare(spec, run_dir):
        captured['spec'] = spec
        captured['run_dir'] = Path(run_dir)
        return Path(run_dir)

    monkeypatch.setattr(validation_module, 'prepare', fake_prepare)
    assert validation_main([
        'prepare-composite-from-smiles', 'C', '/synthetic/methane-full',
        '--max-nodes', '2', '--partition', 'day-long-cpu',
        '--anl0-only']) == 0
    spec = captured['spec']
    assert spec['name'] == 'anl-composite-validation'
    assert spec['intent']['requested_ladder_head'] == 'ANL0-F12'
    assert spec['limits']['max_nodes'] == 2
    assert captured['run_dir'] == Path('/synthetic/methane-full')
    assert all(task['resources']['partition'] == 'day-long-cpu'
               for task in spec['tasks'])


def test_prepare_post_geometry_cli_builds_one_shared_graph(monkeypatch):
    molecule = deepcopy(_molecule())
    molecule['source'] = {
        'run_dir': '/synthetic/current', 'task_id': 'l3_geometry',
        'geometry_sha256': '5' * 64, 'artifact_sha256': '6' * 64}
    captured = {}
    monkeypatch.setattr(
        validation_module, 'molecule_from_completed_run',
        lambda source_run, geometry_task: molecule)

    def fake_prepare(spec, run_dir):
        captured['spec'] = spec
        captured['run_dir'] = Path(run_dir)
        return Path(run_dir)

    monkeypatch.setattr(validation_module, 'prepare', fake_prepare)
    assert validation_main([
        'prepare-post-geometry-from-run', '/synthetic/current',
        '/synthetic/post', '--max-nodes', '3',
        '--partition', 'day-long-cpu']) == 0
    assert captured['spec']['name'] == 'anl-post-geometry-validation'
    assert captured['spec']['limits']['max_nodes'] == 3
    assert len(captured['spec']['tasks']) == 11
    assert captured['run_dir'] == Path('/synthetic/post')

    assert validation_main([
        'prepare-post-geometry-from-run', '/synthetic/current',
        '/synthetic/post-anl0', '--max-nodes', '3', '--anl0-only']) == 0
    assert {task['id'] for task in captured['spec']['tasks']} == {
        'ccsdt_dz', 'ccsdtq_dz',
        'cv_ae_tz', 'cv_ae_qz', 'cv_fc_tz', 'cv_fc_qz',
        'rel_dkh', 'rel_nonrel'}


def test_prepare_higher_order_cli_can_stage_only_the_mrcc_pair(monkeypatch):
    molecule = deepcopy(_molecule())
    molecule['source'] = {
        'run_dir': '/synthetic/qz', 'task_id': 'l3_geometry',
        'geometry_sha256': '5' * 64, 'artifact_sha256': '6' * 64}
    captured = {}
    monkeypatch.setattr(
        validation_module, 'molecule_from_completed_run',
        lambda source_run, geometry_task: molecule)

    def fake_prepare(spec, run_dir):
        captured['spec'] = spec
        return Path(run_dir)

    monkeypatch.setattr(validation_module, 'prepare', fake_prepare)
    assert validation_main([
        'prepare-higher-order-from-run', '/synthetic/qz',
        '/synthetic/targeted', '--task', 'ccsdtq_dz',
        '--task', 'ccsdtqp_dz']) == 0
    assert [task['id'] for task in captured['spec']['tasks']] == [
        'ccsdtq_dz', 'ccsdtqp_dz']


def test_kinbot_gate_requires_accepted_reaction_hir_and_rotdpy():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        reaction = 'parent_hom_sci_1_2'
        (root / 'kinbot_monitor.out').write_text(
            f'-1\t0\t{reaction}\tch3 ch3\n')
        (root / 'kinbot.log').write_text('Reaction generation done!\n')
        database = connect(root / 'kinbot.db')
        for point in range(4):
            database.write(
                Atoms('H', positions=[[0., 0., float(point)]]),
                name=f'hir/parent_hir_0_0{point}',
                data={'status': 'normal'})
        (root / 'vrctst').mkdir()
        correction = {
            'dist': [30.], 'e_samp': [0.], 'e_high': [0.],
            'scan_ref': [[0, 0]], 'ra': [[0], [0]],
            'e_inf_samp': -1., 'e_inf_high': -1.,
            'frags_atom': [['C'], ['C']],
            'frags_geom': [[[0., 0., 0.]], [[1., 0., 0.]]],
            'frags_mult': [2, 2],
        }
        (root / 'vrctst' / f'corr_{reaction}.json').write_text(
            json.dumps(correction))
        (root / 'rotdPy').mkdir()
        rotdpy_input = root / 'rotdPy' / f'{reaction}.py'
        rotdpy_input.write_text('# rotdPy input\n')

        input_only = audit_kinbot_run(
            root, reaction, parent='parent', hir_points=4,
            require_rotdpy=True)
        assert input_only['status'] == 'kinbot_reaction_complete'
        assert input_only['rotdpy_input'] == str(rotdpy_input.resolve())
        assert input_only['rotdpy_execution'] is None
        assert input_only['rotdpy_surfaces'] == 0

        result_root = root / 'rotdPy' / f'kb_{reaction}'
        surface = result_root / 'output' / 'surface_0.dat'
        states = result_root / 'Ne_0.out'
        surface.parent.mkdir(parents=True)
        surface.write_text('surface flux\n')
        states.write_text('0.0 1.0\n')
        result_files = [str(states.relative_to(root / 'rotdPy')),
                        str(surface.relative_to(root / 'rotdPy'))]
        (root / 'rotdPy' / f'{reaction}.rotdpy.json').write_text(
            json.dumps({'schema': 2, 'status': 'complete',
                        'reaction': reaction, 'surface_count': 1,
                        'result_files': result_files,
                        'result_sha256': {
                            name: hashlib.sha256(
                                (root / 'rotdPy' / name).read_bytes()).hexdigest()
                            for name in result_files}}))
        (root / 'rotdPy' / f'{reaction}.execution.json').write_text(
            json.dumps({'schema': 1, 'status': 'complete', 'returncode': 0,
                        'input_sha256': hashlib.sha256(
                            rotdpy_input.read_bytes()).hexdigest()}))

        result = audit_kinbot_run(
            root, reaction, parent='parent', hir_points=4,
            require_rotdpy_execution=True)
        assert result['status'] == 'kinbot_reaction_complete'
        assert result['products'] == ['ch3', 'ch3']
        assert result['normal_hir_points'] == 4
        assert result['rotdpy_surfaces'] == 1
