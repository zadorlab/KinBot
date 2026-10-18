"""Recovery preserves live ownership and rejects wrong-geometry reuse."""
from copy import deepcopy
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from ase import Atoms

from kinbot.anl import dispatch, salvage, automatic
from kinbot.anl.validation import automatic_composite_spec
from kinbot import run_salvage
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def small_spec():
    spec = dispatch_spec()
    first, geometry = deepcopy(spec['tasks'][:2])
    first['id'] = 'accepted_geometry'
    geometry['geometry_from'] = 'accepted_geometry'
    geometry.update(failure_policy='skip_component', component_families=['geometry'])
    sp = molpro_task('energy', 'rhf\n')
    sp['geometry_fallback'] = 'initial'
    vpt2 = molpro_task('l2_property', 'rhf\n', geometry_from='initial')
    spec['tasks'] = [first, geometry, sp, vpt2]
    spec['intent'] = {}
    return spec


def test_wrong_geometry_invalidates_descendants_but_preserves_l2_property():
    spec = small_spec()
    state = {'tasks': {
        'accepted_geometry': {'status': 'complete'},
        'l3_geometry': {'status': 'unavailable'},
        'energy': {'status': 'complete', 'geometry_from': 'initial'},
        'l2_property': {'status': 'complete', 'geometry_from': 'initial'}}}
    assert salvage.geometry_invalidations(spec, state) == {'l3_geometry', 'energy'}
    state['tasks']['l3_geometry'] = {'status': 'complete', 'final_geometry_sha256': 'a'}
    state['tasks']['energy'] = {'status': 'complete', 'geometry_from': 'l3_geometry', 'geometry_sha256': 'b'}
    assert salvage.geometry_invalidations(spec, state) == {'energy'}


def test_repair_archives_fallback_energy_and_blocks_it_until_new_geometry(tmp_path):
    graph = dispatch.prepare(small_spec(), tmp_path/'graph')
    _mock_execution(graph, 'accepted_geometry', geometry='final.xyz')
    _mock_execution(graph, 'l2_property')
    dispatch.advance(graph)
    path = graph/'state.json'
    state = json.loads(path.read_text())
    state['tasks']['l3_geometry']['status'] = 'unavailable'
    atomic_json(path, state)
    dispatch.advance(graph)
    _mock_execution(graph, 'energy')
    dispatch.advance(graph)
    original = (graph/'tasks/l2_property/task.json').read_bytes()
    report = salvage.repair_graph(graph, tmp_path/'archive')
    state = dispatch._load(graph)[2]
    assert state['tasks']['l3_geometry']['status'] == 'staged'
    assert state['tasks'].get('energy', {}).get('status', 'waiting') == 'waiting'
    assert not (graph/'tasks/energy').exists()
    assert (Path(report['archive'])/'energy/execution.json').exists()
    assert (graph/'tasks/l2_property/task.json').read_bytes() == original
    assert state['tasks']['l2_property']['status'] == 'complete'
    _mock_execution(graph, 'l3_geometry', geometry='final.xyz')
    state = dispatch.advance(graph)
    assert state['tasks']['energy']['geometry_from'] == 'l3_geometry'
    assert state['tasks']['energy']['geometry_sha256'] == state['tasks']['l3_geometry']['final_geometry_sha256']
    assert salvage.repair_graph(graph, tmp_path/'archive') == report


def test_live_l3_is_adopted_without_changing_its_input_or_submitting(tmp_path):
    graph = dispatch.prepare(small_spec(), tmp_path/'graph')
    _mock_execution(graph, 'accepted_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    state = dispatch._load(graph)[2]
    state['tasks']['l3_geometry'].update(status='submitted', job_id='101')
    atomic_json(graph/'state.json', state)
    hashes = {p.name: file_sha256(p) for p in (graph/'tasks/l3_geometry').iterdir() if p.is_file()}
    with patch.object(dispatch, '_job_active', return_value=True), patch.object(dispatch, 'subprocess') as commands:
        report = salvage.repair_graph(graph, tmp_path/'archive')
    commands.run.assert_not_called()
    assert report['adopted_jobs'] == {'l3_geometry': '101'}
    assert hashes == {p.name: file_sha256(p) for p in (graph/'tasks/l3_geometry').iterdir() if p.is_file()}


def test_incompatible_live_task_is_never_moved_by_graph_repair(tmp_path):
    graph = dispatch.prepare(small_spec(), tmp_path/'graph')
    _mock_execution(graph, 'accepted_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    state = dispatch._load(graph)[2]
    state['tasks']['l3_geometry']['status'] = 'unavailable'
    atomic_json(graph/'state.json', state)
    dispatch.advance(graph)
    state = dispatch._load(graph)[2]
    state['tasks']['energy'].update(status='submitted', job_id='102')
    atomic_json(graph/'state.json', state)
    with patch.object(dispatch, '_job_active', return_value=True):
        with pytest.raises(RuntimeError, match='incompatible live job'):
            salvage.repair_graph(graph, tmp_path/'archive')
    assert (graph/'tasks/energy').exists()
    assert not (tmp_path/'archive').exists()


def test_repair_transaction_rolls_back_on_staging_failure(tmp_path, monkeypatch):
    graph = dispatch.prepare(small_spec(), tmp_path/'graph')
    before = {name: (graph/name).read_bytes() for name in ('state.json', 'workflow.json')}
    stage = dispatch._stage_task
    monkeypatch.setattr(dispatch, '_stage_task', lambda *a, **kw: (_ for _ in ()).throw(RuntimeError('test staging failure')))
    with pytest.raises(RuntimeError, match='test staging failure'):
        salvage.repair_graph(graph, tmp_path/'archive')
    for name, data in before.items():
        assert (graph/name).read_bytes() == data
    monkeypatch.setattr(dispatch, '_stage_task', stage)
    dispatch.advance(graph)


def test_l3_retry_is_bounded_and_does_not_change_settings(tmp_path):
    spec = small_spec()
    spec['intent'] = {'require_l3_geometry': True, 'l3_max_retries': 1}
    graph = dispatch.prepare(spec, tmp_path/'graph')
    _mock_execution(graph, 'accepted_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    path = graph/'state.json'
    state = dispatch._load(graph)[2]
    state['tasks']['l3_geometry'].update(status='failed', failure_phase='execution')
    atomic_json(path, state)
    task = (graph/'tasks/l3_geometry/task.json').read_bytes()
    result = automatic._retry_l3_execution_failures(graph, spec, state)
    assert result['tasks']['l3_geometry']['attempt'] == 2
    assert (graph/'tasks/l3_geometry/task.json').read_bytes() == task
    assert (graph/'attempts/l3_geometry/1/task.json').exists()
    result['tasks']['l3_geometry']['status'] = 'failed'
    result['tasks']['l3_geometry']['failure_phase'] = 'execution'
    assert automatic._retry_l3_execution_failures(graph, spec, result) == result
    assert not (graph/'attempts/l3_geometry/2').exists()


def test_anl_default_requires_primary_l3_but_vpt2_keeps_l2():
    mol = small_spec()['molecule']
    mol['source'] = {'geometry_sha256': dispatch._geometry_hash(Atoms(mol['symbols'], positions=mol['positions'])),
                     'record_sha256': '2'*64, 'energy_hartree': -40.}
    spec = automatic_composite_spec(mol, highest_tier='ANL1-F12')
    assert spec['intent']['require_l3_geometry'] is True
    assert spec['intent']['l3_max_retries'] == 1
    assert all('geometry_fallback' not in t for t in spec['tasks'])
    tasks = {t['id']: t for t in spec['tasks']}
    assert tasks['gaussian_vpt2']['geometry_from'] == 'initial'
    assert tasks['anl1_6z']['geometry_from'] == 'l3_geometry'
    assert 'failure_policy' not in tasks['l3_geometry']


def test_exchange_retains_open_directory_inode_and_old_path(tmp_path, monkeypatch):
    old, new = tmp_path/'old', tmp_path/'new'
    old.mkdir()
    (old/'result').write_text('running')
    inode = old.stat().st_ino
    # Portable stand-in for the Linux atomic primitive; test the path contract.
    def exchange(a, b):
        temporary = tmp_path/'exchange'
        a.rename(temporary)
        b.rename(a)
        temporary.rename(b)
    monkeypatch.setattr(run_salvage, '_exchange', exchange)
    assert run_salvage.rehome(old, new) == new
    assert old.is_symlink() and old.resolve() == new
    assert new.stat().st_ino == inode
    (old/'result').write_text('complete')
    assert (new/'result').read_text() == 'complete'
    assert run_salvage.rehome(old, new) == new


def test_failed_exchange_leaves_original_tree(tmp_path, monkeypatch):
    old, new = tmp_path/'old', tmp_path/'new'
    old.mkdir()
    def fail(*args):
        raise OSError('unsupported filesystem')
    monkeypatch.setattr(run_salvage, '_exchange', fail)
    with pytest.raises(OSError):
        run_salvage.rehome(old, new)
    assert old.is_dir() and not old.is_symlink()
    assert not new.is_symlink()


def test_cancellation_selects_only_incompatible_owned_job(tmp_path, monkeypatch):
    root = tmp_path/'run'
    graph = root/'case/123/anl_automatic/states/a/composite'
    graph.parent.mkdir(parents=True)
    dispatch.prepare(small_spec(), graph)
    _mock_execution(graph, 'accepted_geometry', geometry='final.xyz')
    state = dispatch.advance(graph)
    state['tasks']['l3_geometry']['status'] = 'unavailable'
    atomic_json(graph/'state.json', state)
    state = dispatch.advance(graph)
    state['tasks']['energy'].update(status='submitted', job_id='102')
    atomic_json(graph/'state.json', state)
    jobs = [{'job_id': '102', 'directory': str(graph/'tasks/energy')},
            {'job_id': '101', 'directory': str(root/'other')}]
    cancelled = []
    monkeypatch.setattr(run_salvage, '_inventory', lambda root: {'jobs': jobs, 'processes': []})
    def cancel(root, selected):
        cancelled.extend(j['job_id'] for j in selected)
        jobs[:] = [j for j in jobs if j not in selected]
    monkeypatch.setattr('kinbot.run_cleanup._cancel', cancel)
    run_salvage._cancel_incompatible(root, {})
    assert cancelled == ['102']
    assert jobs[0]['job_id'] == '101'


def test_shared_filesystem_without_exchange_keeps_old_tree_and_new_view(tmp_path, monkeypatch):
    import errno
    old, new = tmp_path/'old', tmp_path/'new'
    old.mkdir()
    inode = old.stat().st_ino
    def unsupported(*args):
        raise OSError(errno.EOPNOTSUPP, 'shared filesystem')
    monkeypatch.setattr(run_salvage, '_exchange', unsupported)
    assert run_salvage.rehome(old, new) == old
    assert old.stat().st_ino == inode and not old.is_symlink()
    assert new.is_symlink() and new.resolve() == old
    (old/'native.out').write_text('still writing')
    assert (new/'native.out').read_text() == 'still writing'
    assert run_salvage.rehome(old, new) == old
