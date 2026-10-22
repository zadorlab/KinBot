"""Concurrent leaves preserve orbital chains, graph ownership and refreshes."""
from copy import deepcopy
import fcntl
import json
from pathlib import Path
import threading
import time
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot import pes_phase, vrc_parallel, vrc_multistate, vrc_graph
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from kinbot.vrc_tst_scan import VTS
from tests.test_vrc_graph import graph, spec
from tests.test_vrc_multistate import model
from tests.test_rotdpy_profiles import sulfur_correction


def surfaces(tmp_path):
    par, _, _ = sulfur_correction(tmp_path)
    par['vrc_tst_parallel_surfaces'] = 4
    vts = VTS(SimpleNamespace(), par, None)
    channel = vrc_multistate.validate_surface_model(model())['channels'][0]
    channel['bottlenecks'] = channel['bottlenecks'][:1]
    channel['combination'] = 'single'
    choices = {}
    for name in ('cs', 'ss'):
        vts.scan_reac[name] = SimpleNamespace(species=SimpleNamespace(mult=2),
            products=[SimpleNamespace(mult=2, geom=np.zeros((2, 3))),
                      SimpleNamespace(mult=3, geom=np.ones((1, 3)))])
        choices[name] = (deepcopy(channel), {'source': 'review.json'})
    return vts, choices


def test_bonds_and_roots_overlap_but_each_scan_precedes_its_corrections(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    vts, choices = surfaces(tmp_path)
    rendezvous = threading.Barrier(4, timeout=5)
    seen, done = {}, set()
    mutex = threading.Lock()
    def scan(child, names):
        name = names[0]
        leaf = child.scan_reac[name]
        with mutex:
            seen[name] = (threading.get_ident(), leaf.species.mult, deepcopy(child.par))
        # Every leaf reaches its scan before ANY leaf can finish.
        rendezvous.wait()
        leaf.products[0].geom[:] = 7
        assert np.all(leaf.products[1].geom == 1)
        with mutex:
            done.add(name)
    def energies(child, names):
        assert names[0] in done
        Path(f'vrctst/corr_{names[0]}.json').write_text('{}')
    monkeypatch.setattr(VTS, 'do_scan', scan)
    monkeypatch.setattr(VTS, 'energies', energies)
    paths = vrc_multistate.calculate_channels(vts, choices, limit=4)
    assert set(paths) == {'cs', 'ss'}
    assert len({v[0] for v in seen.values()}) == 4
    assert set(seen) == {f'{name}__inner_m2_r{root}' for name in choices for root in (1, 2)}
    for name, reaction in vts.scan_reac.items():
        assert reaction.species.mult == 2
        assert np.all(reaction.products[0].geom == 0)
        payload = json.loads(paths[name].read_text())
        assert [s['root'] for s in payload['bottlenecks'][0]['states']] == [1, 2]
        for leaf in payload['bottlenecks'][0]['states']:
            assert file_sha256(leaf['correction']) == leaf['correction_sha256']
    status = json.loads(Path('vrctst/parallel_surfaces.json').read_text())
    assert {s['status'] for s in status['surfaces'].values()} == {'complete'}


def test_failed_surface_does_not_stop_other_channel_and_plotting_stays_on_caller(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    vts, choices = surfaces(tmp_path)
    caller = threading.get_ident()
    plots = []
    def scan(child, names):
        if names[0].startswith('cs'):
            raise RuntimeError('Native state validation failed')
    def energies(child, names):
        child._vrc_plot = {'name': names[0]}
        Path(f'vrctst/corr_{names[0]}.json').write_text('{}')
    monkeypatch.setattr(VTS, 'do_scan', scan)
    monkeypatch.setattr(VTS, 'energies', energies)
    monkeypatch.setattr('kinbot.utils.create_matplotlib_graph',
                        lambda **kw: plots.append(threading.get_ident()))
    with pytest.raises(RuntimeError, match='Native state validation'):
        vrc_multistate.calculate_channels(vts, choices, limit=4)
    assert not Path('vrctst/corr_cs.json').exists()
    assert Path('vrctst/corr_ss.json').exists()
    assert plots == [caller, caller]
    status = json.loads(Path('vrctst/parallel_surfaces.json').read_text())['surfaces']
    assert sum(v['status'] == 'failed' for v in status.values()) == 2
    assert sum(v['status'] == 'complete' for v in status.values()) == 2


def test_parallel_default_requires_shared_admission(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    vts, choices = surfaces(tmp_path)
    seen = []
    monkeypatch.setattr(vrc_parallel, 'run', lambda jobs, **kw: seen.append(kw['limit']))
    monkeypatch.setattr('kinbot.run_policy.find', lambda: None)
    vrc_multistate.calculate_channels(vts, choices)
    monkeypatch.setattr('kinbot.run_policy.find', lambda: tmp_path/'policy.json')
    vrc_multistate.calculate_channels(vts, choices)
    vts.par['vrc_tst_parallel_surfaces'] = 1
    vrc_multistate.calculate_channels(vts, choices)
    assert seen == [1, 4, 1]


def test_epoch_refresh_joins_all_drivers_releases_locks_and_keeps_allocations(tmp_path, monkeypatch):
    paths = [graph(tmp_path/f'leaf{i}', spec(), 'submitted', str(100+i)) for i in range(2)]
    before = {p: (p/'state.json').read_bytes() for p in paths}
    epoch = [1]
    monkeypatch.setattr(pes_phase, '_EPOCH', 1)
    monkeypatch.setattr(pes_phase, '_read', lambda: {'epoch': epoch[0]})
    rendezvous = threading.Barrier(2, action=lambda: epoch.__setitem__(0, 2), timeout=5)
    def advance(path, **kwargs):
        rendezvous.wait()
        return dispatch._load(path)[2]
    monkeypatch.setattr(dispatch, 'advance', advance)
    results = []
    with pytest.raises(pes_phase.ExplorationReopened):
        vrc_parallel.run([(str(p), lambda p=p: dispatch.drive_run(p, interval=.001)) for p in paths],
            limit=2, completed=lambda *args: results.append(args))
    assert len(results) == 2
    for p in paths:
        assert (p/'state.json').read_bytes() == before[p]
        with (p/'drive.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not any(t.name.startswith('kinbot-vrc-') for t in threading.enumerate())


def test_thread_start_failure_quiesces_already_started_controller(monkeypatch):
    original = threading.Thread.start
    calls = []
    def start(thread):
        calls.append(thread)
        if len(calls) == 2:
            raise RuntimeError('Cannot start another thread')
        original(thread)
    monkeypatch.setattr(threading.Thread, 'start', start)
    def poll():
        for _ in range(1000):
            pes_phase.checkpoint()
            time.sleep(.001)
        pytest.fail('Controller did not quiesce')
    with pytest.raises(RuntimeError, match='Cannot start another thread'):
        vrc_parallel.run([('a', poll), ('b', poll)], limit=2,
                         completed=lambda *args: None)
    assert not any(t.name.startswith('kinbot-vrc-') for t in threading.enumerate())


@pytest.mark.parametrize('limit', [0, -1, True, 1.5])
def test_invalid_parallel_limit_is_rejected_before_starting(limit):
    with pytest.raises(ValueError, match='positive integer'):
        vrc_parallel.run([('a', lambda: pytest.fail('Should not run'))], limit=limit,
                         completed=lambda *args: None)


def test_concurrent_identical_graph_selection_creates_one_legacy_binding(tmp_path, monkeypatch):
    calls = []
    original = dispatch.prepare
    def prepare(*args):
        calls.append(args)
        return original(*args)
    monkeypatch.setattr(dispatch, 'prepare', prepare)
    results = []
    vrc_parallel.run([(str(i), lambda: vrc_graph.select_graph(tmp_path, 'pivot', spec())) for i in range(4)],
        limit=4, completed=lambda name, value, error: results.append((value, error)))
    assert len(calls) == 1
    assert all(error is None for _, error in results)
    assert len({value for value, _ in results}) == 1


def test_shared_pivot_runs_once_and_is_reused_by_two_surfaces(tmp_path, monkeypatch):
    from kinbot import vrc_orbitals
    monkeypatch.chdir(tmp_path)
    value = spec('pivot')
    value['tasks'][0]['id'] = 'pivot'
    value['intent'] = {'orbital_only': True}
    monkeypatch.setattr(vrc_orbitals, 'fixed_pivot_spec', lambda *args: deepcopy(value))
    monkeypatch.setattr(dispatch, 'preflight_until_scheduler_ready', lambda *args, **kw: None)
    monkeypatch.setattr(vrc_orbitals, 'read_orbital_cube', lambda *args: {})
    monkeypatch.setattr(vrc_orbitals, 'validate_cube_geometry', lambda *args: None)
    allocations = []
    def drive(path, ids):
        state = dispatch._load(path)[2]
        if state['tasks']['pivot']['status'] != 'complete':
            allocations.append('123')
            state['tasks']['pivot'].update(status='complete', job_id='123')
            (path/'tasks/pivot/homo.cube').write_text('verified native cube')
            atomic_json(path/'state.json', state)
        return 0, {'pivot': 'complete'}
    monkeypatch.setattr('kinbot.vrc_tst_scan._drive_vrc_dispatch', drive)
    fragment = SimpleNamespace(atom=['H', 'H'], geom=np.array([[0, 0, 0], [0, 0, .74]]))
    results = []
    vrc_parallel.run([(str(i), lambda i=i: vrc_orbitals.prepare_fixed_pivot(fragment, {},
        tmp_path/f'leaf{i}/fragment_0', {'leaf': i})) for i in range(2)], limit=2,
        completed=lambda name, value, error: results.append((value, error)))
    assert all(error is None for _, error in results), results
    assert allocations == ['123']
    assert len(list(Path('vrctst/molpro_pivots').glob('*/workflow.json'))) == 1
    assert all(Path(value).read_text() == 'verified native cube' for value, _ in results)


def test_partial_correction_batch_reattaches_existing_jobs_without_rewriting(tmp_path):
    from tests.anl_fixture import molpro_task
    original = spec('corrections')
    original['tasks'] = [molpro_task(n, 'rhf\n', geometry_from='initial') for n in ('a', 'b')]
    path = tmp_path/'dispatch'
    dispatch.prepare(original, path)
    state = dispatch._load(path)[2]
    state['tasks']['a'].update(status='complete', job_id='101')
    state['tasks']['b'].update(status='submitted', job_id='102')
    atomic_json(path/'state.json', state)
    before = {p: file_sha256(p) for p in path.rglob('*') if p.is_file()}
    pending = deepcopy(original)
    pending['tasks'] = pending['tasks'][1:]
    pending['tasks'][0]['resources']['walltime'] = '31-00:00:00'
    assert vrc_graph.select_correction_graph(tmp_path, pending) == path
    assert {p: file_sha256(p) for p in before} == before
    assert len(list(tmp_path.glob('dispatch*/workflow.json'))) == 1
    pending['tasks'][0]['input_template'] += 'different scientific input\n'
    assert vrc_graph.select_correction_graph(tmp_path, pending) != path


@pytest.mark.parametrize('change', ['command', 'charge', 'initial_geometry'])
def test_correction_reuse_does_not_ignore_scientific_context(tmp_path, change):
    original = spec('corrections')
    original['tasks'][0]['input_template'] = 'geometry={\n{{CARTESIAN}}\n}\nrhf\n'
    old = vrc_graph.select_correction_graph(tmp_path, original)
    modified = deepcopy(original)
    if change == 'command':
        modified['tasks'][0]['command'].append('--different-method')
    elif change == 'charge':
        modified['molecule']['charge'] = 2
    else:
        modified['molecule']['positions'][1][2] += .1
    assert vrc_graph.select_correction_graph(tmp_path, modified) != old
