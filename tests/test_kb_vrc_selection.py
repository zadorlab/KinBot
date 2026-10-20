"""VRC auto-selection uses public ReactionFinder identities."""

import json
import os
from pathlib import Path
import threading
from types import SimpleNamespace

import pytest

import kinbot.kb as kb
import kinbot.anl.automatic as automatic
from kinbot.kb import (_accepted_homolytic_channels,
                       _check_initial_minimum_frequencies,
                       _configure_automatic_queue_admission,
                       _post_reaction_branches)


def test_completed_hom_sci_with_internal_hs_label_is_selected():
    scission = SimpleNamespace(
        instance_name='parent_hom_sci_1_2', family_name='hs',
        products=[object(), object()])
    tight = SimpleNamespace(
        instance_name='parent_h2_elim_4_7', family_name='h2elim',
        products=[object(), object()])
    species = SimpleNamespace(
        reac_obj=[scission, tight], reac_ts_done=[-1, -1],
        reac_type=['hom_sci', 'h2_elim'])

    assert _accepted_homolytic_channels(species) == [
        'parent_hom_sci_1_2']


def test_incomplete_or_nonbimolecular_scission_is_not_selected():
    incomplete = SimpleNamespace(
        instance_name='incomplete', family_name='hs',
        products=[object(), object()])
    one_product = SimpleNamespace(
        instance_name='one_product', family_name='hs', products=[object()])
    species = SimpleNamespace(
        reac_obj=[incomplete, one_product], reac_ts_done=[0, -1],
        reac_type=['hom_sci', 'hom_sci'])

    assert _accepted_homolytic_channels(species) == []


def test_initial_frequency_check_supports_atoms_and_two_electron_diatomics():
    assert _check_initial_minimum_frequencies([], 'H')
    h2 = [-12.]
    assert _check_initial_minimum_frequencies(h2, 'H2')
    assert h2 == [12.]


def test_initial_frequency_check_rejects_two_imaginary_modes():
    modes = [-5., -2., 100.]
    assert not _check_initial_minimum_frequencies(modes, 'test')
    assert modes[0] == 5.


def _branch_parameters(tmp_path):
    return {
        'queue_job_limit': 8,
        'vrc_tst_max_nodes': 8,
        'rotdpy_max_jobs': 8,
        'vrc_tst_scan': {'root': ['root_hom_sci_1_2']},
        'vrc_tst_noscan': {},
        'anl_automatic_directory': str(tmp_path / 'anl_automatic'),
    }


def test_automatic_run_configures_one_recursive_queue_scope(
        monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    for name in ('KINBOT_QUEUE_ROOT', 'KINBOT_QUEUE_ROOT_RECURSIVE',
                 'KINBOT_QUEUE_LOCK'):
        monkeypatch.delenv(name, raising=False)

    root = _configure_automatic_queue_admission({
        'queuing': 'slurm', 'queue_job_limit': 8}, True)

    assert root == tmp_path.resolve()
    assert os.environ['KINBOT_QUEUE_ROOT'] == str(tmp_path.resolve())
    assert os.environ['KINBOT_QUEUE_ROOT_RECURSIVE'] == '1'
    assert os.environ['KINBOT_QUEUE_LOCK'] == str(
        tmp_path.resolve() / '.kinbot_submission.lock')


def test_vrc_and_composite_branches_run_concurrently_with_disjoint_slots(
        monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    barrier = threading.Barrier(2)
    observed = {}

    def prepare(well, parameters, directory):
        observed['prepare_slots'] = parameters['queue_job_limit']
        return {'schema': 1, 'status': 'planned'}

    def run_vrc(well, parameters, qc):
        observed['vrc_slots'] = (
            parameters['vrc_tst_max_nodes'], parameters['rotdpy_max_jobs'])
        barrier.wait(timeout=2)
        return {'status': 'complete'}

    def run_auto(well, parameters, input_file, directory, plan):
        observed['auto_slots'] = parameters['queue_job_limit']
        barrier.wait(timeout=2)
        return {'status': 'complete', 'handoff_manifest': 'handoff.json',
                'thermochemistry_models': 'models.json'}

    monkeypatch.setattr(automatic, 'prepare_automatic_cbh_anl', prepare)
    monkeypatch.setattr(kb, '_execute_vrc_branch', run_vrc)
    monkeypatch.setattr(kb, '_execute_automatic_branch', run_auto)

    result = _post_reaction_branches(
        object(), _branch_parameters(tmp_path), object(), 'input.json', True)

    assert result['status'] == 'complete'
    assert observed == {
        'prepare_slots': 4, 'vrc_slots': (4, 4), 'auto_slots': 4}
    state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
    assert state['status'] == 'complete'
    assert state['mess_gate'] == 'ready'
    assert state['branches']['vrc_tst']['status'] == 'complete'
    assert state['branches']['cbh_anl']['status'] == 'complete'


def test_shared_pool_preserves_frozen_plan_but_reports_no_branch_reservations(monkeypatch, tmp_path):
    from kinbot import run_policy
    monkeypatch.chdir(tmp_path)
    run_policy.update(tmp_path, max_jobs=23)
    observed = []
    def prepare(well, parameters, directory):
        observed.append(parameters['queue_job_limit'])
        return {'schema': 1, 'status': 'planned'}
    monkeypatch.setattr(automatic, 'prepare_automatic_cbh_anl', prepare)
    monkeypatch.setattr(kb, '_execute_vrc_branch', lambda *args: {'status': 'complete'})
    monkeypatch.setattr(kb, '_execute_automatic_branch', lambda *args: {
        'status': 'complete', 'handoff_manifest': 'handoff.json', 'thermochemistry_models': 'models.json'})
    _post_reaction_branches(object(), _branch_parameters(tmp_path), object(), 'input.json', True)
    assert observed == [4]  # Resource contract in persisted chemistry plans stays compatible.
    state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
    assert state['slot_cap'] == 'live shared pool'
    assert all(b['slots'] == 'shared' for b in state['branches'].values())


def test_branch_restart_replaces_stale_failure_and_rejoins_completed_work(
        monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    Path('.kinbot_post_reaction_branches.json').write_text(json.dumps({
        'schema': 1, 'status': 'failed', 'mess_gate': 'blocked',
        'branches': {
            'vrc_tst': {'status': 'failed', 'error': 'old outage'},
            'cbh_anl': {'status': 'failed', 'error': 'old validator'},
        },
    }))
    monkeypatch.setattr(
        automatic, 'prepare_automatic_cbh_anl',
        lambda *args, **kwargs: {'schema': 1, 'status': 'planned'})
    monkeypatch.setattr(
        kb, '_execute_vrc_branch', lambda *args: {'status': 'complete'})
    monkeypatch.setattr(
        kb, '_execute_automatic_branch', lambda *args: {
            'status': 'complete', 'handoff_manifest': 'handoff.json',
            'thermochemistry_models': 'models.json'})

    result = _post_reaction_branches(
        object(), _branch_parameters(tmp_path), object(), 'input.json', True)
    state = json.loads(Path(
        '.kinbot_post_reaction_branches.json').read_text())
    assert result['status'] == 'complete'
    assert state['status'] == 'complete'
    assert state['mess_gate'] == 'ready'
    assert 'error' not in state['branches']['vrc_tst']
    assert 'error' not in state['branches']['cbh_anl']


@pytest.mark.parametrize('fails', [True, False])
def test_composite_outcome_is_visible_while_vrc_still_runs(monkeypatch, tmp_path, fails):
    import time
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(automatic, 'prepare_automatic_cbh_anl',
                        lambda *args, **kwargs: {'status': 'planned'})

    def auto(*args):
        if fails:
            raise RuntimeError('reference L2 failed')
        return {'status': 'complete'}

    def vrc(*args):
        deadline = time.monotonic() + 2.
        while time.monotonic() < deadline:
            state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
            if state['branches']['cbh_anl']['status'] != 'running':
                break
            time.sleep(.01)
        assert state['branches']['cbh_anl']['status'] == ('failed' if fails else 'complete')
        assert state['branches']['vrc_tst']['status'] == 'running'
        assert state['mess_gate'] == ('blocked' if fails else 'waiting')
        if fails:
            branch = state['branches']['cbh_anl']
            assert branch['error_type'] == 'RuntimeError'
            assert 'in auto' in branch['traceback']
            assert 'RuntimeError: reference L2 failed' in branch['traceback']
        return {'status': 'complete'}

    monkeypatch.setattr(kb, '_execute_automatic_branch', auto)
    monkeypatch.setattr(kb, '_execute_vrc_branch', vrc)
    if fails:
        with pytest.raises(RuntimeError, match='MESS was not started'):
            _post_reaction_branches(object(), _branch_parameters(tmp_path),
                                    object(), 'input.json', True)
    else:
        assert _post_reaction_branches(object(), _branch_parameters(tmp_path),
                                      object(), 'input.json', True)['status'] == 'complete'


def test_vrc_traceback_is_published_before_composite_finishes(monkeypatch, tmp_path):
    import time
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(automatic, 'prepare_automatic_cbh_anl',
                        lambda *args, **kwargs: {'status': 'planned'})

    def auto(*args):
        deadline = time.monotonic() + 2.
        while time.monotonic() < deadline:
            state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
            branch = state['branches']['vrc_tst']
            if branch['status'] == 'failed':
                assert state['branches']['cbh_anl']['status'] == 'running'
                assert branch['error_type'] == 'RecursionError'
                assert 'in fail_vrc' in branch['traceback']
                assert 'RecursionError: fixture VRC failure' in branch['traceback']
                return {'status': 'complete'}
            time.sleep(.01)
        pytest.fail('VRC traceback was hidden until composite completion')

    def fail_vrc(*args):
        raise RecursionError('fixture VRC failure')

    monkeypatch.setattr(kb, '_execute_automatic_branch', auto)
    monkeypatch.setattr(kb, '_execute_vrc_branch', fail_vrc)
    with pytest.raises(RuntimeError, match='MESS was not started'):
        _post_reaction_branches(object(), _branch_parameters(tmp_path),
                                object(), 'input.json', True)
    state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
    assert state['branches']['cbh_anl']['status'] == 'complete'


def test_mess_gate_stays_blocked_until_both_independent_branches_succeed(
        monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    completed = []
    monkeypatch.setattr(
        automatic, 'prepare_automatic_cbh_anl',
        lambda *args, **kwargs: {'schema': 1, 'status': 'planned'})

    def fail_vrc(*args):
        raise RuntimeError('VRC failed')

    def finish_auto(*args):
        completed.append('cbh_anl')
        return {'status': 'complete'}

    monkeypatch.setattr(kb, '_execute_vrc_branch', fail_vrc)
    monkeypatch.setattr(kb, '_execute_automatic_branch', finish_auto)

    with pytest.raises(RuntimeError, match='MESS was not started'):
        _post_reaction_branches(
            object(), _branch_parameters(tmp_path), object(), 'input.json', True)

    assert completed == ['cbh_anl']
    state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
    assert state['status'] == 'failed'
    assert state['mess_gate'] == 'blocked'
    assert state['branches']['vrc_tst']['status'] == 'failed'
    assert state['branches']['cbh_anl']['status'] == 'complete'


def test_missing_mess_does_not_prevent_either_independent_branch(
        monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    completed = []
    monkeypatch.setattr(
        automatic, 'validate_mess_environment',
        lambda parameters: (_ for _ in ()).throw(
            FileNotFoundError('mess unavailable')))
    monkeypatch.setattr(
        automatic, 'prepare_automatic_cbh_anl',
        lambda *args, **kwargs: {'schema': 1, 'status': 'planned'})
    monkeypatch.setattr(
        kb, '_execute_vrc_branch',
        lambda *args: completed.append('vrc'))
    monkeypatch.setattr(
        kb, '_execute_automatic_branch',
        lambda *args: completed.append('cbh_anl'))

    with pytest.raises(RuntimeError, match='MESS was not started'):
        _post_reaction_branches(
            object(), _branch_parameters(tmp_path), object(), 'input.json', True)

    assert sorted(completed) == ['cbh_anl', 'vrc']
    state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
    assert state['branches']['vrc_tst']['status'] == 'complete'
    assert state['branches']['cbh_anl']['status'] == 'complete'
    assert state['final_solver']['status'] == 'failed'
    assert state['mess_gate'] == 'blocked'


def test_automatic_preparation_failure_does_not_launch_unprepared_branch(
        monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    completed = []

    def fail_prepare(*args, **kwargs):
        raise RuntimeError('planning failed')

    monkeypatch.setattr(
        automatic, 'prepare_automatic_cbh_anl', fail_prepare)
    monkeypatch.setattr(
        kb, '_execute_vrc_branch',
        lambda *args: completed.append('vrc') or {'status': 'complete'})
    monkeypatch.setattr(
        kb, '_execute_automatic_branch',
        lambda *args: pytest.fail('unprepared composite branch was launched'))

    with pytest.raises(RuntimeError, match='MESS was not started'):
        _post_reaction_branches(
            object(), _branch_parameters(tmp_path), object(), 'input.json', True)

    assert completed == ['vrc']
    state = json.loads(Path('.kinbot_post_reaction_branches.json').read_text())
    assert state['status'] == 'failed'
    assert state['mess_gate'] == 'blocked'
    assert state['branches']['vrc_tst']['status'] == 'complete'
    assert state['branches']['cbh_anl']['status'] == 'failed'
    assert state['branches']['cbh_anl']['error'] == 'planning failed'
