"""Native-diagnostic VPT2 recovery is shared by every graph owner."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot.anl import dispatch, live_control
from kinbot import run_recover, run_control
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def graph(tmp_path, *, symmetry='Symmetry=None', freq='Freq=Anharmonic', route='#p', linear=False):
    spec = dispatch_spec()
    task = deepcopy(next(t for t in spec['tasks'] if t['id'] == 'gaussian_vpt2'))
    task.update(id='custom_anharmonic')
    task.pop('depends_on')
    task['input_template'] = task['input_template'].replace('NoSymm', symmetry).replace('Freq=Anharmonic', freq).replace('#p', route)
    spec['tasks'] = [spec['tasks'][0], task, molpro_task('live', 'rhf\n', geometry_from='initial'),
                     molpro_task('done', 'rhf\n', geometry_from='initial')]
    if linear:
        spec['molecule'] = dict(symbols=['H', 'C', 'N'], positions=[[0,0,-1.06],[0,0,0],[0,0,1.15]], charge=0, multiplicity=1)
    root = dispatch.prepare(spec, tmp_path/'graph')
    _mock_execution(root, 'l2_geometry', geometry='final.xyz')
    dispatch.advance(root)
    return root


def fail(root, text=None, *, legacy_linear=False, status='failed'):
    d = root/'tasks/custom_anharmonic'
    state = dispatch._load(root)[2]
    entry = state['tasks']['custom_anharmonic']
    entry.update(status=status, job_id='42')
    if legacy_linear:
        p = d/'vpt2.com'
        p.write_text(p.read_text().replace('NoSymm', 'Symmetry=None'))
        entry['input_sha256'] = file_sha256(p)
    text = text or ('ERROR: Inconsistency found in framework group definition:\n'
                    'New: CS - Old:\nError termination via Lnk1e in l717.exe\n')
    (d/'vpt2.log').write_text(text)
    atomic_json(d/'execution.json', dict(schema=1, task_id='custom_anharmonic',
        geometry_sha256=entry['geometry_sha256'], status='failed', slurm_job_id='42',
        failure_phase='execution', error='Program exited with status 1.'))
    atomic_json(root/'state.json', state)
    return d


@pytest.mark.parametrize('group', ['CS', 'C1', 'C2V', 'TD'])
@pytest.mark.parametrize('syntax', ['Symmetry=None', 'Symm=(None)'])
def test_common_dispatch_recovers_any_task_and_native_group(tmp_path, monkeypatch, group, syntax):
    root = graph(tmp_path, symmetry=syntax, freq='Frequencies=(Anharmonic,VibRot)')
    diagnostic = f'ERROR: Inconsistency found in framework group definition:\nNew: {group} - Old:\nError termination via Lnk1e in l717.exe\n'
    d = fail(root, diagnostic)
    geometry = (d/'geometry.xyz').read_bytes()
    _mock_execution(root, 'done')
    state = dispatch._load(root)[2]
    state['tasks']['live'].update(status='submitted', job_id='43')
    atomic_json(root/'state.json', state)
    hashes = {p:file_sha256(p) for p in (root/'tasks/live').iterdir() if p.is_file()}
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job == '43')
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['attempt'] == 2
    assert state['tasks']['custom_anharmonic']['status'] == 'staged'
    assert f'Symmetry=(PG={group})' in (d/'vpt2.com').read_text()
    assert 'Frequencies=(Anharmonic,VibRot)' in (d/'vpt2.com').read_text()
    assert (d/'geometry.xyz').read_bytes() == geometry
    assert state['tasks']['done']['status'] == 'complete'
    assert state['tasks']['live']['job_id'] == '43'
    assert hashes == {p:file_sha256(p) for p in hashes}
    # The shared C1 fallback is tried once if the native assignment was not C1.
    fail(root, diagnostic)
    if group != 'C1':
        assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
        assert 'Symmetry=(PG=C1)' in (d/'vpt2.com').read_text()
        fail(root, diagnostic)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'NoSymm' in (d/'vpt2.com').read_text()
    fail(root, diagnostic)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert len(list((root/'attempts').glob('live_*'))) == (2 if group == 'C1' else 3)
    assert (d/'geometry.xyz').read_bytes() == geometry
    assert hashes == {p:file_sha256(p) for p in hashes}


@pytest.mark.parametrize('route', ['#', '#P', '#n', '#t'])
def test_compound_frequency_without_symmetry_gets_route_keyword(tmp_path, monkeypatch, route):
    root = graph(tmp_path, symmetry='', freq='Freq=(Anharmonic,VibRot)', route=route)
    d = fail(root, 'ERROR: Inconsistency found in framework group definition:\nNew: CS - Old: C1\nError termination via Lnk1e\n')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    dispatch.advance(root)
    assert 'Symmetry=None' in (d/'vpt2.com').read_text()
    assert 'Freq=(Anharmonic,VibRot)' in (d/'vpt2.com').read_text()


@pytest.mark.parametrize('stage', ['failed', 'staged'])
def test_linear_native_fallback_reaches_common_dispatcher(tmp_path, monkeypatch, stage):
    root = graph(tmp_path, linear=True, freq='Freq=(Anharmonic,VibRot)')
    d = fail(root, 'Molecule is non-linear but NTrRo=5\nError termination via Lnk1e in l717.exe\n', legacy_linear=True, status=stage)
    if stage == 'staged':
        (d/'execution.json').unlink()
        state = dispatch._load(root)[2]
        state['tasks']['custom_anharmonic'].pop('job_id')
        atomic_json(root/'state.json', state)
    geometry = (d/'geometry.xyz').read_bytes()
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'NoSymm' in (d/'vpt2.com').read_text()
    assert (d/'geometry.xyz').read_bytes() == geometry
    assert 'gaussian_vpt2_linear_framework_v1' in state['native_interface_repairs']['custom_anharmonic']


@pytest.mark.parametrize('diagnostic', ['Error: segmentation violation', 'SCF convergence failure', 'Molecule is non-linear but NTrRo=5'])
def test_unmatched_or_nonlinear_failure_is_not_relabelled(tmp_path, monkeypatch, diagnostic):
    root = graph(tmp_path)
    fail(root, diagnostic+'\nError termination via Lnk1e\n')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert not (root/'.live_control_transaction.json').exists()


def test_recovery_resumes_interrupted_manifest_transaction(tmp_path, monkeypatch):
    root = graph(tmp_path); fail(root)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    original = run_recover.atomic_json
    def interrupt(path, value):
        if Path(path) == root/'state.json':
            raise OSError('interrupt')
        return original(path, value)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='interrupt'):
        dispatch.advance(root)
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['attempt'] == 2
    assert 'Symmetry=(PG=CS)' in (root/'tasks/custom_anharmonic/vpt2.com').read_text()
    assert len(list((root/'attempts').glob('live_*'))) == 1


@pytest.mark.parametrize('block', ['active', 'legacy', 'unknown'])
def test_allocation_fences_prevent_recovery_race(tmp_path, monkeypatch, block):
    root = graph(tmp_path); d = fail(root)
    before = (d/'vpt2.com').read_bytes()
    if block == 'legacy':
        state = dispatch._load(root)[2]
        state['tasks']['live'].update(status='submitted', job_id='43')
        atomic_json(root/'state.json', state)
        p = root/'tasks/live/task.json'; record=json.loads(p.read_text())
        record.pop('live_control_schema'); atomic_json(p, record)
        state['tasks']['live']['task_sha256']=file_sha256(p)
        atomic_json(root/'state.json', state)
    def active(job):
        if block == 'unknown':
            raise dispatch.SchedulerUnavailableError('queue unavailable')
        return block == 'active' or job == '43'
    monkeypatch.setattr(dispatch, '_job_active', active)
    if block == 'unknown':
        with pytest.raises(dispatch.SchedulerUnavailableError):
            dispatch.advance(root)
    else:
        dispatch.advance(root)
    assert (d/'vpt2.com').read_bytes() == before
    assert not (root/'.live_control_transaction.json').exists()


def test_reconcile_is_previewed_non_submitting_and_respects_live_owner(tmp_path, monkeypatch):
    import fcntl
    root = graph(tmp_path); fail(root)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    before = (root/'state.json').read_bytes()
    assert run_control.reconcile(tmp_path)['graphs'][0]['status'] == 'reviewed'
    assert (root/'state.json').read_bytes() == before
    with (root/'drive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert run_control.reconcile(tmp_path, apply=True)['graphs'][0]['status'] == 'owned_by_controller'
        assert (root/'state.json').read_bytes() == before
    def no_submit(*a, **kw):
        raise AssertionError('No scheduler mutation permitted')
    monkeypatch.setattr(dispatch.subprocess, 'run', no_submit)
    report = run_control.reconcile(tmp_path, ['graph'], apply=True)['graphs'][0]
    assert report['status'] == 'reconciled'
    assert report['changes']['custom_anharmonic']['after'] == 'staged'
    assert run_control.reconcile(tmp_path, ['no_match'], apply=True)['graphs'] == []


@pytest.mark.parametrize('frequency', ['Freq=(Anharmonic,VibRot)', 'Frequencies=(VibRot,Anharmonic)', 'Freq(Anharmonic)', 'Freq=Anharm'])
def test_output_reader_accepts_same_compound_options_as_input(frequency):
    from tests.test_anl_results import _GAUSSIAN_VPT2
    from kinbot.anl.results import parse_gaussian_vpt2
    expected = parse_gaussian_vpt2(_GAUSSIAN_VPT2, method='B3LYP', basis='cc-pVTZ')
    result = parse_gaussian_vpt2(_GAUSSIAN_VPT2.replace('Freq=Anharmonic', frequency), method='B3LYP', basis='cc-pVTZ')
    assert result == expected
    with pytest.raises(ValueError, match='requested Freq=Anharmonic'):
        parse_gaussian_vpt2(_GAUSSIAN_VPT2.replace('Freq=Anharmonic', 'Freq=Harmonic'), method='B3LYP', basis='cc-pVTZ')


def test_general_framework_recovery_ladder_is_bounded(tmp_path, monkeypatch):
    root = graph(tmp_path, symmetry='Symmetry=PG=C1')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    d = fail(root, 'ERROR: Inconsistency found in framework group definition:\nNew: TD - Old: T\nError termination via Lnk1e\n')
    dispatch.advance(root)
    assert 'Symmetry=None' in (d/'vpt2.com').read_text()
    fail(root)
    dispatch.advance(root)
    assert 'Symmetry=(PG=CS)' in (d/'vpt2.com').read_text()
    fail(root)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'NoSymm' in (d/'vpt2.com').read_text()
    fail(root)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'failed'
    assert state['tasks']['custom_anharmonic']['attempt'] == 4


def test_unrecognized_user_symmetry_requires_manual_review_without_abandoning_siblings(tmp_path, monkeypatch):
    root = graph(tmp_path, symmetry='Symmetry=(Loose,Follow)')
    d = fail(root)
    _mock_execution(root, 'done')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'failed'
    assert state['tasks']['done']['status'] == 'complete'
    assert not (root/'.live_control_transaction.json').exists()


def test_reconcile_isolates_one_graph_error(tmp_path, monkeypatch):
    bad = graph(tmp_path/'bad'); fail(bad)
    good = graph(tmp_path/'good'); fail(good)
    original = dispatch.advance
    def advance(path, **kwargs):
        if path == bad:
            raise RuntimeError('bad provenance')
        assert kwargs == {'submit': False}
        return original(path, **kwargs)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    monkeypatch.setattr(dispatch, 'advance', advance)
    records = run_control.reconcile(tmp_path, apply=True)['graphs']
    assert [(r['graph'], r['status']) for r in records] == [('bad/graph', 'blocked'), ('good/graph', 'reconciled')]


def test_mismatched_allocation_cannot_trigger_recovery(tmp_path, monkeypatch):
    root = graph(tmp_path); d = fail(root)
    p = d/'execution.json'; outcome = json.loads(p.read_text())
    outcome['slurm_job_id'] = '99'; atomic_json(p, outcome)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    with pytest.raises(RuntimeError, match='different allocation'):
        dispatch.advance(root)
    assert not (root/'.live_control_transaction.json').exists()


@pytest.mark.parametrize('new,old', [('C2V', 'C2'), ('C3V', 'C3'), ('C5V', 'C5')])
def test_native_detection_disagreement_gets_one_loose_retry(tmp_path, monkeypatch, new, old):
    root = graph(tmp_path)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    d = fail(root, f'ERROR: Inconsistency found in framework group definition:\nNew: {new} - Old:\nError termination via Lnk1e\n')
    dispatch.advance(root)
    geometry = (d/'geometry.xyz').read_bytes()
    diagnostic = (f'Full point group                 {old}      NOp 2\n'
                  f'Framework definition based on atomic masses: {new}\n'
                  'ERROR: Inconsistency found in framework group definition:\n'
                  f'New: {new} - Old: {old}\nError termination via Lnk1e\n')
    fail(root, diagnostic)
    state = dispatch.advance(root)
    assert state['tasks']['custom_anharmonic']['status'] == 'staged'
    assert f'Symmetry=(PG={new},Loose)' in (d/'vpt2.com').read_text()
    assert (d/'geometry.xyz').read_bytes() == geometry
    task = next(t for t in dispatch._load(root)[1]['tasks'] if t['id'] == 'custom_anharmonic')
    assert task['recovery']['symmetry_detection'] == 'loose'
    fail(root, diagnostic)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'Symmetry=(PG=C1)' in (d/'vpt2.com').read_text()
    assert (d/'geometry.xyz').read_bytes() == geometry
    fail(root, diagnostic)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'NoSymm' in (d/'vpt2.com').read_text()
    fail(root, diagnostic)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert len(list((root/'attempts').glob('live_*'))) == 4


def test_loose_retry_requires_both_native_framework_assignments(tmp_path, monkeypatch):
    root = graph(tmp_path)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    fail(root)
    dispatch.advance(root)
    fail(root, 'ERROR: Inconsistency found in framework group definition:\nNew: CS - Old: C1\nError termination via Lnk1e\n')
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    text = (root/'tasks/custom_anharmonic/vpt2.com').read_text()
    assert 'Loose' not in text and 'Symmetry=(PG=C1)' in text


def test_legacy_loose_receipt_uses_untried_general_fallback():
    from kinbot.anl.vpt2_recovery import plan, VPT2RecoveryUnavailable
    task = dict(id='any_name', input_template='#p Freq=Anharmonic Symmetry=(PG=C3V,Loose)',
                recovery=dict(assigned_native_framework=True,
                              attempted_framework_groups=['CS', 'C3V']))
    output = ('Full point group C3 NOp 3\nFramework definition based on atomic masses: C3V\n'
              'ERROR: Inconsistency found in framework group definition:\n'
              'New: C3V - Old: C3\nError termination via Lnk1e\n')
    assert plan(task, output)[0] == 'Symmetry=(PG=C1)'
    assert plan(task, output, history=['PG=C1'])[0] == 'NoSymm'
    with pytest.raises(VPT2RecoveryUnavailable, match='already attempted'):
        plan(task, output, history=['PG=C1', 'NOSYMM'])


def test_manual_retry_archive_prevents_automatic_repeat(tmp_path, monkeypatch):
    from kinbot import live_control as control
    from tests.test_live_control import allocation
    from kinbot.anl.vpt2_recovery import context
    control.configure(tmp_path)
    root = graph(tmp_path, symmetry='Symmetry=(PG=C1)')
    d = fail(root)
    geometry = (d/'geometry.xyz').read_bytes()
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    monkeypatch.setattr(live_control, 'slurm_job_record', lambda job: allocation(d, job))
    # An explicit operator retry predates the new automatic strategy receipt.
    request = control.request(tmp_path, d, reason='Reviewed symmetry retry', fallback=dict(schema=1, input_edits=[
        dict(old='Symmetry=(PG=C1)', new='Symmetry=None')]))
    live_control.process(root)
    dispatch.advance(root)
    fail(root)
    dispatch.advance(root)  # Empty old group -> native CS.
    _, spec, state = dispatch._load(root)
    task = next(t for t in spec['tasks'] if t['id'] == 'custom_anharmonic')
    history, _ = context(root, task, state['tasks']['custom_anharmonic'])
    assert {'PG=C1', 'NONE', 'PG=CS'} <= set(history)
    assert 'Symmetry=(PG=CS)' in (d/'vpt2.com').read_text()
    # Automatic edits preserve the manual fallback's original scientific binding.
    from kinbot.manual_fallback import logical_task
    assert 'Symmetry=(PG=C1)' in logical_task(task)['input_template']
    fail(root)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'NoSymm' in (d/'vpt2.com').read_text()
    fail(root)
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert (d/'geometry.xyz').read_bytes() == geometry
    assert json.loads((tmp_path/'.kinbot_controls/requests'/f'{request["id"]}.json').read_text())['status'] == 'applied'
    assert len(list((root/'attempts').glob('live_*'))) == 3


def test_native_success_after_general_fallback_cannot_trigger_another_retry(tmp_path, monkeypatch):
    from tests.test_anl_results import _GAUSSIAN_VPT2
    root = graph(tmp_path)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    d = fail(root)
    dispatch.advance(root)  # native CS
    fail(root)
    dispatch.advance(root)  # C1
    # A native-success parser recovery must run before any further symmetry retry.
    fail(root, _GAUSSIAN_VPT2)
    assert not dispatch.can_recover_gaussian_vpt2_symmetry(
        root, next(t for t in dispatch._load(root)[1]['tasks'] if t['id'] == 'custom_anharmonic'),
        dispatch._load(root)[2]['tasks']['custom_anharmonic'])
    assert len(list((root/'attempts').glob('live_*'))) == 2


def test_linear_framework_failure_does_not_cycle_none_or_c1(tmp_path, monkeypatch):
    root = graph(tmp_path, linear=True)
    d = fail(root, 'ERROR: Inconsistency found in framework group definition:\n'
             'New: C*V - Old: C1\nError termination via Lnk1e\n')
    monkeypatch.setattr(dispatch, '_job_active', lambda _: False)
    assert 'NoSymm' in (d/'vpt2.com').read_text()
    # The unchanged template says None, but rendering has already tried NoSymm.
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'staged'
    assert 'Symmetry=(PG=C*V)' in (d/'vpt2.com').read_text()
    fail(root, 'ERROR: Inconsistency found in framework group definition:\n'
         'New: C*V - Old: C1\nError termination via Lnk1e\n')
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert len(list((root/'attempts').glob('live_*'))) == 1


def test_route_edit_preserves_title_and_refuses_unrecognized_controls():
    from kinbot.anl.vpt2_recovery import replace_symmetry, VPT2RecoveryUnavailable
    source = '#p Freq=(Anharmonic,VibRot)\n Symm=(PG=C3V, Loose)\n\nTitle Symmetry=None\n\n0 1\n'
    changed = replace_symmetry(source, 'Symmetry=(PG=C1)')
    assert changed == source.replace('Symm=(PG=C3V, Loose)', 'Symmetry=(PG=C1)')
    for value in ('Symmetry=(PG=C3V,Follow)', 'NoSymm Symmetry=None'):
        with pytest.raises(VPT2RecoveryUnavailable):
            replace_symmetry('#p Freq=Anharmonic ' + value, 'Symmetry=None')


@pytest.mark.parametrize('new,old', [('C2V', 'C2'), ('C3V', 'C3')])
def test_persistent_eckart_mismatch_orientation_retries_are_bounded(tmp_path, monkeypatch, new, old):
    root = graph(tmp_path, symmetry='Symmetry=(PG=C1)')
    _, spec, state = dispatch._load(root)
    task = next(t for t in spec['tasks'] if t['id'] == 'custom_anharmonic')
    task['recovery'] = dict(attempted_symmetry_settings=[
        'NONE', f'PG={new}', f'PG={new},LOOSE', 'PG=C1'])
    # Restage the synthetic pre-existing history through the same transaction.
    live_control.commit(root, spec, state, {'custom_anharmonic'}, {'custom_anharmonic': 'fixture'})
    dispatch.advance(root)
    _mock_execution(root, 'done')
    state = dispatch._load(root)[2]
    state['tasks']['live'].update(status='submitted', job_id='43')
    atomic_json(root/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job == '43')
    native = (f'Symm: InNGrp=C1 InAxis=0\nFull point group {old} NOp 2\n'
              f'Framework definition based on atomic masses: {new}\n'
              'NOTE: The system is set in Eckart orientation for the anharmonic treatment.\n'
              'ERROR: Inconsistency found in framework group definition:\n'
              f'New: {new} - Old: {old}\nError termination via Lnk1e in l717.exe\n')
    d = fail(root, native)
    geometry = (d/'geometry.xyz').read_bytes()
    siblings = {p: file_sha256(p) for n in ('live','done')
                for p in (root/'tasks'/n).iterdir() if p.is_file()}
    for keyword in ('NoSymm', f'Symmetry=(PG={new},Loose,CenterOfMass)'):
        state = dispatch.advance(root)
        assert state['tasks']['custom_anharmonic']['status'] == 'staged'
        assert keyword in (d/'vpt2.com').read_text()
        assert (d/'geometry.xyz').read_bytes() == geometry
        assert siblings == {p: file_sha256(p) for p in siblings}
        fail(root, native)
    before = len(list((root/'attempts').glob('live_*')))
    assert dispatch.advance(root)['tasks']['custom_anharmonic']['status'] == 'failed'
    assert len(list((root/'attempts').glob('live_*'))) == before


@pytest.mark.parametrize('missing', ['mass', 'full', 'eckart', 'success', 'linear'])
def test_orientation_retry_requires_corroborated_native_mismatch(missing):
    from kinbot.anl.vpt2_recovery import plan, VPT2RecoveryUnavailable
    task = dict(id='generic', input_template='#p Freq=Anharmonic Symmetry=(PG=C1)',
                recovery=dict(attempted_symmetry_settings=['NONE', 'PG=C2V', 'PG=C2V,LOOSE']))
    lines = dict(full='Full point group C2 NOp 2\n',
                 mass='Framework definition based on atomic masses: C2V\n',
                 eckart='The system is set in Eckart orientation for the anharmonic treatment.\n')
    lines.pop(missing, None)
    native = ''.join(lines.values()) + ('ERROR: Inconsistency found in framework group definition:\n'
                                       'New: C2V - Old: C2\nError termination via Lnk1e\n')
    if missing == 'success':
        native += 'Normal termination of Gaussian\n'
    with pytest.raises(VPT2RecoveryUnavailable):
        plan(task, native, history=['NOSYMM'], linear=missing == 'linear')


@pytest.mark.parametrize('new,old', [('CS', 'C1'), ('C2V', 'C2'), ('TD', '')])
def test_nonlinear_nosymm_does_not_require_mass_or_eckart_details(new, old):
    from kinbot.anl.vpt2_recovery import plan, VPT2RecoveryUnavailable
    task = dict(id='any_anharmonic', input_template='#p Freq=Anharmonic Symmetry=(PG=C1)',
                recovery=dict(attempted_symmetry_settings=['NONE', f'PG={new}']))
    output = ('ERROR: Inconsistency found in framework group definition:\n'
              f'New: {new} - Old: {old}\nError termination via Lnk1e\n')
    assert plan(task, output)[0] == 'NoSymm'
    task['input_template'] = '#p Freq=Anharmonic NoSymm'
    task['recovery']['attempted_symmetry_settings'].extend(['PG=C1', 'NOSYMM'])
    with pytest.raises(VPT2RecoveryUnavailable, match='exhausted'):
        plan(task, output)
    with pytest.raises(VPT2RecoveryUnavailable, match='not the framework-group'):
        plan(task, 'Error: segmentation violation\nError termination via Lnk1e\n')
