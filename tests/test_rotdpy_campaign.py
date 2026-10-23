import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import rotdpy, rotdpy_campaign as campaign
from kinbot.execution import atomic_json, file_sha256


@pytest.fixture
def inp(tmp_path, monkeypatch):
    p = tmp_path/'surface.py'
    correction = tmp_path/'corr_surface.json'
    correction.write_text('{}')
    p.write_text('validation_metadata = '+repr({'name': 'production', 'correction': {
        'source': str(correction), 'source_sha256': file_sha256(correction)}})+'\n')
    monkeypatch.setattr('kinbot.vrc_scan_recovery.require_ready_for_rates', lambda p: None)
    return p


def test_launch_once_adopt_and_live_retry_preserve_review_and_contract(inp, monkeypatch):
    calls = []
    monkeypatch.setattr(campaign.subprocess, 'Popen', lambda *a, **kw:
                        calls.append(a) or SimpleNamespace(pid=42))
    monkeypatch.setattr('kinbot.pes_phase.process_identity', lambda pid: {'pid': pid, 'start': 'same'})
    record = campaign.launch(inp)
    assert len(calls) == 1
    assert campaign.launch(inp) == record
    assert campaign.launch(inp, retry=True) == record  # live owner never duplicated
    campaign.accept(inp, 'Native curve and state continuity reviewed')
    campaign.require_review(inp)
    monkeypatch.setattr(campaign, 'active', lambda *args: False)
    assert campaign.launch(inp)['launches'] == 1  # exited owner is not a retry loop
    assert campaign.launch(inp, retry=True)['launches'] == 2
    campaign.require_review(inp)
    with pytest.raises(RuntimeError, match='differs'):
        rotdpy.write_input(inp, inp.read_text()+'# changed\n')
    (inp.parent/'corr_surface.json').write_text('{"changed": true}')
    with pytest.raises(RuntimeError, match='correction changed'):
        campaign.launch(inp, retry=True)


def test_sampling_completion_still_required_without_a_second_correction_vote(inp, monkeypatch):
    atomic_json(campaign.receipt(inp), dict(schema=1, **campaign.contract(inp)))
    campaign.require_review(inp)
    with pytest.raises(RuntimeError, match='did not write'):
        rotdpy.number_of_states_file(inp)
    campaign.accept(inp, 'Reviewed')
    campaign.require_review(inp)
    # Approval is not a QC completion verdict.
    with pytest.raises(RuntimeError, match='did not write'):
        rotdpy.number_of_states_file(inp)


def test_accept_fences_the_reviewed_contract_under_its_lock(inp):
    current = campaign.contract(inp)
    atomic_json(campaign.receipt(inp), dict(schema=1, **current))
    before = campaign.receipt(inp).read_bytes()
    with pytest.raises(RuntimeError, match='changed after review'):
        campaign.accept(inp, 'Reviewed an older input', expected={**current, 'input_sha256':'old'})
    assert campaign.receipt(inp).read_bytes() == before
    campaign.accept(inp, 'Reviewed exact input', expected=current)
    campaign.require_review(inp)


def test_launch_failure_is_visible_and_requires_explicit_retry(inp, monkeypatch):
    def fail(*args, **kwargs):
        raise OSError('Unable to start executor')
    monkeypatch.setattr(campaign.subprocess, 'Popen', fail)
    with pytest.raises(OSError, match='Unable'):
        campaign.launch(inp)
    record = json.loads(campaign.receipt(inp).read_text())
    assert record['status'] == 'failed' and record['error'] == 'Unable to start executor'
    assert campaign.launch(inp) == record
    before = {p: p.read_bytes() for p in inp.parent.iterdir() if p.is_file()}
    assert not campaign.active(inp, record)
    assert before == {p: p.read_bytes() for p in inp.parent.iterdir() if p.is_file()}


def test_join_keeps_retry_window_open_without_resubmitting(inp, monkeypatch):
    atomic_json(campaign.receipt(inp), dict(schema=1, **campaign.contract(inp)))
    atomic_json(rotdpy.execution_path(inp), {'status': 'failed', 'error': 'native failure'})
    monkeypatch.setattr(campaign, 'active', lambda *args: False)
    monkeypatch.setattr('kinbot.live_control.read', lambda *args: {'hold_task_failures': True})
    monkeypatch.setattr(rotdpy, 'read_result', lambda p: {'status': 'complete'})
    waits = []
    def finish(seconds):
        waits.append(seconds)
        atomic_json(rotdpy.execution_path(inp), {'status': 'complete'})
        campaign.accept(inp, 'Reviewed while sampling')
    monkeypatch.setattr(campaign.time, 'sleep', finish)
    assert campaign.ensure(inp, interval=.01)['status'] == 'complete'
    assert len(waits) == 1


def test_input_write_is_idempotent_and_preserves_bound_samples(inp):
    original = inp.read_text()
    atomic_json(rotdpy.execution_path(inp), {'status': 'running'})
    rotdpy.write_input(inp, original)
    with pytest.raises(RuntimeError, match='differs'):
        rotdpy.write_input(inp, '# different\n')
    assert inp.read_text() == original


def test_scan_hold_defers_one_resume_and_owner_staging_consumes_it(tmp_path, monkeypatch):
    from kinbot import vrc_scan_recovery as scans
    inp = tmp_path/'surface.py'
    correction = tmp_path/'corr_surface.json'
    atomic_json(correction, {})
    inp.write_text('validation_metadata = '+repr(dict(name='production', correction={
        'source':str(correction), 'source_sha256':file_sha256(correction)}))+'\n')
    saved = dict(schema=1, input=str(inp), launches=2, **campaign.contract(inp))
    atomic_json(campaign.receipt(inp), saved)
    review = tmp_path/'scan_review/surface/scan_review.json'
    evidence = dict(context_sha256='context', points={'point':dict(status='accepted')})
    record = dict(reaction='surface', **evidence, evidence_sha256=scans.digest(evidence),
                  status='running', complete_grid=False)
    atomic_json(review, record)
    calls = []
    monkeypatch.setattr(campaign.subprocess, 'Popen', lambda *args, **kw:
                        calls.append(args) or SimpleNamespace(pid=42))
    monkeypatch.setattr(campaign, 'active', lambda *args: False)
    monkeypatch.setattr('kinbot.pes_phase.process_identity', lambda pid: {'pid':pid})
    campaign.validate(inp, saved, require_ready=False)
    with pytest.raises(scans.ScanReviewRequired):
        campaign.validate(inp, saved)
    with pytest.raises(scans.ScanReviewRequired):
        campaign.launch(inp, retry=True)  # An ordinary retry never waives a hold.
    held = campaign.launch(inp, retry=True, defer_scan_review=True)
    assert held['resume_requested'] and held['launches'] == 2 and not calls
    assert campaign.launch(inp, retry=True, defer_scan_review=True) == held
    with pytest.raises(scans.ScanReviewRequired):
        campaign.require_review(inp)
    record.update(status='resolved', complete_grid=True)
    atomic_json(review, record)
    # A finished scan without its bound correction is not permission to launch.
    with pytest.raises(ValueError, match='not bound'):
        campaign.launch(inp)
    assert not calls
    scans.bind_correction_review(correction)
    resumed = campaign.launch(inp)  # Replacement owner's normal staging call.
    assert len(calls) == 1 and resumed['launches'] == 3
    assert not resumed.get('resume_requested') and not resumed.get('resume_hold')
    campaign.launch(inp)
    assert len(calls) == 1
    correction.write_text('{"changed": true}')
    with pytest.raises(RuntimeError, match='correction changed'):
        campaign.validate(inp, resumed, require_ready=False)


def test_operator_rejection_is_possible_during_scan_hold(inp, tmp_path, monkeypatch):
    from kinbot import run_policy, run_walltime, vrc_scan_recovery as scans
    saved = dict(schema=1, input=str(inp), **campaign.contract(inp))
    atomic_json(campaign.receipt(inp), saved)
    def held(*args):
        raise scans.ScanReviewRequired([{'reaction':inp.stem}])
    monkeypatch.setattr(scans, 'require_ready_for_rates', held)
    monkeypatch.setattr(campaign, 'active', lambda *args: False)
    monkeypatch.setattr(run_walltime, 'jobs', lambda *args: [])
    run_policy.update(tmp_path, max_jobs=2)
    assert campaign.stop_for_rejection(inp, tmp_path/'archive', dict(id='decision', reason='Reject path'))
    assert json.loads(campaign.receipt(inp).read_text())['operator_rejected']['request'] == 'decision'
    campaign.launch(inp, retry=True)
    with pytest.raises(RuntimeError, match='rejected by the operator'):
        campaign.require_review(inp)


def test_early_leaf_and_final_pes_generate_one_identical_campaign(tmp_path, monkeypatch):
    from kinbot.pes import create_rotdpy_inputs
    from tests.test_rotdpy_profiles import sulfur_correction
    worker = tmp_path/'parent'
    worker.mkdir()
    par, data, correction = sulfur_correction(worker)
    channel = 'parent_hom_sci_1_2'
    leaf_name = channel+'__single_m2_r1'
    correction = correction.with_name('corr_'+leaf_name+'.json')
    data['frags_geom'] = [[[0., 0., 0.], [1.08, 0., 0.],
                          [-.54, .935, 0.], [-.54, -.935, 0.]], [[0., 0., 0.]]]
    data['surface_parameters'] = {'vrc_tst_state_model_file': str(worker/'ch3s_state_model.json')}
    atomic_json(correction, data)
    state = dict(reaction=leaf_name, correction=str(correction), correction_sha256=file_sha256(correction))
    par.update(pes=1, rotdpy_run=1, rotdpy_directory='sampling')
    species = lambda name: SimpleNamespace(name=name, wellorts=1)
    vts = SimpleNamespace(par=par, well=species('parent'), scan_reac={channel:
        SimpleNamespace(products=[species('methyl'), species('sulfur')])})
    calls = []
    monkeypatch.setattr(campaign.subprocess, 'Popen', lambda *args, **kw:
                        calls.append((args, kw)) or SimpleNamespace(pid=42))
    monkeypatch.setattr('kinbot.pes_phase.process_identity', lambda pid: {'pid': pid, 'start': 'same'})
    monkeypatch.setenv('KINBOT_PES_PHASE_ROOT', str(tmp_path))
    monkeypatch.chdir(worker)
    campaign.ready_leaf(vts, channel, state)
    inp = tmp_path/'sampling'/f'{leaf_name}.py'
    original = inp.read_bytes()
    assert len(calls) == 1 and calls[0][1]['cwd'] == tmp_path/'sampling'
    campaign.accept(inp, 'Native correction reviewed without changing it')
    atomic_json(worker/'vrctst'/f'corr_{channel}.json', dict(kind='multistate',
        reaction=channel, bottlenecks=[dict(name='single', states=[state])]))
    # The final aggregator runs one directory above the early worker.
    monkeypatch.chdir(tmp_path)
    joined = []
    monkeypatch.setattr(campaign, 'ensure', lambda p: joined.append(Path(p).resolve()))
    create_rotdpy_inputs(par, [['parent', channel, ['methyl', 'sulfur'], 0.]], [])
    assert joined == [inp] and len(calls) == 1
    assert inp.read_bytes() == original
    campaign.require_review(inp)
