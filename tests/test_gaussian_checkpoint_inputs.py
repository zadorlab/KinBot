"""Exercise real worker scratch staging and owner-driven attempt replacement."""
from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from kinbot import live_control, qc_scratch, run_recover
from kinbot.anl import checkpoint_inputs as checkpoints, dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.test_anl_dispatch import _mock_execution
from tests.test_r8_validation_panel import load, spectrum
from kinbot.anl.results import parse_gaussian_harmonic


@pytest.fixture
def graph(tmp_path, monkeypatch, request):
    panel = load('campaign')
    native = tmp_path/'native'
    native.write_text(f'#!{sys.executable}\n'
        'from pathlib import Path\nimport re,sys\n'
        'text=sys.stdin.read()\n'
        'old=re.search(r"(?m)^%oldchk=(.*)$",text)[1]\n'
        'assert Path(old).name == old\n'
        'assert Path(old).read_text() == "mock output\\n"\n'
        'assert "IRC=(RCFC,GradientOnly,EulerPC,MassWeighted,Forward" in text\n'
        'new=re.search(r"(?m)^%chk=(.*)$",text)[1]\n'
        'Path(new).write_text("new branch checkpoint")\n'
        'print("Normal termination of Gaussian")\n')
    native.chmod(0o755)
    root = tmp_path/'panel'; root.mkdir(); directory = root/'example'
    point = dict(symbols=['O','H','H'],positions=[[0.,0,0],[1.,0,0],[0,1.,0]],charge=0,multiplicity=1)
    spec = panel.native_irc('example', point, directory, step_sizes=(10,))
    for task in spec['tasks']:
        task['resources'] = dict(cores=2, memory_mb=8192, walltime='01:00:00', partition='test')
        task['command'] = [str(native)]
    monkeypatch.setattr(dispatch, 'assign_partitions', lambda spec: None)
    monkeypatch.setattr(dispatch, 'render_site_setup', lambda _: '# simulated site\n')
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job) == '77')
    monkeypatch.setattr(qc_scratch, 'candidates', lambda *a: iter([('test', tmp_path/'scratch')]))
    dispatch.prepare(spec, directory)
    _mock_execution(directory, 'hessian')
    p = directory/'tasks/hessian/execution.json'; e = json.loads(p.read_text())
    e['details']['parsed_result'] = parse_gaussian_harmonic(spectrum([-350.,800.,1000.]),
        method='B2PLYP',basis='cc-pVTZ',dispersion='GD3BJ',natoms=3,stationary_point='transition_state')
    atomic_json(p, e)
    dispatch.advance(directory)
    if getattr(request, 'param', True):
        live_control.configure(root)
    return directory


def old_failure(graph, *, job='42'):
    directory = graph/'tasks/forward_10'
    _, _, state = dispatch._load(graph)
    entry = state['tasks']['forward_10']; entry.update(status='failed', job_id=job)
    source = str(graph/'tasks/hessian/hessian.chk')
    atomic_json(directory/'execution.json', dict(schema=1, task_id='forward_10',
        geometry_sha256=entry['geometry_sha256'], slurm_job_id=job, status='failed',
        failure_phase='validation', error=f'Native QC file must be a relative basename: {source!r}'))
    state['tasks']['reverse_10'].update(status='submitted', job_id='77')
    atomic_json(graph/'state.json', state)


def test_real_worker_executes_local_checkpoint_and_keeps_logical_input(graph):
    directory = graph/'tasks/forward_10'
    parent = graph/'tasks/hessian/hessian.chk'
    original = (directory/'forward_10.com').read_bytes()
    digest = file_sha256(parent)
    result = dispatch.run_task(directory/'task.json')
    assert result['status'] == 'executed'
    proof = result['details']['checkpoint_inputs'][0]
    audit = directory/proof['executed_input']
    assert file_sha256(audit) == proof['executed_sha256']
    assert audit.read_text() == original.decode().replace(str(parent), proof['checkpoints'][0]['local'])
    assert (directory/'forward_10.com').read_bytes() == original
    assert file_sha256(parent) == digest
    assert proof['executed_input'] in result['artifacts']
    assert proof['checkpoints'][0]['local'] in result['artifacts']
    state = dispatch.advance(graph)
    assert state['tasks']['forward_10']['status'] == 'complete'
    assert dispatch.advance(graph) == state


def test_dependency_inspection_preserves_published_imports_and_staging_guard(graph):
    directory = graph/'tasks/forward_10'
    text = (directory/'forward_10.com').read_text()
    expected = checkpoints.inspect_dependencies(directory, 'forward_10.com', text)
    dispatch.run_task(directory/'task.json')
    before = {p: file_sha256(p) for p in directory.iterdir() if p.is_file()}
    assert checkpoints.inspect_dependencies(directory, 'forward_10.com', text) == expected
    with pytest.raises(ValueError, match='filename already exists'):
        checkpoints.plan(directory, 'forward_10.com', text)
    assert before == {p: file_sha256(p) for p in directory.iterdir() if p.is_file()}
    parent = graph/'tasks/hessian/hessian.chk'
    parent.write_text('changed accepted parent')
    with pytest.raises(RuntimeError, match='artifact hessian.chk changed'):
        checkpoints.inspect_dependencies(directory, 'forward_10.com', text)


def test_dangling_import_destination_is_not_overwritten(graph):
    directory = graph/'tasks/forward_10'
    text = (directory/'forward_10.com').read_text()
    binding, = checkpoints.inspect_dependencies(directory, 'forward_10.com', text)
    local = directory/binding['local']
    local.symlink_to(directory/'missing.chk')
    with pytest.raises(ValueError, match='filename already exists'):
        checkpoints.plan(directory, 'forward_10.com', text)
    assert local.is_symlink() and not local.exists()


def test_queued_repair_archives_failure_preserves_live_work_and_repeats(graph, monkeypatch):
    from kinbot.anl import live_control as owner
    old_failure(graph)
    monkeypatch.setattr(owner, 'slurm_job_record', lambda _: {})
    stable = {p:file_sha256(p) for ident in ('hessian', 'reverse_10')
              for p in (graph/'tasks'/ident).iterdir() if p.is_file()}
    logical = (graph/'workflow.json').read_bytes()
    assert checkpoints.retry_requests(graph)[0]['status'] == 'eligible'
    assert checkpoints.retry_requests(graph, apply=True)[0]['status'] == 'queued'
    assert checkpoints.retry_requests(graph, apply=True)[0]['status'] == 'already_requested'
    owner.process(graph)
    state = dispatch.advance(graph)
    retry = state['tasks']['forward_10']
    assert retry['status'] == 'staged' and retry['attempt'] == 2
    assert (Path(retry['previous_attempt'])/'execution.json').is_file()
    assert state['tasks']['reverse_10']['status'] == 'submitted'
    assert state['tasks']['reverse_10']['job_id'] == '77'
    assert stable == {p:file_sha256(p) for p in stable}
    assert (graph/'workflow.json').read_bytes() == logical
    assert checkpoints.retry_requests(graph, apply=True) == []
    dispatch.run_task(graph/'tasks/forward_10/task.json')
    state = dispatch.advance(graph)
    assert state['tasks']['forward_10']['status'] == 'complete'
    assert dispatch.advance(graph) == state
    assert len(list((graph/'attempts').glob('live_*'))) == 1


@pytest.mark.parametrize('defect', ['modified_parent', 'unaccepted_parent', 'wrong_source', 'symlink', 'wrong_task', 'native_failure'])
def test_invalid_checkpoint_or_failure_is_not_retried(graph, defect):
    old_failure(graph)
    parent = graph/'tasks/hessian/hessian.chk'
    if defect == 'modified_parent': parent.write_text('changed')
    elif defect == 'symlink':
        replacement = parent.with_suffix('.saved'); parent.rename(replacement); parent.symlink_to(replacement)
    elif defect == 'unaccepted_parent':
        _, _, state = dispatch._load(graph); state['tasks']['hessian']['status'] = 'submitted'
        atomic_json(graph/'state.json', state)
    elif defect in ('wrong_source', 'wrong_task'):
        p = graph/'tasks/forward_10/forward_10.com'
        p.write_text(p.read_text().replace(str(parent), str(graph/'tasks/other/hessian.chk')))
    else:
        p = graph/'tasks/forward_10/execution.json'; e=json.loads(p.read_text())
        e['error'] = 'Program exited with status 1.'; e['failure_phase']='execution'; atomic_json(p,e)
    report = checkpoints.retry_requests(graph, apply=True)
    assert not report or report[0]['status'] == 'blocked'
    assert not list((graph.parent/'.kinbot_controls/requests').glob('*.json'))


def test_retry_budget_and_queue_outage_preserve_failed_attempt(graph, monkeypatch):
    from kinbot.anl import live_control as owner
    old_failure(graph)
    checkpoints.retry_requests(graph, apply=True)
    monkeypatch.setattr(owner, 'slurm_job_record', lambda _: {})
    monkeypatch.setattr(dispatch, '_job_active', lambda _: (_ for _ in ()).throw(
        dispatch.SchedulerUnavailableError('temporary scheduler outage')))
    owner.process(graph)
    assert dispatch._load(graph)[2]['tasks']['forward_10']['status'] == 'failed'
    assert not (graph/'attempts').exists()
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job)=='77')
    owner.process(graph)
    dispatch.advance(graph)
    old_failure(graph, job='43')
    assert checkpoints.retry_requests(graph, apply=True)[0]['status'] == 'exhausted'


def test_interrupted_owner_transaction_resumes_without_second_retry(graph, monkeypatch):
    from kinbot.anl import live_control as owner
    old_failure(graph)
    monkeypatch.setattr(owner, 'slurm_job_record', lambda _: {})
    checkpoints.retry_requests(graph, apply=True)
    save = run_recover.atomic_json
    def interrupted(path, value):
        if Path(path) == graph/'state.json': raise OSError('transient publication race')
        return save(path, value)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupted)
    with pytest.raises(OSError, match='publication race'):owner.process(graph)
    monkeypatch.setattr(run_recover, 'atomic_json', save)
    state = dispatch.advance(graph)
    assert state['tasks']['forward_10']['attempt'] == 2
    assert state['tasks']['reverse_10']['job_id'] == '77'
    assert len(list((graph/'attempts').glob('live_*'))) == 1


@pytest.fixture
def resume_panel(graph, monkeypatch):
    resume = load('resume'); root = graph.parent
    spec = deepcopy(dispatch._load(graph)[1])
    monkeypatch.setattr(resume.campaign, 'graphs', lambda root: {'example':deepcopy(spec)})
    atomic_json(graph/'requested.json', spec)
    atomic_json(root/'stage.json', dict(sources=resume.campaign.fingerprint(), workflows=['example']))
    return resume, root, graph


@pytest.mark.parametrize('graph', [False], indirect=True)
def test_replacement_owner_applies_requests_and_adopts_allocations(resume_panel, monkeypatch):
    from kinbot.anl import live_control as owner, protocol
    from kinbot import run_control, run_policy
    from types import SimpleNamespace
    resume, root, graph = resume_panel
    old_failure(graph)
    # The handoff enables only the panel, even with an unrelated run tree next
    # to it. Exercise the actual CLI before queueing owner-driven repairs.
    unrelated = root.parent/'unrelated-checkpoints'; unrelated.mkdir()
    (unrelated/'retained').write_text('unrelated live work')
    scan = live_control.os.scandir
    visited = []
    def scoped_scan(path):
        visited.append(Path(path))
        assert Path(path).is_relative_to(root)
        return scan(path)
    with monkeypatch.context() as patch:
        patch.setattr(live_control.os, 'scandir', scoped_scan)
        run_control.main(['enable', str(root)])
        run_control.main(['enable', str(root)])
    assert visited
    assert not (root.parent/live_control.NAME).exists()
    checkpoints.retry_requests(root, apply=True)
    monkeypatch.setattr(owner, 'slurm_job_record', lambda _: {})
    monkeypatch.setattr(protocol, 'preflight_ready', lambda _: True)
    submitted = []
    def submit(command, **kwargs):
        kwargs['before_submit'](); submitted.append(str(kwargs['cwd']))
        return SimpleNamespace(returncode=0, stdout='88', stderr='')
    monkeypatch.setattr(run_policy, 'submit', submit)
    class Tick(Exception): pass
    monkeypatch.setattr(resume.time, 'sleep', lambda _: (_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick):resume.run(root)
    state = dispatch._load(graph)[2]
    assert state['tasks']['forward_10']['job_id'] == '88'
    assert state['tasks']['reverse_10']['job_id'] == '77'
    monkeypatch.setattr(dispatch, '_job_active', lambda _: True)
    with pytest.raises(Tick):resume.run(root)
    assert len(submitted) == 1
    monkeypatch.setenv('SLURM_JOB_ID', '88')
    dispatch.run_task(graph/'tasks/forward_10/task.json')
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job)=='77')
    with pytest.raises(Tick):resume.run(root)
    state = dispatch._load(graph)[2]
    assert state['tasks']['forward_10']['status'] == 'complete'
    assert state['tasks']['reverse_10']['job_id'] == '77'
    assert len(submitted) == 1
    assert not json.loads((root/'progress.json').read_text())['scientific_approval']
    assert (unrelated/'retained').read_text() == 'unrelated live work'


@pytest.mark.parametrize('race', ['none','submitting','staging','receipt','archival','unowned_job','scheduler_outage'])
def test_refresh_checks_full_snapshot_before_replacing_owner(resume_panel, monkeypatch, race):
    from contextlib import contextmanager, nullcontext
    import subprocess
    resume, root, graph = resume_panel
    old_failure(graph)
    state = dispatch._load(graph)[2]
    if race == 'submitting':
        state['tasks']['forward_10']['status'] = 'submitting'; atomic_json(graph/'state.json',state)
    elif race == 'staging':(graph/'tasks/new_task').mkdir()
    elif race == 'receipt':(graph/'requested.json').unlink()
    elif race == 'archival':atomic_json(graph/'.live_control_transaction.json', {'status':'applying'})
    events = []
    records = [dict(pid=123, starttime='1234')]
    monkeypatch.setattr(resume, 'inventory', lambda _:records)
    @contextmanager
    def freeze(*a, **kw):
        events.append('stopped')
        try:yield [(records[0], 12)]
        finally:events.append('resumed_if_alive')
    monkeypatch.setattr(resume.run_refresh, '_freeze', freeze)
    monkeypatch.setattr(resume.run_refresh, '_submission_guard', lambda _:nullcontext())
    monkeypatch.setattr(resume.run_refresh, '_pidfd_send_signal', lambda *a:events.append('kill_owner'))
    def jobs(_):
        if race == 'scheduler_outage':raise subprocess.CalledProcessError(1,'squeue')
        return [dict(job_id='99' if race=='unowned_job' else '77', directory=str(graph/'tasks/reverse_10'))]
    monkeypatch.setattr(resume.run_cleanup, '_jobs', jobs)
    monkeypatch.setattr(resume, 'launch', lambda _:events.append('launch_owner') or {'pid':456})
    before = (graph/'state.json').read_bytes()
    if race == 'none':
        assert resume.refresh(root)['pid'] == 456
        assert events == ['stopped','kill_owner','resumed_if_alive','launch_owner']
        assert len(list((root/'controller_handoffs').glob('*.json'))) == 1
        atomic_json(root/'native_owner.json', dict(pid=123, code=file_sha256(resume.__file__)))
        assert resume.refresh(root)['status'] == 'already_running'
        assert len(events) == 4
    else:
        with pytest.raises((RuntimeError, ValueError, subprocess.CalledProcessError)):
            resume.refresh(root)
        assert events == ['stopped','resumed_if_alive']
    assert (graph/'state.json').read_bytes() == before
