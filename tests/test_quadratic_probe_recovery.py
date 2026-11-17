"""Bounded third-wave electronic audit through the real panel owner."""
from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from kinbot import run_policy
from kinbot.anl import dispatch, path_review, probe_recovery, protocol
from kinbot.execution import atomic_json, file_sha256
from tests.test_gaussian_checkpoint_inputs import graph, resume_panel
from tests.test_native_path_review import failure, make_plan
from tests.test_native_probe_recovery import audit, scf_native

pytestmark = pytest.mark.parametrize('audit', ['scf'], indirect=True)


def quadratic_plan(resume, root, graph, audit, monkeypatch, base_factory=None):
    monkeypatch.setattr(dispatch, '_job_active', lambda j: str(j) == '77')
    if base_factory:
        base = base_factory(resume, root, graph, monkeypatch)
    else:
        failure(graph)
        base = make_plan(resume, root, monkeypatch)
    source, _ = audit
    source.rename(root/'audit_old'); source = root/'audit_old'
    spec = dispatch._load(source)[1]; atomic_json(source/'requested.json', spec)
    proof = path_review.source(source, 'gradient', ['gradient.out'])[3]
    base['curvature']['audit_old'] = dict(spec=spec, proof=proof)
    base['sha256'] = path_review.digest({k:v for k,v in base.items() if k != 'sha256'})
    atomic_json(root/resume.reconciliation.DECISION, base)
    second = resume.reconciliation.prepare_repairs(root)
    resume.reconciliation.apply(root, second)
    # The source spec has name=audit, while its retained directory is audit_old.
    # prepare_repairs keys by source directory, so use the declared decision key.
    target = root/'audit_old_electronic_v2'
    _, _, state = dispatch._load(target)
    entry = state['tasks']['electronic']; entry.update(status='failed', job_id='210')
    d = target/'tasks/electronic'
    native = 'Version 2024.1\n'+scf_native()+'Final alpha occupancy: 5\nFinal beta  occupancy: 4\n'
    (d/'electronic.out').write_text(native)
    atomic_json(d/'execution.json', dict(schema=1, task_id='electronic', slurm_job_id='210',
        geometry_sha256=entry['geometry_sha256'], status='failed', failure_phase='execution'))
    for ident in probe_recovery.DERIVATIVES: state['tasks'].setdefault(ident, {})['status'] = 'blocked'
    atomic_json(target/'state.json', state)
    prior = {p:p.read_bytes() for p in (root/resume.reconciliation.DECISION, root/resume.reconciliation.REPAIRS)}
    plan = resume.reconciliation.prepare_quadratic(root)
    assert prior == {p:p.read_bytes() for p in prior}
    assert resume.reconciliation.prepare_quadratic(root) == plan
    return plan, target


@pytest.mark.parametrize('native_success', [True, False])
def test_quadratic_owner_worker_preserves_sources_and_never_releases_derivatives(audit, resume_panel, monkeypatch, tmp_path, native_success):
    resume, root, graph = resume_panel
    plan, source = quadratic_plan(resume, root, graph, audit, monkeypatch)
    name = next(n for n in plan['curvature'] if n.endswith('_quadratic_v3'))
    spec = plan['curvature'][name]['spec']; task = spec['tasks'][0]
    assert len(spec['tasks']) == 1 and task['id'] == 'electronic'
    assert 'rhf,so,maxit=60,energy=1.d-12;wf,9,1,1;occ,5;closed,4;orbital,2100.2' in task['input_template']
    assert spec['molecule'] == dispatch._load(source)[1]['molecule']
    assert 'thrden=1.d-11,thrvar=1.d-11' in task['input_template']
    assert 'shift,' not in task['input_template']
    before = {p:file_sha256(p) for p in root.glob('*/tasks/*/*') if p.is_file()}
    monkeypatch.setattr(protocol, 'preflight_ready', lambda _:True)
    submitted = []
    def submit(command, **kwargs):
        kwargs['before_submit'](); submitted.append(Path(kwargs['cwd']))
        return SimpleNamespace(returncode=0, stdout=str(300+len(submitted)), stderr='')
    monkeypatch.setattr(run_policy, 'submit', submit)
    class Tick(Exception): pass
    monkeypatch.setattr(resume.time, 'sleep', lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick): resume.run(root)
    assert [p.name for p in submitted] == ['forward_10_dvv_v1', 'electronic']
    assert dispatch._load(graph)[2]['tasks']['reverse_10']['job_id'] == '77'
    assert before == {p:file_sha256(p) for p in before}
    executable = tmp_path/'molpro_so.py'
    executable.write_text(f'#!{sys.executable}\nfrom pathlib import Path\nimport sys\n'
        'p=Path(sys.argv[-1]);s=p.read_text()\n'
        'assert "rhf,so,maxit=60,energy=1.d-12;wf,9,1,1;occ,5;closed,4" in s\n'
        'assert "forces," not in s and "frequencies," not in s\n'
        'p.with_suffix(".out").write_text("Molpro calculation terminated\\n")\n'
        'Path("electronic.molden").write_text("[Molden Format]\\n")\n')
    executable.chmod(0o755)
    if not native_success:
        executable.write_text(executable.read_text().replace('Molpro calculation terminated', 'No convergence in MULTI'))
    monkeypatch.setattr(dispatch, '_external_command', lambda t,*a,**kw:[str(executable), t['input_name']])
    monkeypatch.setenv('SLURM_JOB_ID', '302')
    target = root/name
    if native_success: dispatch.run_task(target/'tasks/electronic/task.json')
    else:
        with pytest.raises(RuntimeError): dispatch.run_task(target/'tasks/electronic/task.json')
    monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','301','302'))
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(target)[2]['tasks']['electronic']['status'] == 'submitted'
    monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','301'))
    with pytest.raises(Tick): resume.run(root)
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(target)[2]['tasks']['electronic']['status'] == ('complete' if native_success else 'failed')
    assert len(submitted) == 2 and before == {p:file_sha256(p) for p in before}
    assert resume.reconciliation.prepare_quadratic(root) == plan
    assert not json.loads((root/'progress.json').read_text())['scientific_approval']
    with pytest.raises(ValueError): probe_recovery.quadratic_rohf_probe(target)


def test_third_wave_decision_rejects_wrong_ancestry(audit, resume_panel, monkeypatch):
    resume, root, graph = resume_panel
    quadratic_plan(resume, root, graph, audit, monkeypatch)
    p=root/resume.reconciliation.ELECTRONIC; record=json.loads(p.read_text())
    record['parent_sha256']='not_the_parent'
    record['sha256']=path_review.digest({k:v for k,v in record.items() if k!='sha256'})
    atomic_json(p,record)
    with pytest.raises(ValueError, match='parent decision'): resume.reconciliation.load(root)


@pytest.mark.parametrize('defect', ['active','wrong_job','version','occupations','contracting','started_derivative','changed_source','geometry','second_order_again'])
def test_quadratic_rejects_inapplicable_or_changed_evidence(audit, resume_panel, monkeypatch, defect):
    resume, root, graph = resume_panel
    plan, source = quadratic_plan(resume, root, graph, audit, monkeypatch)
    if defect == 'active': monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','210'))
    elif defect == 'wrong_job':
        p = source/'tasks/electronic/execution.json'; e=json.loads(p.read_text());e['slurm_job_id']='211';atomic_json(p,e)
    elif defect in ('version','occupations','contracting'):
        p=source/'tasks/electronic/electronic.out'; text=p.read_text()
        if defect=='version': text=text.replace('2024.1','2019.1')
        elif defect=='occupations': text=text.replace('occupancy: 5','occupancy: 6')
        else: text=text.replace('0.39D-04','0.39D-10')
        p.write_text(text)
    elif defect=='started_derivative':
        p=source/'state.json';s=json.loads(p.read_text());s['tasks']['gradient']['job_id']='211';atomic_json(p,s)
    elif defect=='changed_source': (root/'audit_old/tasks/gradient/gradient.out').write_text('changed')
    elif defect=='geometry': (source/'tasks/electronic/geometry.xyz').write_text('3\nchanged\nO 0 0 0\nH 1.1 0 0\nH 0 1 0\n')
    else:
        trial = next(x['spec'] for n,x in plan['curvature'].items() if n.endswith('_quadratic_v3'))
        resume.reconciliation.stage_extra(root, 'already_so', trial); source=root/'already_so'
    with pytest.raises((ValueError, RuntimeError)): probe_recovery.quadratic_rohf_probe(source)


@pytest.mark.parametrize('race', ['none','submitting','staging','archival','unowned_job','scheduler_outage'])
def test_third_wave_refresh_preserves_live_work_and_transient_races(audit, resume_panel, monkeypatch, race):
    from tests import test_native_path_review as original
    def third_wave(resume, root, patcher):
        return quadratic_plan(resume, root, resume_panel[2], audit, patcher)[0]
    monkeypatch.setattr(original, 'make_plan', third_wave)
    original.test_decision_refresh_checks_before_retiring_controller(resume_panel, monkeypatch, race)
