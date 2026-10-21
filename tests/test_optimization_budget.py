import json
from types import SimpleNamespace
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from kinbot.optimization_budget import run_optimizer
from kinbot.anl.optimization_recovery import _last_point, load_seed
from kinbot.execution import atomic_json, file_sha256


def test_extensions_retain_optimizer_and_do_not_repeat_after_success(tmp_path):
    calls=[]; answers=iter([False,False,True])
    def run(**kw): calls.append(kw);return next(answers)
    report=tmp_path/'budget.json'
    assert run_optimizer(SimpleNamespace(run=run),fmax=.001,steps=100,chunks=3,report=report)
    assert len(calls)==3 and all(c==dict(fmax=.001,steps=100) for c in calls)
    assert json.loads(report.read_text())['chunks_used']==3


def test_exhaustion_is_finite_and_calculator_errors_not_retried():
    calls=[]
    def run(**kw):calls.append(1);return False
    assert not run_optimizer(SimpleNamespace(run=run),fmax=.01,steps=100,chunks=3)
    assert len(calls)==3
    def bad(**kw):raise RuntimeError('native SCF failed')
    with pytest.raises(RuntimeError,match='native SCF'):
        run_optimizer(SimpleNamespace(run=bad),fmax=.01,steps=100,chunks=3)


def test_saved_point_requires_matching_evaluated_trajectory_and_seed_hash(tmp_path):
    first=Atoms('SH',positions=[[0,0,0],[0,0,1.4]])
    last=first.copy();last.positions[1,2]=1.36
    for atom,e in [(first,-10.),(last,-10.1)]:
        atom.calc=SinglePointCalculator(atom,energy=e,forces=np.zeros((2,3)))
    write(tmp_path/'optimization.traj',[first,last]);atomic_json(tmp_path/'task.json',{})
    (tmp_path/'optimization.log').write_text('BFGS: 100 12:00:00 -10.100000 .001\n')
    accepted=_last_point(tmp_path,{'optimizer':{'steps':100}},first)
    assert accepted.positions==pytest.approx(last.positions)
    write(tmp_path/'continuation.xyz',accepted)
    record={'optimizer_continuation':{'sha256':file_sha256(tmp_path/'continuation.xyz'),'chunks_used':1}}
    restored=first.copy();assert load_seed(tmp_path,record,restored)==1
    assert restored.positions==pytest.approx(last.positions)
    (tmp_path/'continuation.xyz').write_text('changed')
    with pytest.raises(RuntimeError,match='seed changed'):
        load_seed(tmp_path,record,first.copy())
    (tmp_path/'optimization.log').write_text('BFGS: 100 12:00:00 -10.200000 .001\n')
    with pytest.raises(RuntimeError,match='energies disagree'):
        _last_point(tmp_path,{'optimizer':{'steps':100}},first)


def test_sella_continuation_uses_logged_accepted_frame_not_hessian_probe(tmp_path):
    atoms=Atoms('SH',positions=[[0,0,0],[0,0,1.4]])
    points=[]
    for distance,energy in [(1.4,-10.),(1.36,-10.1),(1.37,-10.09)]:
        a=atoms.copy();a.positions[1,2]=distance
        a.calc=SinglePointCalculator(a,energy=energy,forces=np.zeros((2,3)));points.append(a)
    write(tmp_path/'optimization.traj',points)
    atomic_json(tmp_path/'task.json',{})
    (tmp_path/'optimization.log').write_text('Sella 100 12:00:00 -10.100000 .001 0.0 .01 1.0 1\n')
    accepted=_last_point(tmp_path,{'optimizer':{'steps':100}},atoms)
    assert accepted.positions[1,2]==pytest.approx(1.36)


def failed_graph(tmp_path,monkeypatch):
    from kinbot.anl import dispatch
    from tests.anl_fixture import dispatch_spec,molpro_task
    spec=dispatch_spec();spec['tasks']=[spec['tasks'][0],molpro_task('live','rhf\n',geometry_from='initial')]
    graph=dispatch.prepare(spec,tmp_path/'graph');state=dispatch._load(graph)[2]
    task=spec['tasks'][0];name=task['id'];directory=graph/'tasks'/name
    state['tasks'][name].update(status='unavailable',job_id='42')
    state['tasks']['live'].update(status='submitted',job_id='43')
    atomic_json(graph/'state.json',state)
    atomic_json(directory/'execution.json',dict(schema=1,task_id=name,status='failed',
        slurm_job_id='42',geometry_sha256=state['tasks'][name]['geometry_sha256'],
        error='ASE geometry optimization did not converge.'))
    from ase.io import read
    first=read(directory/'geometry.xyz');last=first.copy();last.positions[0,0]+=.01
    for a,e in [(first,-10.),(last,-10.1)]:
        a.calc=SinglePointCalculator(a,energy=e,forces=np.zeros((len(a),3)))
    write(directory/'optimization.traj',[first,last])
    (directory/'optimization.log').write_text('Sella 120 12:00:00 -10.100000 .001 0 .01 1 1\n')
    monkeypatch.setattr(dispatch,'_job_active',lambda j:j=='43')
    return graph,name,last


def test_graph_continuation_is_idempotent_and_preserves_live_sibling(tmp_path,monkeypatch):
    from kinbot.anl import dispatch
    graph,name,last=failed_graph(tmp_path,monkeypatch)
    original=(graph/'workflow.json').read_bytes()
    sibling={p:file_sha256(p) for p in (graph/'tasks/live').iterdir() if p.is_file()}
    state=dispatch.advance(graph)
    assert state['tasks'][name]['status']=='staged'
    assert state['tasks']['live']['job_id']=='43'
    assert (graph/'workflow.json').read_bytes()==original
    assert sibling=={p:file_sha256(p) for p in sibling}
    record=json.loads((graph/'tasks'/name/'task.json').read_text())
    from ase.io import read
    atoms=read(graph/'tasks'/name/'geometry.xyz')
    assert load_seed(graph/'tasks'/name,record,atoms)==1
    assert atoms.positions==pytest.approx(last.positions)
    assert dispatch.advance(graph)['optimizer_continuations']==state['optimizer_continuations']


def test_graph_continuation_recovers_after_archive_write_interruption(tmp_path,monkeypatch):
    from kinbot import run_recover
    from kinbot.anl import dispatch
    graph,name,_=failed_graph(tmp_path,monkeypatch)
    original=run_recover.atomic_json
    def crash(path,data):
        if Path(path)==graph/'state.json':raise OSError('interrupted')
        return original(path,data)
    from pathlib import Path
    monkeypatch.setattr(run_recover,'atomic_json',crash)
    with pytest.raises(OSError,match='interrupted'):dispatch.advance(graph)
    monkeypatch.setattr(run_recover,'atomic_json',original)
    state=dispatch.advance(graph)
    assert state['tasks'][name]['status']=='staged'
    assert state['optimizer_continuations'][name]['chunks_used']==1
