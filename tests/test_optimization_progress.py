"""Accepted-point continuation survives abrupt exits without resetting budgets."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from kinbot.anl.optimization_progress import Progress, ownership, restore_attempt_start
from kinbot.execution import file_sha256


def atoms(distance=1.4):
    a=Atoms('SH',positions=[[0,0,0],[0,0,distance]])
    a.calc=SinglePointCalculator(a,energy=-10.,forces=np.zeros((2,3)))
    return a


def optimizer(name='BFGS',nsteps=0):
    obj=type(name,(),{})();obj.nsteps=nsteps
    obj.attach=lambda callback,interval:setattr(obj,'save',callback)
    return obj


def test_sella_saves_logged_accepted_frame_not_last_hessian_probe(tmp_path):
    a=atoms();record={'task':{'kind':'ase_optimize'}}
    progress=Progress(tmp_path,record,a,used_chunks=0,chunk_steps=100,max_chunks=3)
    accepted=atoms(1.35);probe=atoms(1.36)
    write(tmp_path/'optimization.traj',[a,accepted,probe])
    (tmp_path/'optimization.log').write_text('Sella 7 12:00:00 -10.000000 .001 0 .01 1 1\n')
    opt=optimizer('Sella',7);progress.attach(opt,probe);opt.save()
    fresh=atoms();restart=Progress(tmp_path,record,fresh,used_chunks=0,chunk_steps=100,max_chunks=3)
    assert fresh.positions==pytest.approx(accepted.positions)
    assert restart.steps==7 and restart.remaining()==293
    old=(tmp_path/'optimization.traj').read_bytes()
    archive=restart.archive_previous_attempt('geom');restart.receipt(archive)
    assert (tmp_path/archive/'optimization.traj').read_bytes()==old
    expected=atoms();restore_attempt_start(tmp_path,record,expected)
    assert expected.positions==pytest.approx(accepted.positions)
    assert not (tmp_path/'optimization.traj').exists()


def test_task_change_cannot_reuse_and_corruption_cannot_overwrite_checkpoint(tmp_path):
    a=atoms();record={'task':{'profile':{'method':'CCSD(T)'}}}
    p=Progress(tmp_path,record,a,used_chunks=0,chunk_steps=100,max_chunks=3)
    opt=optimizer(nsteps=2);p.attach(opt,a);opt.save()
    altered=deepcopy(record);altered['task']['profile']['method']='other'
    assert Progress(tmp_path,altered,a.copy(),used_chunks=0,chunk_steps=100,max_chunks=3).loaded is None
    pointer=json.loads((p.root/'latest.json').read_text());path=p.root/pointer['generation']
    path.write_text('corrupt');before=path.read_bytes()
    with pytest.raises(RuntimeError,match='checkpoint changed'):
        Progress(tmp_path,record,a,used_chunks=0,chunk_steps=100,max_chunks=3)
    assert path.read_bytes()==before


def test_pointer_publication_failure_keeps_prior_generation(tmp_path,monkeypatch):
    import kinbot.anl.optimization_progress as module
    a=atoms();p=Progress(tmp_path,{},a,used_chunks=0,chunk_steps=100,max_chunks=3)
    opt=optimizer(nsteps=1);p.attach(opt,a);opt.save()
    pointer=(p.root/'latest.json').read_bytes();original=module.publish
    def fail(path,data):
        if path.name=='latest.json':raise OSError('temporary persistent filesystem outage')
        original(path,data)
    monkeypatch.setattr(module,'publish',fail);opt.nsteps=2
    with pytest.raises(OSError):opt.save()
    assert (p.root/'latest.json').read_bytes()==pointer
    assert Progress(tmp_path,{},a.copy(),used_chunks=0,chunk_steps=100,max_chunks=3).steps==1
    assert len(list(p.root.glob('point-*.json')))==2


def test_step_budget_not_reset_by_repeated_process_restarts(tmp_path):
    from kinbot.optimization_budget import run_optimizer
    a=atoms();record={};p=Progress(tmp_path,record,a,used_chunks=0,chunk_steps=100,max_chunks=3)
    opt=optimizer(nsteps=180);p.attach(opt,a);opt.save()
    p=Progress(tmp_path,record,a.copy(),used_chunks=0,chunk_steps=100,max_chunks=3)
    calls=[];opt=optimizer();p.attach(opt,a)
    def run(fmax,steps):
        calls.append(steps);opt.nsteps+=steps;opt.save();return False
    opt.run=run
    assert not run_optimizer(opt,fmax=.03,steps=100,chunks=3,report=tmp_path/'budget.json',
                             remaining_steps=p.remaining(),used_steps=p.steps)
    assert calls==[100,20] and p.steps==300
    assert json.loads((tmp_path/'budget.json').read_text())['chunks_used']==3
    assert Progress(tmp_path,record,a.copy(),used_chunks=0,chunk_steps=100,max_chunks=3).remaining()==0


def test_concurrent_optimizer_refused_and_requeue_does_not_rollback_generations(tmp_path,monkeypatch):
    from kinbot.slurm_retry import _artifacts
    p=Progress(tmp_path,{},atoms(),used_chunks=0,chunk_steps=100,max_chunks=3)
    opt=optimizer();p.attach(opt,atoms());opt.save()
    with ownership(tmp_path):
        with pytest.raises(RuntimeError,match='Another process'):
            with ownership(tmp_path):pass
    paths=_artifacts([], [tmp_path], tmp_path/'.kinbot_slurm_io_retries')
    assert not any('.kinbot_checkpoints' in x.parts for x in paths)


def test_dispatcher_resumes_real_bfgs_after_abrupt_calculator_exit(tmp_path,monkeypatch):
    from kinbot.anl import dispatch
    from kinbot.ase_modules.calculators import factory
    record=dict(molecule=dict(charge=0,multiplicity=1),task=dict(id='geom',kind='ase_optimize',
        profile={},optimizer=dict(steps=20,step_retries=0,fmax=.001)))
    original=Atoms('H2',positions=[[0,0,0],[0,0,2.]])
    write(tmp_path/'geometry.xyz',original)
    starts=[];crash=[True]
    class Harmonic(Calculator):
        implemented_properties=['energy','forces']
        def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
            super().calculate(atoms,properties,system_changes)
            starts.append(atoms.positions.copy())
            if crash[0] and len(starts)==3:raise KeyboardInterrupt('simulated node death')
            delta=atoms.positions-original.positions/2
            self.results=dict(energy=float((delta**2).sum()/2),forces=-delta)
    monkeypatch.setattr(factory,'build_calculator',lambda *a,**k:Harmonic())
    monkeypatch.setattr(dispatch,'_runtime_profile',lambda t:SimpleNamespace(optimizer='bfgs'))
    with pytest.raises(KeyboardInterrupt):dispatch._run_ase_optimize(tmp_path,record)
    points=list((tmp_path/'.kinbot_checkpoints').rglob('point-*.json'));assert len(points)>=2
    old=(tmp_path/'optimization.traj').read_bytes();crash[0]=False;starts.clear()
    result=dispatch._run_ase_optimize(tmp_path,record)
    receipt=json.loads((tmp_path/'optimization_resume.json').read_text())
    assert starts[0]==pytest.approx(np.array(receipt['checkpoint']['positions']))
    assert not np.allclose(starts[0],original.positions)
    assert (tmp_path/receipt['archive']/'optimization.traj').read_bytes()==old
    assert result['optimization_restart']['checkpoint'] is not None
    assert result['optimization_restart']['steps_used']<=20
    assert np.max(np.abs(read(tmp_path/'final.xyz').positions-original.positions/2))<.001


def test_real_sella_recovers_accepted_point_after_hessian_probe_crash(tmp_path,monkeypatch):
    from kinbot.anl import dispatch
    from kinbot.ase_modules.calculators import factory
    record=dict(molecule=dict(charge=0,multiplicity=2),task=dict(id='geom',kind='ase_optimize',
        profile={},optimizer=dict(steps=40,step_retries=0,fmax=.001)))
    original=Atoms('H3',positions=[[0,0,0],[0,0,1.8],[1.5,0,.4]])
    write(tmp_path/'geometry.xyz',original)
    crash=[True];calls=[]
    class PairSprings(Calculator):
        implemented_properties=['energy','forces']
        def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
            super().calculate(atoms,properties,system_changes);calls.append(atoms.positions.copy())
            for pointer in (tmp_path/'.kinbot_checkpoints').glob('optimization/*/latest.json'):
                p=json.loads(pointer.read_text());point=json.loads((pointer.parent/p['generation']).read_text())
                if crash[0] and point['steps_used']>=1:raise KeyboardInterrupt('node lost during Hessian probe')
            e=0.;forces=np.zeros_like(atoms.positions)
            for i,j in ((0,1),(0,2),(1,2)):
                delta=atoms.positions[i]-atoms.positions[j];d=np.linalg.norm(delta);f=-(d-1.1)*delta/d
                e+=(d-1.1)**2/2;forces[i]+=f;forces[j]-=f
            self.results=dict(energy=e,forces=forces)
    monkeypatch.setattr(factory,'build_calculator',lambda *a,**k:PairSprings())
    monkeypatch.setattr(dispatch,'_runtime_profile',lambda t:SimpleNamespace(optimizer='sella'))
    with pytest.raises(KeyboardInterrupt):dispatch._run_ase_optimize(tmp_path,record)
    crash[0]=False;calls.clear()
    result=dispatch._run_ase_optimize(tmp_path,record)
    saved=result['optimization_restart']['checkpoint']
    assert saved['steps_used']>=1
    assert calls[0]==pytest.approx(np.array(saved['positions']))
    assert not np.allclose(calls[0],original.positions)
    final=read(tmp_path/'final.xyz')
    assert [final.get_distance(i,j) for i,j in ((0,1),(0,2),(1,2))]==pytest.approx([1.1]*3,abs=.002)
