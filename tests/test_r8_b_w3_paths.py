"""Reuse a real native Hessian and exercise the entire paths-only owner."""
from copy import deepcopy
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import shutil

from ase.io import read
import pytest
from kinbot import saddle_recheck as checks
from kinbot.anl.results import parse_gaussian_harmonic
from kinbot.execution import atomic_json,file_sha256
from kinbot.qc import QuantumChemistry

SCRIPT=Path('examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/review_20261120_resolution/stage_b_w3_paths.py')
FIXTURE=Path('tests/reference/mirror_even_transition.log')

@pytest.fixture
def imported(tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('b_w3_paths',SCRIPT)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=tmp_path/'source';task=source/'tasks/saddle';task.mkdir(parents=True)
    shutil.copy2(FIXTURE,task/'saddle.log')
    (task/'saddle.com').write_text('original input retained\n')
    atomic_json(task/'execution.json',dict(status='executed',slurm_job_id='123'))
    atomic_json(source/'workflow.json',dict(molecule=dict(charge=0,multiplicity=2)))
    spectrum=parse_gaussian_harmonic(FIXTURE.read_text(),method='B2PLYP',basis='cc-pVTZ',
        dispersion='GD3BJ',natoms=6,stationary_point='transition_state')
    monkeypatch.setattr(module,'SOURCE_HASH',file_sha256(task/'saddle.log'))
    monkeypatch.setattr(module,'EXPECTED_FREQUENCIES',spectrum['frequencies_cm_inverse'])
    return module,source,tmp_path/'paths'

def test_native_import_and_complete_paths_owner(imported,monkeypatch):
    module,source,root=imported
    source_before={p:p.read_bytes() for p in source.rglob('*') if p.is_file()}
    module.stage(source,root);module.stage(source,root)
    assert {p:p.read_bytes() for p in source.rglob('*') if p.is_file()}==source_before
    campaign=json.loads((root/'campaign.json').read_text())
    case=campaign['cases'][0]
    assert len(campaign['cases'])==1 and case['mode']=='reuse'
    assert case['native']['profile']['optimizer']=='gaussian'
    assert campaign['parameters']['l1_profile']['optimizer']=='sella'
    jobs=[];polls={};atoms=read(source/'tasks/saddle/saddle.log',format='gaussian-out')
    def submit(qc,job,nproc,**kwargs):
        assert job in ('saddle_IRC_F','saddle_IRC_R') and job not in jobs
        jobs.append(job);qc.job_ids[job]=str(100+len(jobs));return 1
    def status(qc,job,**kwargs):
        if qc.db.count(name=job):return 'normal'
        if job not in jobs:return 0
        polls[job]=polls.get(job,0)+1
        if polls[job]<3:return 'running'
        qc.db.write(atoms,name=job,data=dict(status='normal',energy=-1.))
        qc.db.write(atoms,name=job+'_prod',data=dict(status='normal',energy=-1.))
        atomic_json(job+'.irc_path.json',dict(converged=True,moved=True))
        return 'normal'
    monkeypatch.setattr(QuantumChemistry,'submit_qc',submit)
    monkeypatch.setattr(QuantumChemistry,'check_qc',status)
    monkeypatch.setattr(QuantumChemistry,'qc_freq',lambda *a,**kw:pytest.fail('Unrequested frequency'))
    monkeypatch.setattr(checks.time,'sleep',lambda _:None)
    checks.run(root,paths_only=True)
    assert jobs==['saddle_IRC_F','saddle_IRC_R']
    assert json.loads((root/'results.json').read_text())['b_w3']['outcome']=='irc_paths_available'
    module.stage(source,root);checks.run(root,paths_only=True)
    assert len(jobs)==2
    assert {p:p.read_bytes() for p in source.rglob('*') if p.is_file()}==source_before

def test_source_hash_and_execution_are_enforced(imported):
    module,source,root=imported
    with pytest.raises(ValueError,match='differs'):
        module.build_manifest(source,expected_hash='wrong')
    atomic_json(source/'tasks/saddle/execution.json',dict(status='failed'))
    with pytest.raises(ValueError,match='successfully'):module.stage(source,root)
    assert not root.exists()

@pytest.mark.parametrize('change',[{'method':'B3LYP'}, {'basis':'cc-pVDZ'},
    {'calculator_kwargs':{'EmpiricalDispersion':'GD3'}}])
def test_optimizer_reuse_does_not_weaken_electronic_theory_check(imported,monkeypatch,change):
    module,source,root=imported;module.stage(source,root)
    campaign=json.loads((root/'campaign.json').read_text());case=deepcopy(campaign['cases'][0])
    case['native']['profile'].update(change)
    with checks.in_directory(root/'b_w3'):
        parameters=checks.Parameters('input.json',show_warnings=False).par
        qc=QuantumChemistry(parameters)
        with pytest.raises(ValueError,match='theory differs'):checks.adopt(case,qc,parameters)
        assert not qc.db.count(name='saddle')
