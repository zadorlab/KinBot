"""Exercise controller replacement together with real cached QC handoffs."""
from copy import deepcopy
import json
from pathlib import Path
import pytest
from ase import Atoms
from kinbot import run_refresh, run_salvage
from kinbot.execution import atomic_json
from tests.test_run_refresh import live_tree, sampling_tree
from tests.test_qc_reuse import backend


def test_full_refresh_reconciles_retired_hir_without_touching_live_work(sampling_tree,monkeypatch):
    t=sampling_tree;worker=t.root/'ch3ss/123'
    (worker/'kinbot.db').unlink()
    monkeypatch.chdir(worker)
    q=backend(worker)
    job='parent_well_high';atoms=Atoms('HH',positions=[[0,0,0],[0,0,1]])
    old=q.db.write(atoms,name=job,data=dict(status='normal',energy=-1.))
    (worker/(job+'.exitcode')).write_text('42 0\n')
    marker=q.db.write(atoms,name=job,data=dict(status=0,reason='HIR restart',supersedes_row_id=old))
    (worker/(job+'.log.restart_'+str(marker))).write_text('retained original')
    # Retain unrelated running allocations and sampler processes through the
    # real freeze/archive/launch orchestration, including repeated invocation.
    before={p:p.read_bytes() for p in t.root.rglob('execution.json')}
    jobs=deepcopy(t.jobs)
    result=run_refresh.refresh_when_idle(t.root,apply=True,wait_idle=30)
    assert result['phase']=='launched' and t.jobs==jobs
    assert {p:p.read_bytes() for p in before}==before
    assert {p['pid'] for p in t.processes}=={700,701,702}
    assert run_refresh.refresh_when_idle(t.root,apply=True,wait_idle=30)==result
    restored=backend(worker)
    assert restored.check_qc(job)==0
    assert not (worker/(job+'.exitcode')).exists()
    assert restored.db.get(id=old).data.energy==-1.
    assert restored.check_qc(job)==0


def test_documented_single_pes_launch_can_refresh_without_campaign_receipts(live_tree,monkeypatch):
    t=live_tree;root=t.root/'ch3s'
    # Site fixture has a second independent PES: keep its controllers/jobs out
    # of this single-PES inventory, without changing them on disk.
    outside=deepcopy([p for p in t.processes if not Path(p['cwd']).is_relative_to(root)])
    t.jobs[:]=[j for j in t.jobs if Path(j['directory']).is_relative_to(root)]
    atomic_json(root/'.kinbot_accepted_input.json',{})
    (root/'pes_exploration.json').write_text('')
    assert not (root/'fresh_run.json').exists() and not (root/'prepared.json').exists()
    assert run_refresh.refresh(root)['phase']=='reviewed'
    assert not (root/'.kinbot_ordinary_run.json').exists()
    before=deepcopy(t.jobs)
    result=run_refresh.refresh(root,apply=True)
    assert result['phase']=='launched' and t.jobs==before
    assert outside==t.processes
    assert json.loads((root/'.kinbot_ordinary_run.json').read_text())['inputs']==[str(root/'input.json')]
    assert run_refresh.refresh(root,apply=True)==result
    assert len(t.launches)==1


def test_refresh_from_batch_does_not_detach_controllers(live_tree,monkeypatch):
    t=live_tree;monkeypatch.setenv('SLURM_JOB_ID','700')
    with pytest.raises(RuntimeError,match='original python -m kinbot.pes'):
        run_refresh.refresh(t.root,apply=True)
    assert not t.signals and not t.launches


def test_refresh_preserves_automatic_orbital_probes_and_waits_for_submission(live_tree,monkeypatch):
    from kinbot.execution import file_sha256
    t=live_tree
    folder=t.root/'ch3s/123/vrctst/auto_radical_pair/contract/panel/singlet'
    job=folder/'job';job.mkdir(parents=True)
    inp=job/'probe.inp';inp.write_text('synthetic native input')
    (job/'probe.sbatch').write_text('synthetic retained batch')
    entry=dict(name='probe',directory=str(job),input_sha256=file_sha256(inp),job_id='901',
               submission_status='submitting')
    path=atomic_json(folder/'reconnaissance.json',dict(jobs=[entry]))
    t.jobs.append(dict(job_id='901',directory=str(job),state='RUNNING'))
    with pytest.raises(run_refresh.RefreshTransientError,match='probe submission'):
        run_refresh.refresh(t.root,apply=True)
    assert not t.signals
    entry['submission_status']='submitted';atomic_json(path,dict(jobs=[entry]))
    before={p:p.read_bytes() for p in (path,inp,job/'probe.sbatch')}
    result=run_refresh.refresh(t.root,apply=True)
    assert result['phase']=='launched'
    assert {p:p.read_bytes() for p in before}==before
    assert any(j['job_id']=='901' for j in t.jobs)
    assert run_refresh.refresh(t.root,apply=True)==result
