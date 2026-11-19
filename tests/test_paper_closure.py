"""Necessity-limited scan panel: real staging/dispatch, simulated Slurm only."""
from copy import deepcopy
import fcntl
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot import run_policy
from kinbot.anl import dispatch
from kinbot.execution import atomic_json
from tests.test_anl_dispatch import _mock_execution
from tests.test_r8_validation_panel import load, site


@pytest.fixture
def panel(tmp_path, monkeypatch):
    module = load('paper_closure')
    monkeypatch.setattr(dispatch, 'assign_partitions', site)
    monkeypatch.setattr(dispatch, 'render_site_setup', lambda _: '# simulated site\n')
    return module, tmp_path/'paper_closure_20261119'


def test_exact_plan_stages_and_repeats_without_any_submission(panel, monkeypatch):
    module, root = panel
    monkeypatch.setattr(run_policy, 'submit', lambda *a, **k: pytest.fail('stage submits'))
    assert module.execute(root) == 0
    before = {p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    assert module.execute(root) == 0
    assert all(p.read_bytes() == data for p,data in before.items())
    assert set(module.graphs()) == {'w3_cs_out','complex_a_cs_in','complex_b_cs_in','w4_ch_out'}
    for name,spec in module.graphs().items():
        _,saved,state = dispatch._load(root/name)
        dispatch.validate_spec(saved)
        task=saved['tasks'][0]
        dispatch._verify_stage_files(root/name,task,state['tasks']['scan'])
        text=(root/name/'tasks/scan/scan.com').read_text()
        assert module.ROUTE in text and 'Opt=(ModRedundant,CalcFC,Tight,MaxCycle=150)' in text
        scan=spec['intent']['scan']; i,j=scan['bond']
        assert f'B {i+1} {j+1} S {scan["steps"]} {scan["step_angstrom"]:.12f}' in text
        assert saved['molecule'] == spec['molecule']
        assert saved['intent']['scientific_approval'] is False


def test_full_owner_submission_admission_race_and_restart(panel,monkeypatch):
    module,root=panel
    module.execute(root)
    names=list(module.graphs())
    first=root/names[0]
    _,_,state=dispatch._load(first)
    state['tasks']['scan'].update(status='submitted',job_id='100')
    atomic_json(first/'state.json',state)
    submitted=[];active={'100'};defer={names[1]}
    def submit(*args,**kwargs):
        name=Path(kwargs['cwd']).parents[1].name
        if name in defer:
            defer.remove(name)
            return None  # shared capacity unavailable on this pass
        assert name not in submitted and name != names[0]
        kwargs['before_submit']()
        job=str(101+len(submitted));submitted.append(name);active.add(job)
        return SimpleNamespace(returncode=0,stdout=job+'\n',stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    monkeypatch.setattr(dispatch,'_job_active',lambda job: str(job) in active)
    class Interrupted(Exception):pass
    monkeypatch.setattr(module.time,'sleep',lambda _:(_ for _ in ()).throw(Interrupted()))
    with pytest.raises(Interrupted):module.execute(root,submit=True)
    assert submitted == names[2:]
    assert dispatch._load(root/names[1])[2]['tasks']['scan']['status']=='staged'
    # Receipt arrives while its allocation remains live: preserve submitted state.
    _mock_execution(first,'scan')
    with pytest.raises(Interrupted):module.execute(root,submit=True)
    assert dispatch._load(first)[2]['tasks']['scan']['status']=='submitted'
    assert submitted == names[2:]+names[1:2]
    for name in names[1:]:_mock_execution(root/name,'scan')
    active.clear()
    assert module.execute(root,submit=True)==0
    assert module.execute(root,submit=True)==0
    assert len(submitted)==3  # the existing allocation was adopted
    assert all(dispatch._load(root/name)[2]['tasks']['scan']['status']=='complete' for name in names)
    progress=json.loads((root/'progress.json').read_text())
    assert progress['execution_finished'] and not progress['scientific_approval']
    # Simulated native receipts cannot manufacture a scientifically complete scan.
    report=json.loads((root/'scan_observations.json').read_text())
    assert not any(x['complete_constrained_scan'] for x in report['scans'].values())


def test_owner_lock_and_scheduler_outage_preserve_work(panel,monkeypatch):
    module,root=panel;module.execute(root)
    before={p:p.read_bytes() for p in root.glob('*/state.json')}
    with (root/'controller.lock').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert module.execute(root,submit=True)==0
    monkeypatch.setattr(dispatch,'advance',lambda *a,**k:(_ for _ in ()).throw(dispatch.SchedulerUnavailableError('offline')))
    class Interrupted(Exception):pass
    monkeypatch.setattr(module.time,'sleep',lambda _:(_ for _ in ()).throw(Interrupted()))
    with pytest.raises(Interrupted):module.execute(root,submit=True)
    progress=json.loads((root/'progress.json').read_text())
    assert len(progress['infrastructure'])==4 and not progress['execution_finished']
    assert all(p.read_bytes()==data for p,data in before.items())


def test_partial_stage_and_changed_inputs_stop_without_overwriting(panel,monkeypatch):
    module,root=panel
    name=next(iter(module.graphs())); folder=root/name
    folder.mkdir(parents=True);(folder/'workflow.json').write_text('interrupted')
    with pytest.raises(RuntimeError,match='incomplete staging'):module.execute(root)
    assert (folder/'workflow.json').read_text()=='interrupted'
    # Test mutations in a distinct complete staged panel.
    other=root.parent/'complete';module.execute(other)
    changed=deepcopy(module.load_inputs());changed['scans'][0]['scan']['steps']-=1
    monkeypatch.setattr(module,'load_inputs',lambda:changed)
    with pytest.raises(ValueError,match='inputs differ'):module.execute(other)


def test_staged_native_input_mutation_blocks_restart(panel):
    module,root=panel;module.execute(root)
    path=root/next(iter(module.graphs()))/'tasks/scan/scan.com'
    path.write_text(path.read_text().replace('B2PLYP','B3LYP'))
    with pytest.raises(RuntimeError,match='staged scan.com changed'):module.execute(root,submit=True)


def native_scan(case, *, skip_energy=False, wrong_atoms=False):
    from ase.data import atomic_numbers
    scan=case['scan'];i,j=scan['bond'];origin=np.array(case['molecule']['positions'])
    vec=origin[j]-origin[i];vec/=np.linalg.norm(vec)
    output=[]
    for step in range(scan['steps']+1):
        xyz=origin.copy();xyz[j]=xyz[i]+vec*(scan['initial_distance_angstrom']+step*scan['step_angstrom'])
        output += [' Input orientation:', ' -----',' Center Atomic Atomic Coordinates',' Number Number Type X Y Z',' -----']
        for n,(sym,pos) in enumerate(zip(case['molecule']['symbols'],xyz),1):
            z=8 if wrong_atoms else atomic_numbers[sym]
            output.append(f'{n} {z} 0 '+ ' '.join(f'{x:.8f}' for x in pos))
        output+=[' -----',' SCF Done: E(UB2PLYP) = -999.000000']
        if not skip_energy:
            e=-835.+.001*min(step,2)-.001*max(0,step-2)
            output += [f' E2(B2PLYP) = -0.1D+00 E(B2PLYP) = {e:.10f}D+00']
        output += [' Optimization completed.']
    return '\n'.join(output)+'\nNormal termination of Gaussian\n'


def test_scan_parser_uses_correlated_energy_and_requires_full_coordinate_sequence(panel):
    module,root=panel;root.mkdir()
    case=deepcopy(module.load_inputs()['scans'][0]);case['scan']['steps']=4
    path=root/'scan.log';path.write_text(native_scan(case))
    result=module.scan_points(path,case)
    assert result['complete_constrained_scan'] and len(result['points'])==5
    assert result['points'][2]['energy_hartree']==pytest.approx(-834.998)
    assert result['candidate_maximum_indices']==[2] and not result['scientific_approval']
    path.write_text(native_scan(case)+'Error termination\n')
    assert not module.scan_points(path,case)['complete_constrained_scan']
    path.write_text(native_scan(case).replace('Normal termination of Gaussian',''))
    assert not module.scan_points(path,case)['complete_constrained_scan']
    path.write_text(native_scan(case));case['scan']['step_angstrom']=.05
    assert not module.scan_points(path,case)['coordinate_targets_met']
    for flags in ({'skip_energy':True},{'wrong_atoms':True}):
        path.write_text(native_scan(case,**flags))
        with pytest.raises(ValueError,match='matching atoms or correlated energy'):
            module.scan_points(path,case)


def test_export_includes_closure_native_evidence_not_checkpoints(panel,tmp_path):
    module,root=panel;module.execute(root)
    task=root/next(iter(module.graphs()))/'tasks/scan'
    (task/'scan.log').write_text('evidence')
    (task/'scan.chk').write_text('binary placeholder')
    collector=load('collect')
    paths=collector.selected(root.parent)
    assert task/'scan.log' in paths and root/'plan.json' in paths
    assert task/'scan.chk' not in paths
    assert len(paths)==len(set(paths))
    plan=json.loads((root/'plan.json').read_text());plan['scans'][0]['name']='../outside'
    atomic_json(root/'plan.json',plan)
    with pytest.raises(ValueError,match='graph name'):collector.selected(root.parent)
