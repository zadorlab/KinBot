"""Validate the single follow-up through real dispatch, with simulated Slurm."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import json
import pytest
from kinbot import run_policy
from kinbot.anl import dispatch
from kinbot.anl.results import parse_result
from tests.test_anl_dispatch import _mock_execution
from tests.test_r8_validation_panel import site,spectrum

SPEC=Path('examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/review_20261119_closure/w3_complex_b_saddle.json')

@pytest.fixture
def staged(tmp_path,monkeypatch):
    spec=json.loads(SPEC.read_text())
    monkeypatch.setattr(dispatch,'assign_partitions',site)
    monkeypatch.setattr(dispatch,'render_site_setup',lambda _:'# simulated site\n')
    root=tmp_path/'candidate'
    dispatch.prepare(spec,root)
    return root

def test_single_task_preserves_theory_geometry_and_staging(staged,monkeypatch):
    root=staged;_,spec,state=dispatch._load(root)
    assert list(state['tasks'])==['saddle']
    assert spec['molecule']==json.loads(SPEC.read_text())['molecule']
    task=spec['tasks'][0]
    assert task['result_parser']['stationary_point']=='transition_state'
    text=(root/'tasks/saddle/saddle.com').read_text()
    assert 'B2PLYP/cc-pVTZ' in text and 'Opt=(TS,CalcFC,Tight,MaxCycle=150) Freq=Numer' in text
    assert not task.get('depends_on') and len(spec['tasks'])==1
    for _ in range(2):dispatch._verify_stage_files(root,task,state['tasks']['saddle'])

def test_deferred_admission_live_receipt_and_repeated_drive(staged,monkeypatch):
    root=staged;attempts=[];active=set()
    def submit(*args,**kwargs):
        attempts.append('admission')
        if len(attempts)==1:return None
        assert len(attempts)==2
        kwargs['before_submit']()
        active.add('800001')
        return SimpleNamespace(returncode=0,stdout='800001\n',stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j) in active)
    waits=[]
    def tick(_):
        status=dispatch._load(root)[2]['tasks']['saddle']['status']
        waits.append(status)
        assert len(waits)<=3
        if len(waits)==2:_mock_execution(root,'saddle')
        if len(waits)==3:active.clear()
    monkeypatch.setattr(dispatch.time,'sleep',tick)
    assert dispatch.main(['drive',str(root),'--interval','1'])==0
    assert waits==['staged','submitted','submitted']
    assert dispatch.main(['drive',str(root),'--interval','1'])==0
    assert dispatch._load(root)[2]['tasks']['saddle']['status']=='complete'
    assert len(attempts)==2

def test_spectrum_gate_does_not_accept_minimum_or_higher_order_saddle():
    request=deepcopy(json.loads(SPEC.read_text())['tasks'][0]['result_parser']);request.pop('file')
    freqs=[-30.,65.,175.,190.,345.,1008.,1029.,1090.,1506.,2714.,3099.,3192.]
    assert len(parse_result(spectrum(freqs),request)['frequencies_cm_inverse'])==12
    for bad in ([30.]+freqs[1:],[-30.,-10.]+freqs[2:]):
        with pytest.raises(ValueError):parse_result(spectrum(bad),request)
