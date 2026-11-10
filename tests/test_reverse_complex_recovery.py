import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import Mock

import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

from kinbot import reverse_complex_recovery, qc_optimizer_recovery, qc_timeout, run_policy
from kinbot.execution import atomic_json, file_sha256
from tests.test_timeout_recovery import pool, ordinary, allocation


def stalled_complex(ordinary):
    qc,calls=ordinary
    qc.qc='gauss'
    source="""from ase import Atoms
mol=Atoms(symbols=['C','H','H','H','H'],positions=[[0,0,0],[0,0,1],[0,1,0],[1,0,0],[0,0,-1]])
kwargs={'method':'B2PLYP','basis':'cc-pVTZ','mult':2,'charge':0}
order=0
sella_kwargs={'internal':True,'gamma':0.0}
opt=Sella(mol,order=order,trajectory='well.traj',**sella_kwargs)
converged=run_optimizer(opt,fmax=0.0005,steps=250)
if converged: validate_frequency(mol)
"""
    Path('well.py').write_text(source)
    token=hashlib.sha256(b'well').hexdigest()[:20]
    guard=Path('.kinbot_slurm_io_retries')/token/'42/execution.json'
    value=json.loads(guard.read_text())
    (guard.parent/'initial/well.py').write_text(source)
    value['initial_artifacts']['well.py']=file_sha256('well.py');atomic_json(guard,value)
    frames=[]
    for e in [-10.,-10.1]:
        a=Atoms('CH4',positions=[[0,0,0],[0,0,1],[0,1,0],[1,0,0],[0,0,-1]])
        a.calc=SinglePointCalculator(a,energy=e,forces=[[.01,0,0]]*5);frames.append(a)
    write('well.traj',frames)
    Path('well_sella.log').write_text('Sella 750 12:00:00 -10.100000 .01 0 .0001 1 1\nWorker failure: Optimization did not converge within 750 steps.\n')
    atomic_json(Path('well.optimization_budget.json'),dict(chunks_used=3,status='exhausted'))
    atomic_json(Path('.kinbot_timeout_allocations/42.json'),allocation(Path.cwd(),status='RUNNING'))
    qc._reverse_complex_recovery={'well':dict(transition_state='ts_high',row_id=2)}
    return qc,calls


def test_cartesian_recovery_is_bounded_retains_science_and_adopts_once(ordinary,monkeypatch):
    qc,calls=stalled_complex(ordinary)
    monkeypatch.setattr('kinbot.timeout_retry.evidence',lambda *a,**k:None)
    assert qc_timeout.recover(qc,'well')=='running' and len(calls)==1
    record=json.loads(qc_timeout._marker('well').read_text())
    assert record['attempts'][-1]['source']=='optimizer_cartesian'
    archive=Path(record['attempts'][-1]['archive'])
    text=(archive/'retry_worker.py').read_text()
    assert "sella_kwargs['internal'] = False" in text
    assert 'fmax=0.0005' in text and 'steps=250' in text and 'chunks=3' in text
    assert 'order = 0' in text and 'validate_frequency(mol)' in text
    assert "'method': 'B2PLYP'" in text and "'basis': 'cc-pVTZ'" in text
    assert (archive/'files/well.optimization_budget.json').is_file()
    assert Path('other.chk').read_text()=='sibling'
    qc.job_ids['well']='42';qc_timeout.resume(qc,'well')
    assert qc.job_ids['well']=='45' and len(calls)==1
    # Even a newly failed replacement must not run Cartesian again.
    for name in ('well.traj','well_sella.log','well.optimization_budget.json'):
        Path(name).write_bytes((archive/'files'/name).read_bytes())
    assert qc_optimizer_recovery.evidence(qc,'well','45',cartesian_complex=True) is None
    assert len(calls)==1


@pytest.mark.parametrize('cause',['no_ts','saddle','cartesian','frequency','scf','constraints','bad_trajectory'])
def test_nonapplicable_failures_are_not_retried(ordinary,cause):
    qc,calls=stalled_complex(ordinary)
    p=Path('well.py');source=p.read_text()
    if cause=='no_ts':qc._reverse_complex_recovery={}
    elif cause=='saddle':p.write_text(source.replace('order=0','order=1'))
    elif cause=='cartesian':p.write_text(source.replace("'internal':True","'internal':False"))
    elif cause=='constraints':p.write_text(source.replace('order=order,','order=order,constraints=constraints,'))
    elif cause=='bad_trajectory':Path('well_sella.log').write_text(Path('well_sella.log').read_text().replace('-10.100000','-12.100000'))
    else:Path('well_sella.log').write_text('Worker failure: '+('Frequency validation failed' if cause=='frequency' else 'SCF failed')+'\n')
    assert qc_optimizer_recovery.evidence(qc,'well','42',cartesian_complex=True) is None
    assert not calls


def test_staging_race_preserves_single_retry_and_result_adoption(ordinary,monkeypatch):
    qc,calls=stalled_complex(ordinary)
    monkeypatch.setattr('kinbot.timeout_retry.evidence',lambda *a,**k:None)
    submit=run_policy.submit
    monkeypatch.setattr(run_policy,'submit',lambda *a,**k:None)
    assert qc_timeout.recover(qc,'well')=='running' and not calls
    monkeypatch.setattr(run_policy,'submit',submit)
    qc_timeout.resume(qc,'well');qc_timeout.resume(qc,'well')
    assert len(calls)==1
    qc.db.write(Atoms('H'),name='ts_high',data={'status':'normal'})
    qc.check_qc=lambda job:'normal'
    opt=NS(shigh=-999,just_high=True,species=NS(wellorts=0),par={'high_level':1},qc=qc,log_name=lambda n:'well')
    opt.do_optimization=Mock()
    ts=NS(shigh=1,species=NS(wellorts=1),selected_job='ts_high')
    reverse_complex_recovery.drive(opt,ts)
    assert opt.shigh==0
    opt.do_optimization.assert_called_once_with()
    assert qc._reverse_complex_recovery['well']['transition_state']=='ts_high'


@pytest.mark.parametrize('ts_status,row_status',[(0,'normal'),(-999,'normal'),(1,'error'),(1,None)])
def test_unaccepted_ts_does_not_enable_complex_recovery(ordinary,ts_status,row_status):
    qc,calls=stalled_complex(ordinary)
    qc._reverse_complex_recovery={}
    if row_status is not None:
        qc.db.write(Atoms('H'),name='ts_high',data={'status':row_status})
    qc.check_qc=Mock()
    opt=NS(shigh=-999,just_high=True,species=NS(wellorts=0),par={'high_level':1},qc=qc,
           log_name=lambda n:'well',do_optimization=Mock())
    ts=NS(shigh=ts_status,species=NS(wellorts=1),selected_job='ts_high')
    reverse_complex_recovery.drive(opt,ts)
    assert not qc._reverse_complex_recovery and not calls
    qc.check_qc.assert_not_called()
    opt.do_optimization.assert_not_called()


@pytest.mark.parametrize('active',[True,None])
def test_live_or_unconfirmed_allocation_is_preserved(ordinary,monkeypatch,active):
    qc,calls=stalled_complex(ordinary)
    monkeypatch.setattr('kinbot.timeout_retry.evidence',lambda *a,**k:None)
    monkeypatch.setattr('kinbot.qc._slurm_job_active',lambda j:active)
    before=file_sha256('well.traj')
    assert qc_timeout.recover(qc,'well')=='running'
    assert file_sha256('well.traj')==before and not calls
    assert not qc_timeout._marker('well').exists()


def test_explicit_retry_takes_precedence_over_automatic_strategy(ordinary,monkeypatch):
    from kinbot import live_control, qc_live
    qc,calls=stalled_complex(ordinary)
    live_control.configure(Path.cwd())
    live_control.observe_qc(Path.cwd()/'well','42','error')
    live_control.request(Path.cwd(),'well',kind='qc',reason='operator Cartesian choice',
        fallback={'schema':1,'worker_parameters':{'sella_kwargs':{'internal':False}}})
    monkeypatch.setattr('kinbot.anl.live_control.slurm_job_record',lambda j:allocation(Path.cwd(),j))
    assert qc_timeout.recover(qc,'well') is None and not calls
    assert qc_live.process(qc,'well','error')=='running' and len(calls)==1
    record=json.loads(qc_timeout._marker('well').read_text())
    assert record['attempts'][-1]['source']=='operator'


@pytest.mark.parametrize('condition,eligible',[
    ("sella_kwargs['internal'] == True and len(mol.symbols) < 5", True),
    ('len(mol.symbols) < 5', True),
    ('len(mol.symbols) >= 5', False),
    ('some_unknown_condition', False),
    ('True', False),
])
def test_generated_small_molecule_switch_is_evaluated_from_atom_count(ordinary,condition,eligible):
    qc,_=stalled_complex(ordinary)
    p=Path('well.py')
    p.write_text(p.read_text().replace('opt=Sella',
        f"if {condition}:\n    sella_kwargs['internal']=False\nopt=Sella"))
    proof=qc_optimizer_recovery.evidence(qc,'well','42',cartesian_complex=True)
    assert (proof is not None) is eligible
