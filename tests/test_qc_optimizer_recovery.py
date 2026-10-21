import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

from kinbot import qc_optimizer_recovery as recovery, qc_timeout, timeout_retry
from kinbot.execution import atomic_json, file_sha256
from test_timeout_recovery import pool, ordinary, allocation


def failed_worker(ordinary, monkeypatch, *, ode=False):
    qc,calls=ordinary
    monkeypatch.setattr(timeout_retry,'slurm_job_record',lambda j:{})
    atomic_json(Path('.kinbot_timeout_allocations/42.json'),allocation(Path.cwd(),status='RUNNING'))
    worker="""from ase import Atoms
from sella import Sella
mol = Atoms(symbols=['S', 'H'], positions=[[0, 0, 0], [0, 0, 1.4]])
kwargs = {'method':'B3LYP','scf':'xqc,conver=8','mult':2}
sella_kwargs = {'internal': True}
opt = Sella(mol, order=1, trajectory='well.traj', **sella_kwargs)
converged = opt.run(fmax=0.0001, steps=250)
# Independent frequency acceptance remains mandatory.
if converged:
    validate_frequency(mol)
"""
    Path('well.py').write_text(worker)
    token=hashlib.sha256(b'well').hexdigest()[:20]
    baseline=Path('.kinbot_slurm_io_retries')/token/'42'
    (baseline/'initial/well.py').write_text(worker)
    record=json.loads((baseline/'execution.json').read_text())
    record['initial_artifacts']['well.py']=file_sha256('well.py')
    atomic_json(baseline/'execution.json',record)
    points=[]
    for r,e in [(1.4,-10.),(1.36,-10.1),(1.37,-10.09)]:
        a=Atoms('SH',positions=[[0,0,0],[0,0,r]])
        a.calc=SinglePointCalculator(a,energy=e,forces=np.zeros((2,3)));points.append(a)
    write('well.traj',points)
    reason=('Optimization failed: RuntimeError: Geometry update ODE is taking too long to converge!'
            if ode else 'Optimization did not converge within 250 steps.')
    Path('well_sella.log').write_text('Sella 250 12:00:00 -10.100000 .001 0 .01 1 1\nWorker failure: '+reason+'\n')
    return qc,calls


@pytest.mark.parametrize('ode',[False,True])
def test_ordinary_continuation_archives_and_submits_once_with_saved_geometry(ordinary,monkeypatch,ode):
    qc,calls=failed_worker(ordinary,monkeypatch,ode=ode)
    assert qc_timeout.recover(qc,'well')=='running'
    assert len(calls)==1 and qc.job_ids['well']=='45'
    record=json.loads(qc_timeout._marker('well').read_text())
    archive=Path(record['attempts'][-1]['archive'])
    source=(archive/'retry_worker.py').read_text()
    ast.parse(source)
    assert '[0.0, 0.0, 1.36]' in source
    assert 'chunks=2' in source and 'used=1' in source
    assert 'fmax=0.0001' in source and 'order=1' in source
    assert "'scf': 'xqc,conver=8'" in source and 'validate_frequency(mol)' in source
    assert ("sella_kwargs['internal'] = False" in source)==ode
    assert (archive/'files/well.traj').exists()
    assert Path('other.chk').read_text()=='sibling'
    assert qc_timeout.resume(qc,'well') is None and len(calls)==1
    assert qc_timeout.recover(qc,'well') is None and len(calls)==1


def test_frequency_rejection_and_exhausted_budget_are_never_auto_accepted(ordinary,monkeypatch):
    qc,_=failed_worker(ordinary,monkeypatch)
    atomic_json(Path('well.optimization_budget.json'),dict(chunks_used=3))
    assert recovery.evidence(qc,'well','42') is None
    Path('well.optimization_budget.json').unlink()
    Path('well_sella.log').write_text('Worker failure: Frequency validation failed for transition state: no imaginary modes\n')
    assert recovery.evidence(qc,'well','42') is None


def test_unverified_trajectory_and_constraints_prevent_recovery(ordinary,monkeypatch):
    qc,_=failed_worker(ordinary,monkeypatch)
    source=Path('well.py').read_text().replace("order=1,", "order=1, constraints=constraints,")
    Path('well.py').write_text(source)
    assert recovery.evidence(qc,'well','42') is None


def test_legacy_diatomic_failure_gets_one_new_cause_specific_attempt(ordinary,monkeypatch):
    qc,calls=ordinary;qc.qc='gauss'
    monkeypatch.setattr(timeout_retry,'slurm_job_record',lambda j:{})
    atomic_json(Path('.kinbot_timeout_allocations/42.json'),allocation(Path.cwd(),status='RUNNING'))
    atomic_json(Path('well.scf_recovery.json'),dict(status='failed',reason='Gaussian QC-SCF failure in Link 508'))
    Path('well.log').write_text('Convergence criterion not met.\nSCFChk: SCF convergence 3D-6 required 1D-8\nError termination via Lnk1e in /opt/gaussian/g16/l701.exe\n')
    Path('well.com').write_text('#p UB3LYP/cc-pVDZ force\n\nSH\n\n0 2\nS 0.0 0.0 0.0\nH 0.0 0.0 1.36\n\n')
    source="kwargs = {'mult':2,'scf':'xqc,conver=8','method':'B3LYP','basis':'cc-pVDZ'}\n"
    Path('well.py').write_text(source)
    token=hashlib.sha256(b'well').hexdigest()[:20];baseline=Path('.kinbot_slurm_io_retries')/token/'42'
    (baseline/'initial/well.py').write_text(source)
    record=json.loads((baseline/'execution.json').read_text());record['initial_artifacts']['well.py']=file_sha256('well.py')
    atomic_json(baseline/'execution.json',record)
    assert qc_timeout.recover(qc,'well')=='running' and len(calls)==1
    m=json.loads(qc_timeout._marker('well').read_text())
    assert m['attempts'][0]['source']=='diatomic_symmetry'
    assert qc_timeout.resume(qc,'well') is None and len(calls)==1
