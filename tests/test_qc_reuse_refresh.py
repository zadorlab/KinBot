"""Full controller handoff must retain source ownership and pending aliases."""
from copy import deepcopy
import json
from pathlib import Path

from kinbot import qc_reuse as reuse, run_refresh
from kinbot.execution import atomic_json
from tests.test_run_refresh import live_tree, sampling_tree  # pytest fixtures
from tests.test_qc_reuse import worker, backend, launch, complete


def test_refresh_preserves_shared_source_and_repeated_consumer_adoption(sampling_tree, monkeypatch):
    t=sampling_tree
    a=t.root/'ch3s/123';b=t.root/'ch3ss/123'
    # Replace only the simulated fixture's placeholder databases.
    for p in (a,b):(p/'kinbot.db').unlink()
    job=worker(a);worker(b)
    monkeypatch.setattr(reuse,'runtime',lambda:{'test_revision':1})
    monkeypatch.chdir(a);qa=backend(a);reuse.submit(qa,job,1,lambda:launch(qa,job))
    monkeypatch.chdir(b);qb=backend(b)
    reuse.submit(qb,job,1,lambda:(_ for _ in ()).throw(AssertionError('duplicate')))
    atomic_json(a/'123.json', {})
    atomic_json(a/'.kinbot_theory_jobs.json',dict(schema=1,jobs={job:dict(level='l2',job_id='42',generation_tagged=True)}))
    atomic_json(b/'123.json', {})
    atomic_json(b/'.kinbot_theory_jobs.json',dict(schema=1,jobs={job:dict(level='l2',job_id=None,generation_tagged=False)}))
    t.jobs.append(dict(job_id='42',directory=str(a),state='RUNNING'))
    paths=[a/(job+'.py'),a/(job+'.sbatch'),a/(job+'.qc_reuse_producer.json'),
           b/(job+'.py'),b/(job+'.qc_reuse_consumer.json')]
    before={p:p.read_bytes() for p in paths}; jobs=deepcopy(t.jobs)
    result=run_refresh.refresh_when_idle(t.root,apply=True,wait_idle=30)
    assert result['phase']=='launched' and t.jobs==jobs
    assert {p:p.read_bytes() for p in paths}==before
    assert {p['pid'] for p in t.processes}=={700,701,702} # retained sampling
    assert run_refresh.refresh_when_idle(t.root,apply=True,wait_idle=30)==result
    # New controller can adopt old-code source evidence after the refresh.
    monkeypatch.setattr(reuse,'runtime',lambda:{'test_revision':2})
    restored=backend(b);assert reuse.poll(restored,job)=='running'
    complete(a,job)
    assert reuse.poll(restored,job)=='normal'
    assert len(list(restored.db.select(name=job)))==1
    assert job not in restored.job_ids
    assert reuse.poll(restored,job) is None


def test_refresh_preserves_mapped_minima_l2_and_sampling(sampling_tree, monkeypatch):
    import numpy as np
    from types import SimpleNamespace
    from tests.test_minimum_reuse import source, enroll_minima, finish_l2, P, ROT
    from tests.test_qc_reuse import DATA
    t=sampling_tree;a=t.root/'ch3s/123';b=t.root/'ch3ss/123'
    for p in (a,b):(p/'kinbot.db').unlink()
    qa,qb=backend(a),backend(b)
    atomic_json(t.root/'.kinbot_resource_policy.json',dict(schema=1,max_jobs=20,excluded_nodes=[]))
    symbols=DATA[0]['symbols'];xyz=np.asarray(DATA[0]['positions'])
    other=xyz[P]@ROT+3;other[0,0]+=2.e-5
    sa=source(a,qa,symbols,xyz);sb=source(b,qb,[symbols[i] for i in P],other)
    case=SimpleNamespace(root=t.root,a=a,b=b,qa=qa,qb=qb,sa=sa,sb=sb,mp=monkeypatch)
    monkeypatch.setattr(reuse,'runtime',lambda:{'test_revision':1})
    enroll_minima(case);job='species_well_high'
    for p,ident in ((a,'42'),(b,None)):
        atomic_json(p/'123.json',{})
        atomic_json(p/'.kinbot_theory_jobs.json',dict(schema=1,jobs={job:dict(level='l2',job_id=ident,generation_tagged=bool(ident))}))
    t.jobs.append(dict(job_id='42',directory=str(a),state='RUNNING'))
    before={p:p.read_bytes() for root in (a,b) for p in root.rglob('*') if p.is_file()
            and p.suffix in ('.py','.sbatch') or (p.is_file() and 'reuse' in p.name)}
    jobs=deepcopy(t.jobs)
    result=run_refresh.refresh_when_idle(t.root,apply=True,wait_idle=30)
    assert result['phase']=='launched' and t.jobs==jobs
    assert {p:p.read_bytes() for p in before}==before
    assert {p['pid'] for p in t.processes}=={700,701,702}
    assert run_refresh.refresh_when_idle(t.root,apply=True,wait_idle=30)==result
    monkeypatch.setattr(reuse,'runtime',lambda:{'test_revision':2})
    restored=backend(b);assert reuse.poll(restored,job)=='running'
    xyz,h,_=finish_l2(case,stored_hessian=False)
    assert restored.check_qc(job)=='normal'
    np.testing.assert_allclose(restored.read_qc_hess(job,6),reuse.permute_hessian(h,P),atol=1.e-12)
    assert job not in restored.job_ids
    assert restored.check_qc(job)=='normal' and len(list(restored.db.select(name=job)))==1
