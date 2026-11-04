"""Native input regressions and independent-driver reuse lifecycle."""
import itertools
import json
import multiprocessing
import os
from pathlib import Path
import pickle
from types import SimpleNamespace
from unittest.mock import Mock

from ase import Atoms
from ase.db import connect
import numpy as np
import pytest

from kinbot import qc_reuse as reuse
from kinbot.qc import QuantumChemistry
from kinbot.execution import atomic_json
from kinbot.spatial import proper_alignment

DATA = json.loads((Path(__file__).parent/'data/qc_reuse_pairs.json').read_text())
TPL = Path(reuse.__file__).parent/'tpl'


def worker(root, *, job='hir/species_hir_0_00', evidence=None, **changes):
    root.mkdir(parents=True, exist_ok=True)
    e = json.loads(json.dumps(evidence or DATA[0]))
    e.update(changes)
    kwargs = dict(e['kwargs'], label=job, chk=Path(job).name)
    kwargs.update(e.get('extra_kwargs', {}))
    template = 'ase_sella_hir.tpl.py' if job.startswith('hir/') else 'ase_sella_opt_well.tpl.py'
    text = (TPL/template).read_text().format(label=job, checkpoint=Path(job).name,
        atom=e['symbols'], geom=e['positions'], kwargs=kwargs, code='gaussian', Code='Gaussian',
        fix=e.get('fixed_atoms_one_based', [5,1,2,3]), sella_kwargs=e['sella_kwargs'],
        order=e.get('order', 0), fmax=e.get('fmax', .0005), steps=250,
        frequency_mode=e.get('frequency_mode', 'native_hessian'))
    path=root/(job+'.py'); path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text)
    return job


def backend(root):
    q = QuantumChemistry.__new__(QuantumChemistry)
    q.qc='gauss'; q.queuing='slurm'; q.use_sella=True; q.qc_command='g16'
    q.par={}; q.job_ids={}; q.job_generation_tagged={}
    q.db=connect(str(root/'kinbot.db'))
    return q


def launch(q, job, ident='42'):
    q.job_ids[job]=ident
    Path(job+'.sbatch').write_text('submitted '+ident)
    return 0


def complete(root, job, ident='42', *, status='normal'):
    e=DATA[0]; data=dict(status=status,energy=-10.,zpe=.02,frequencies=[100.,200.],
                        hess=np.eye(3*len(e['symbols'])))
    payload=dict(sym=e['symbols'],pos=np.asarray(e['positions'])+[.1,.2,.3],name=job,data=data)
    (root/(job+'.pkl')).write_bytes(pickle.dumps(payload))
    (root/(job+'_sella.log')).write_text('converged\ndone\n')
    (root/(job+'.log')).write_text('Normal termination of Gaussian 16\n')
    (root/(job+'.chk')).write_bytes(b'native checkpoint')
    (root/(job+'.com')).write_text('native original input')
    (root/(job+'.exitcode')).write_text(ident+' 0\n')
    return payload


@pytest.fixture
def pair(tmp_path, monkeypatch):
    monkeypatch.setattr(reuse, 'runtime', lambda: {'test_revision': 1})
    atomic_json(tmp_path/'.kinbot_resource_policy.json', dict(schema=1,max_jobs=20,excluded_nodes=[]))
    a,b=tmp_path/'a',tmp_path/'b'
    job=worker(a); worker(b, evidence=DATA[1])
    qa,qb=backend(a),backend(b)
    return SimpleNamespace(root=tmp_path,a=a,b=b,job=job,qa=qa,qb=qb,mp=monkeypatch)


def enroll(t):
    t.mp.chdir(t.a)
    assert reuse.submit(t.qa,t.job,1,lambda:launch(t.qa,t.job))==0
    t.mp.chdir(t.b)
    duplicate=Mock(side_effect=AssertionError('duplicate physical submission'))
    assert reuse.submit(t.qb,t.job,1,duplicate)==0
    duplicate.assert_not_called()


def test_native_pairs(pair):
    t=pair; t.mp.chdir(t.a); a=reuse.definition(t.qa,t.job)
    t.mp.chdir(t.b); b=reuse.definition(t.qb,t.job)
    assert reuse.equivalent(a,b)
    assert proper_alignment(a['positions'],b['positions']).rmsd==pytest.approx(1.3779413778534112e-7)
    # The other pair is not equivalent even after every elemental mapping.
    a,b=DATA[2:]; candidates=[]
    for perm in itertools.permutations(range(6)):
        if [b['symbols'][i] for i in perm] == a['symbols']:
            candidates.append(proper_alignment(a['positions'],np.asarray(b['positions'])[list(perm)]).rmsd)
    assert min(candidates)==pytest.approx(.8801379054790381)


@pytest.mark.parametrize('change',[{'extra_kwargs':{'mult':4}}, {'extra_kwargs':{'basis':'cc-pVDZ'}},
    {'extra_kwargs':{'charge':1}}, {'extra_kwargs':{'method':'B3LYP'}},
    {'fixed_atoms_one_based':[6,1,2,3]}, {'fmax':.001}, {'order':1},
    {'sella_kwargs':{'internal':False}}, {'extra_kwargs':{'EmpiricalDispersion':'GD3'}},
    {'extra_kwargs':{'guess':'Read'}}])
def test_scientific_changes_do_not_share(pair,change):
    t=pair; t.mp.chdir(t.a); a=reuse.definition(t.qa,t.job)
    worker(t.b, **change);t.mp.chdir(t.b);b=reuse.definition(t.qb,t.job)
    assert b is None or not reuse.equivalent(a,b)


def test_resources_labels_rigid_motion_and_reflection(pair):
    t=pair;t.mp.chdir(t.a);a=reuse.definition(t.qa,t.job)
    xyz=np.asarray(DATA[0]['positions']); rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    job=worker(t.b,job='hir/other_name_2_03',positions=(xyz@rotation+3).tolist(),
               extra_kwargs={'mem':'1000MW','nprocshared':16})
    t.mp.chdir(t.b); assert reuse.equivalent(a,reuse.definition(t.qb,job))
    worker(t.b,job=job,positions=(xyz*[-1,1,1]).tolist())
    assert not reuse.equivalent(a,reuse.definition(t.qb,job))
    xyz[0,0]+=.01;worker(t.b,job=job,positions=xyz.tolist())
    assert not reuse.equivalent(a,reuse.definition(t.qb,job))


def test_complete_native_handoff_and_repeat(pair):
    t=pair;enroll(t)
    assert reuse.poll(t.qb,t.job)=='running'
    payload=complete(t.a,t.job)
    original={p.name:p.read_bytes() for p in (t.a/'hir').iterdir()}
    assert t.qb.check_qc(t.job)=='normal'
    assert t.job not in t.qb.job_ids
    assert not Path(t.job+'.exitcode').exists()
    rows=list(t.qb.db.select(name=t.job));assert len(rows)==1
    np.testing.assert_allclose(rows[0].toatoms().positions,payload['pos'])
    np.testing.assert_allclose(rows[0].data['hess'],payload['data']['hess'])
    assert rows[0].data['shared_qc_job_id']=='42'
    assert {p.name:p.read_bytes() for p in (t.a/'hir').iterdir()}==original
    assert t.qb.check_qc(t.job)=='normal'
    assert len(list(t.qb.db.select(name=t.job)))==1
    # Late consumers use the frozen evidence, not mutable producer artifacts.
    c=t.root/'c';worker(c);qc=backend(c);t.mp.chdir(c)
    assert reuse.submit(qc,t.job,1,Mock(side_effect=AssertionError))==0
    assert qc.check_qc(t.job)=='normal'


def test_producer_ingests_first_follower_uses_readonly_database(pair):
    t=pair;enroll(t);complete(t.a,t.job)
    t.mp.chdir(t.a)
    # Ordinary producer checker consumes the handoff; the cache can still read its DB.
    assert QuantumChemistry._check_qc(t.qa,t.job)=='normal'
    assert not Path(t.job+'.pkl').exists()
    reuse.observe(t.qa,t.job,'normal')
    t.mp.chdir(t.b);assert t.qb.check_qc(t.job)=='normal'


def test_exit_zero_does_not_accept_failed_science(pair):
    t=pair;enroll(t);complete(t.a,t.job,status='error')
    assert reuse.poll(t.qb,t.job)=='running' # owner may still recover
    assert not list(t.qb.db.select(name=t.job))
    t.mp.chdir(t.a);reuse.observe(t.qa,t.job,'error')
    t.mp.chdir(t.b);assert reuse.poll(t.qb,t.job)=='error'
    assert reuse.poll(t.qb,t.job)=='error'
    assert not list(t.qb.db.select(name=t.job))


@pytest.mark.parametrize('artifact',['.sbatch','.exitcode','.log','.pkl','database','job_id'])
def test_activation_preserves_existing_work(pair,artifact):
    t=pair;t.mp.chdir(t.a)
    if artifact=='database':t.qa.db.write(Atoms('H'),name=t.job,data={'status':'error'})
    elif artifact=='job_id':t.qa.job_ids[t.job]='old'
    else:Path(t.job+artifact).write_text('original')
    normal=Mock(return_value=77)
    assert reuse.submit(t.qa,t.job,1,normal)==77
    assert not reuse._marker(t.job,'producer').exists()


def test_crash_after_submission_recovers_id_without_resubmitting(pair):
    t=pair;t.mp.chdir(t.a)
    def interrupted():
        atomic_json(t.root/'.kinbot_policy_submissions/42.json',dict(job_id='42',
            directory=str(t.a),command=['sbatch',t.job+'.sbatch'],created_at='9999'))
        raise RuntimeError('controller interrupted after sbatch')
    with pytest.raises(RuntimeError):reuse.submit(t.qa,t.job,1,interrupted)
    t.mp.chdir(t.b);reuse.submit(t.qb,t.job,1,Mock(side_effect=AssertionError))
    assert reuse.poll(t.qb,t.job)=='running'
    complete(t.a,t.job)
    assert reuse.poll(t.qb,t.job)=='normal'


def test_pre_submission_crash_owner_resumes_and_ambiguous_admission_waits(pair):
    t=pair;t.mp.chdir(t.a)
    with pytest.raises(RuntimeError):
        reuse.submit(t.qa,t.job,1,Mock(side_effect=RuntimeError('before admission')))
    t.mp.chdir(t.b);reuse.submit(t.qb,t.job,1,Mock(side_effect=AssertionError))
    assert reuse.poll(t.qb,t.job)=='running'
    atomic_json(t.root/'.kinbot_policy_admissions.json',{'ambiguous':dict(
        directory=str(t.a),command=['sbatch',t.job+'.sbatch'],created_at='9999',status='submitting')})
    assert reuse.poll(t.qb,t.job)=='running'
    t.mp.chdir(t.a);reuse.submit(t.qa,t.job,1,lambda:launch(t.qa,t.job))
    complete(t.a,t.job);t.mp.chdir(t.b);assert reuse.poll(t.qb,t.job)=='normal'


def test_interrupted_copy_is_idempotent(pair):
    t=pair;enroll(t);complete(t.a,t.job)
    write=t.qb.db.write
    t.mp.setattr(t.qb.db,'write',Mock(side_effect=RuntimeError('interrupt before row')))
    with pytest.raises(RuntimeError):reuse.poll(t.qb,t.job)
    assert not list(t.qb.db.select(name=t.job))
    t.mp.setattr(t.qb.db,'write',write)
    assert reuse.poll(t.qb,t.job)=='normal'
    assert len(list(t.qb.db.select(name=t.job)))==1


def test_changed_source_is_cache_ineligible_not_channel_rejection(pair):
    t=pair;enroll(t);worker(t.a,extra_kwargs={'basis':'cc-pVDZ'})
    complete(t.a,t.job)
    assert reuse.poll(t.qb,t.job) is None
    own=Mock(return_value=0);reuse.submit(t.qb,t.job,1,own);own.assert_called_once()


def test_changed_consumer_and_existing_artifacts_never_overwritten(pair):
    t=pair;enroll(t);complete(t.a,t.job)
    worker(t.b,extra_kwargs={'mult':4})
    with pytest.raises(RuntimeError,match='definition changed'):reuse.poll(t.qb,t.job)
    worker(t.b,evidence=DATA[1]);Path(t.job+'.chk').write_bytes(b'unrelated')
    with pytest.raises(RuntimeError,match='existing artifact'):reuse.poll(t.qb,t.job)
    assert Path(t.job+'.chk').read_bytes()==b'unrelated'
    assert not list(t.qb.db.select(name=t.job))


def test_tampered_cache_never_accepted(pair):
    t=pair;enroll(t);complete(t.a,t.job)
    receipt=reuse._read(reuse._marker(t.job,'consumer'))
    entry=reuse._finish(t.root,receipt['token'])
    (t.root/entry['bundle']/'native.chk').write_bytes(b'changed')
    with pytest.raises(RuntimeError,match='artifact changed'):reuse.poll(t.qb,t.job)
    assert not list(t.qb.db.select(name=t.job))


def _race(root,name,event):
    os.chdir(Path(root)/name); q=backend(Path.cwd());job='hir/species_hir_0_00'
    event.wait(10)
    def physical():
        # O_APPEND record witnesses each actual submission, independent of registry.
        with (Path(root)/'physical_submissions').open('a') as f:f.write(name+'\n')
        return launch(q,job)
    reuse.submit(q,job,1,physical)


def test_simultaneous_workers_submit_once(pair):
    t=pair;ctx=multiprocessing.get_context('spawn');event=ctx.Event()
    ps=[ctx.Process(target=_race,args=(str(t.root),name,event)) for name in ('a','b')]
    for p in ps:p.start()
    event.set()
    for p in ps:
        p.join(30)
        if p.is_alive():p.terminate();p.join()
        assert p.exitcode==0
    assert len((t.root/'physical_submissions').read_text().splitlines())==1


@pytest.mark.parametrize('kwargs',[{'field':'X+10'}, {'basis':'gen'}, {'oldchk':'other.chk'}])
def test_external_frame_or_dependencies_opt_out(pair,kwargs):
    t=pair;worker(t.a,extra_kwargs=kwargs);t.mp.chdir(t.a)
    assert reuse.definition(t.qa,t.job) is None


@pytest.mark.parametrize('change',['frequency_mode','order','isotope','runtime','nonfinite'])
def test_other_definition_boundaries(pair,change):
    t=pair;job='conf/example_0001';worker(t.a,job=job);worker(t.b,job=job)
    t.mp.chdir(t.a);a=reuse.definition(t.qa,job)
    t.mp.chdir(t.b)
    if change=='frequency_mode':worker(t.b,job=job,frequency_mode='finite_difference')
    elif change=='order':worker(t.b,job=job,order=1)
    elif change=='isotope':
        p=Path(job+'.py');p.write_text(p.read_text().replace('kwargs = ', 'mol.set_masses([13.,32.,32.,1.,1.,1.])\nkwargs = ',1))
    elif change=='runtime':t.mp.setattr(reuse,'runtime',lambda:{'test_revision':2})
    else:
        p=Path(job+'.py');p.write_text(p.read_text().replace('1.078494975096588', 'float("nan")'))
    b=reuse.definition(t.qb,job)
    assert b is None or not reuse.equivalent(a,b)


def test_source_replacement_releases_alias_without_cancelling(pair):
    t=pair;enroll(t);t.mp.chdir(t.a)
    t.qa.job_ids[t.job]='43';reuse.observe(t.qa,t.job,'running')
    t.mp.chdir(t.b);assert reuse.poll(t.qb,t.job) is None
    assert t.qa.job_ids[t.job]=='43' and not t.qb.job_ids
    assert list(Path('hir').glob('*.released-*'))


def test_explicit_invalidation_preserves_shared_history(pair):
    t=pair;enroll(t);complete(t.a,t.job);assert reuse.poll(t.qb,t.job)=='normal'
    t.qb.invalidate_qc(t.job)
    assert not reuse._marker(t.job,'consumer').exists()
    assert list(Path('hir').glob('*.qc_reuse_consumer.json.restart_*'))
    assert len(list(t.qb.db.select(name=t.job)))==2
    assert t.qb.check_qc(t.job)==0


def test_profiled_router_reloads_consumer_without_foreign_job_id(pair):
    from kinbot.parameters import Parameters
    from kinbot.profiled_qc import ProfiledQuantumChemistry
    t=pair;enroll(t);t.mp.chdir(t.b)
    inp=Path('input.json');atomic_json(inp,dict(theory_preset='b3lyp-b2plyp-anl',
        queuing='slurm',barrier_threshold=150))
    par=Parameters(str(inp),show_warnings=False).par
    from kinbot.run_format import SCHEMA, rdkit_runtime
    atomic_json(Path('.kinbot_run.json'),dict(schema=SCHEMA,**rdkit_runtime()))
    router=ProfiledQuantumChemistry(par)
    router._record(t.job,'l2',None)
    assert router.check_qc(t.job)=='running'
    complete(t.a,t.job)
    restored=ProfiledQuantumChemistry(par)
    assert restored._level_for(t.job)=='l2'
    assert restored.check_qc(t.job)=='normal'
    assert t.job not in restored.l2.job_ids
    status,energy=restored.get_qc_energy(t.job)
    assert status==0 and np.isfinite(energy)


def test_real_submission_path_shares_one_admission(pair):
    t=pair
    for q in (t.qa,t.qb):
        q.queue_job_limit=0;q.queue_name='test';q.slurm_feature=''
        q.par={'queue_template':'','q_temp_am1':'','q_temp_mp2':'','q_temp_hi':''}
        q._worker_site_setup=lambda:''
    t.mp.setattr('kinbot.qc._ordinary_slurm_walltime',lambda *a:'')
    admission=Mock(return_value=SimpleNamespace(returncode=0,stdout=b'Submitted batch job 42\n',stderr=b''))
    t.mp.setattr('kinbot.run_policy.submit',admission)
    t.mp.chdir(t.a); assert t.qa.submit_qc(t.job,8)==1
    t.mp.chdir(t.b); assert t.qb.submit_qc(t.job,8)==0
    assert admission.call_count==1
    assert t.qa.job_ids=={t.job:'42'} and not t.qb.job_ids
    complete(t.a,t.job)
    assert t.qb.check_qc(t.job)=='normal'


def test_transient_producer_script_during_refresh_does_not_release_alias(pair):
    t=pair;enroll(t)
    p=t.a/(t.job+'.py');text=p.read_text();p.write_text('mol = Atoms(')
    assert reuse.poll(t.qb,t.job)=='running'
    assert reuse._marker(t.job,'consumer').exists()
    complete(t.a,t.job)
    assert reuse.poll(t.qb,t.job)=='running'
    p.write_text(text)
    assert reuse.poll(t.qb,t.job)=='normal'
