"""Two independent minima, mapped L2/HIR and restart; no new native QC claims."""
import json
import pickle
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.data import atomic_numbers
from ase.units import Bohr

from kinbot import qc_reuse as reuse, minimum_reuse as minima
from kinbot.calculation import load_calculation_record
from kinbot.execution import atomic_json
from tests.test_qc_reuse import DATA, worker, backend, launch

P = [3, 4, 0, 1, 2, 5]  # target atom i is native atom P[i]
ROT = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])


def native_log(symbols, xyz, freq):
    """Synthetic normally completed spectrum in Gaussian's printed format."""
    lines = [' Input orientation:', ' -----', ' Center Atomic Atomic Coordinates',
             ' Number Number Type X Y Z', ' -----']
    lines += [f'{i+1} {atomic_numbers[s]} 0 '+ ' '.join(f'{v:.9f}' for v in x)
              for i, (s,x) in enumerate(zip(symbols, xyz))]
    lines += [' -----', ' Harmonic frequencies (cm**-1)']
    modes = np.random.default_rng(14).normal(size=(len(freq), len(symbols), 3))
    for start in range(0, len(freq), 3):
        f = freq[start:start+3]
        lines += [' Frequencies -- '+' '.join(map(str, f)), ' Atom AN X Y Z']
        lines += [f'{i+1} {atomic_numbers[s]} '+ ' '.join(map(str, modes[start:start+len(f),i,:].ravel()))
                  for i,s in enumerate(symbols)]
    lines += [' Normal termination of Gaussian 16']
    return '\n'.join(lines)+'\n', modes


def source(root, q, symbols, xyz, *, job='conf/selected_0001', energy=-100., freq=None):
    freq = np.arange(1, 3*len(symbols)-5)*100. if freq is None else np.asarray(freq)
    worker(root, job=job, symbols=symbols, positions=np.asarray(xyz).tolist(),
           extra_kwargs={'method':'B3LYP', 'basis':'cc-pVDZ'})
    log, _ = native_log(symbols, xyz, freq)
    (root/(job+'.log')).write_text(log)
    (root/(job+'_sella.log')).write_text('converged\ndone\n')
    row = q.db.write(Atoms(symbols, positions=xyz), name=job,
                    data=dict(status='normal', energy=energy, zpe=.02, frequencies=freq))
    return SimpleNamespace(wellorts=0, atom=np.asarray(symbols), geom=np.asarray(xyz),
                           source_job=job, source_row_id=row, charge=0, mult=2)


def pair(tmp_path, monkeypatch):
    monkeypatch.setattr(reuse, 'runtime', lambda: {'test_revision':1})
    atomic_json(tmp_path/'.kinbot_resource_policy.json', dict(schema=1,max_jobs=20,excluded_nodes=[]))
    a,b=tmp_path/'a',tmp_path/'b';a.mkdir();b.mkdir()
    qa,qb=backend(a),backend(b)
    symbols=DATA[0]['symbols'];x=np.asarray(DATA[0]['positions'])
    # Independent optimizer stops differ slightly; not exact-request reuse.
    y=(x[P]@ROT+3);y[0,0]+=2.e-5
    sa=source(a,qa,symbols,x);sb=source(b,qb,[symbols[i] for i in P],y,energy=-100.+1.e-7)
    return SimpleNamespace(root=tmp_path,a=a,b=b,qa=qa,qb=qb,sa=sa,sb=sb,mp=monkeypatch)


def enroll_minima(t):
    for root,q,s in ((t.a,t.qa,t.sa),(t.b,t.qb,t.sb)):
        t.mp.chdir(root)
        worker(root,job='species_well_high',symbols=s.atom.tolist(),positions=s.geom.tolist())
        minima.attach(q,s,'species_well_high')
        assert reuse.definition(q,'species_well_high').get('minimum')
        callback = (lambda:launch(q,'species_well_high')) if root==t.a else Mock(side_effect=AssertionError('duplicate'))
        reuse.submit(q,'species_well_high',1,callback)
    receipt=reuse._read(reuse._marker('species_well_high','consumer'))
    assert receipt['mapping']['permutation']==P
    assert receipt['mapping']['minimum_comparison']
    assert not t.qb.job_ids
    return receipt


def finish_l2(t, *, stored_hessian=True):
    job='species_well_high'; symbols=t.sa.atom.tolist();xyz=t.sa.geom+[.1,.2,.3]
    freq=np.arange(1,3*len(symbols)-5)*100.
    log,modes=native_log(symbols,xyz,freq)
    t.a.joinpath(job+'.log').write_text(log)
    t.a.joinpath(job+'_sella.log').write_text('converged\ndone\n')
    t.a.joinpath(job+'.exitcode').write_text('42 0\n')
    h=np.random.default_rng(4).normal(size=(18,18));h=h.T@h
    data=dict(status='normal',energy=-105.,zpe=.023,frequencies=freq)
    if stored_hessian:data['hess']=h
    t.a.joinpath(job+'.pkl').write_bytes(pickle.dumps(dict(name=job,sym=symbols,pos=xyz,data=data)))
    if not stored_hessian:
        values=[('Atomic numbers',[atomic_numbers[s] for s in symbols]),
                ('Current cartesian coordinates',(xyz/Bohr).ravel()),
                ('Cartesian Force Constants',h[np.tril_indices(18)])]
        lines=[]
        for name, flat in values:
            lines.append(f'{name:<43} R   N= {len(flat):12}')
            lines.extend(' '.join(map(str,flat[i:i+5])) for i in range(0,len(flat),5))
        t.a.joinpath(job+'.fchk').write_text('\n'.join(lines)+'\n')
    return xyz,h,modes


def test_independent_minima_share_mapped_l2_and_hir(tmp_path,monkeypatch):
    t=pair(tmp_path,monkeypatch);enroll_minima(t)
    xyz,h,modes=finish_l2(t)
    assert t.qb.check_qc('species_well_high')=='normal'
    row=list(t.qb.db.select(name='species_well_high'))[-1]
    np.testing.assert_allclose(row.positions,xyz[P])
    expected=h.reshape(6,3,6,3)[P,:,:,:][:,:,P,:].reshape(18,18)
    np.testing.assert_allclose(t.qb.read_qc_hess('species_well_high',6),expected)
    # Quadratic energy and native displacement vectors are invariant under mapping.
    displacement=np.random.default_rng(9).normal(size=(6,3))
    assert displacement[P].ravel()@expected@displacement[P].ravel()==pytest.approx(displacement.ravel()@h@displacement.ravel())
    mapped=reuse.mapped_modes(row,geometry=row.positions@ROT,symbols=row.symbols)
    np.testing.assert_allclose(mapped['modes'],modes[:,P,:]@ROT)
    # Load through the owning pipeline and generate mapped HIR requests.
    load_calculation_record(t.sb,t.qb,'species_well_high')
    proof=minima.capture(t.qb,t.sb);assert proof is not None
    t.mp.chdir(t.a);assert t.qa.check_qc('species_well_high')=='normal'
    load_calculation_record(t.sa,t.qa,'species_well_high')
    fix=[5,1,2,3];mapped_fix=[P.index(i-1)+1 for i in fix]
    for root,q,s,indices in ((t.a,t.qa,t.sa,fix),(t.b,t.qb,t.sb,mapped_fix)):
        t.mp.chdir(root)
        job=worker(root,job='hir/independent_hir_1_03',symbols=s.atom.tolist(),positions=s.geom.tolist(),
                   fixed_atoms_one_based=indices)
        minima.attach(q,s,job)
        callback=(lambda:launch(q,job,'43')) if root==t.a else Mock(side_effect=AssertionError('duplicate HIR'))
        reuse.submit(q,job,1,callback)
    assert reuse.poll(t.qb,job)=='running'
    assert reuse._read(reuse._marker(job,'consumer'))['mapping']['permutation']==P
    # Invalidating a reused result preserves the mapped/native evidence and resets derivatives.
    t.qb.invalidate_qc('species_well_high')
    assert t.qb.check_qc('species_well_high')==0
    assert t.qb.db.get(id=row.id).data['shared_qc_native']['permutation']==P


def test_native_hessian_mapping_and_conventional_publication(tmp_path,monkeypatch):
    t=pair(tmp_path,monkeypatch);enroll_minima(t);xyz,h,_=finish_l2(t,stored_hessian=False)
    assert t.qb.check_qc('species_well_high')=='normal'
    expected=reuse.permute_hessian(h,P)
    np.testing.assert_allclose(t.qb.read_qc_hess('species_well_high',6),expected,atol=1.e-12)
    row=list(t.qb.db.select(name='species_well_high'))[-1]
    t.qb.publish_result(row,'conventional_well_high')
    np.testing.assert_allclose(t.qb.read_qc_hess('conventional_well_high',6),expected,atol=1.e-12)
    assert t.qb.check_qc('conventional_well_high')=='normal'


@pytest.mark.parametrize('change',['energy','zpe','frequency','theory','geometry','isotope','no_native','incomplete','saddle','generation'])
def test_failed_minimum_identity_is_not_qc_rejection(tmp_path,monkeypatch,change):
    t=pair(tmp_path,monkeypatch);t.mp.chdir(t.a)
    a=minima.capture(t.qa,t.sa);assert a is not None
    t.mp.chdir(t.b);b=minima.capture(t.qb,t.sb);assert b is not None
    if change=='energy': b['energy_hartree']+=.001
    elif change=='zpe':b['zpe_hartree']+=.001
    elif change=='frequency':b['frequencies_cm1'][0]+=10
    elif change=='theory':b['theory_program']='different'
    elif change=='geometry':b['positions'][0][0]+=.01
    elif change=='isotope':b['masses'][0]+=1
    else:
        path=Path(t.sb.source_job+'.log')
        if change=='no_native':path.unlink()
        elif change=='incomplete':path.write_text(path.read_text().replace('Normal termination','Error termination'))
        elif change=='saddle':
            row=t.qb.db.get(id=t.sb.source_row_id);data=dict(row.data);data['frequencies']=list(data['frequencies']);data['frequencies'][0]=-100.
            t.qb.db.update(row.id,data=data)
        else:
            t.qb.job_ids[t.sb.source_job]='43';Path(t.sb.source_job+'.exitcode').write_text('42 0\n')
        assert minima.capture(t.qb,t.sb) is None
        return
    assert minima.comparison(b,a,P) is None


def test_loose_coordinate_matching_requires_minimum_proof_and_correct_constraint(tmp_path,monkeypatch):
    t=pair(tmp_path,monkeypatch);t.mp.chdir(t.a)
    ja=worker(t.a,job='species_well_high',symbols=t.sa.atom.tolist(),positions=t.sa.geom.tolist())
    a=reuse.definition(t.qa,ja)
    t.mp.chdir(t.b);jb=worker(t.b,job=ja,symbols=t.sb.atom.tolist(),positions=t.sb.geom.tolist())
    b=reuse.definition(t.qb,jb)
    assert reuse.matching(b,a) is None
    minima.attach(t.qb,t.sb,jb);b=reuse.definition(t.qb,jb)
    t.mp.chdir(t.a);minima.attach(t.qa,t.sa,ja);a=reuse.definition(t.qa,ja)
    assert reuse.matching(b,a)['permutation']==P
    # Same minima do not justify different torsional constraints/scan angles.
    a['constraint_atoms']=[4,0,1,2];b['constraint_atoms']=[0,2,3,5]
    assert reuse.matching(b,a) is None


def test_prior_version_consumers_survive_refresh(tmp_path,monkeypatch):
    t=pair(tmp_path,monkeypatch)
    job=worker(t.a);worker(t.b,evidence=DATA[1])
    for root,q in ((t.a,t.qa),(t.b,t.qb)):
        t.mp.chdir(root)
        old=reuse.definition(q,job,legacy=True)
        # Recreate exact old definition/receipts, not just old runtime values.
        if root==t.a:
            token='a'*32;atomic_json(reuse._record(t.root,token),dict(token=token,source='a/'+job,
                directory='a',job=job,definition=old,state='submitted',job_id='42'))
        else:atomic_json(reuse._marker(job,'consumer'),dict(root=str(t.root),token=token,definition=old,state='waiting'))
    from tests.test_qc_reuse import complete
    complete(t.a,job)
    assert reuse.poll(t.qb,job)=='normal'


def test_real_optimization_and_hir_dispatch_share_after_minimum_selection(tmp_path,monkeypatch):
    from kinbot.parameters import Parameters
    from kinbot.qc import QuantumChemistry
    from kinbot.optimize import Optimize
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.run_format import SCHEMA, rdkit_runtime
    t=pair(tmp_path,monkeypatch)
    # Mock only the external scheduler/native execution. The profiled router,
    # Optimize acceptance, HIR generation, native readers and reuse are real.
    admissions=[]
    def admission(*args,**kwargs):
        ident=str(42+len(admissions));admissions.append((Path.cwd(),args,kwargs))
        return SimpleNamespace(returncode=0,stdout=f'Submitted batch job {ident}\n'.encode(),stderr=b'')
    monkeypatch.setattr('kinbot.run_policy.submit',admission)
    monkeypatch.setattr('kinbot.qc._ordinary_slurm_walltime',lambda *a:'')
    monkeypatch.setattr('kinbot.qc._slurm_job_active',lambda *a:True)
    monkeypatch.setattr('kinbot.qc._slurm_job_record',lambda *a:{'JobState':'RUNNING'})
    opts=[]
    for root,s in ((t.a,t.sa),(t.b,t.sb)):
        monkeypatch.chdir(root)
        Path('hir').mkdir(exist_ok=True)
        atomic_json(Path('input.json'),dict(theory_preset='b3lyp-b2plyp-anl',queuing='slurm',
            barrier_threshold=150,conformer_search=1,rotor_scan=1,nrotation=4,
            high_level=1,L3_calc=0,multi_conf_tst=0,plot_hir_profiles=0))
        atomic_json(Path('.kinbot_run.json'),dict(schema=SCHEMA,**rdkit_runtime()))
        par=Parameters('input.json',show_warnings=False).par
        par['conformer_search']=0  # resume after the selected L1 rows staged above
        q=QuantumChemistry(par);q.l1._worker_site_setup=lambda:'';q.l2._worker_site_setup=lambda:''
        q._record(s.source_job,'l1')
        point=StationaryPoint('selected',0,2,atom=s.atom.tolist(),geom=s.geom)
        point.characterize();load_calculation_record(point,q,s.source_job)
        opt=Optimize(point,par,q)
        opt.do_optimization();assert opt.shigh==0
        opts.append(opt)
    assert len(admissions)==1
    ja,jb=(o.log_name(1) for o in opts)
    # Refresh/reconstruct the reader while L2 remains in flight; preserve
    # both scripts, source allocation and consumer's absent Slurm ID.
    script=Path(jb+'.py').read_bytes();opts[1].qc.qc_opt(opts[1].species,opts[1].species.geom,high_level=1)
    assert Path(jb+'.py').read_bytes()==script and len(admissions)==1
    # Publish the simulated L2 result under the real generated producer name.
    xyz,h,modes=finish_l2(t)
    for suffix in ('.pkl','.log','_sella.log','.exitcode'):
        p=t.a/('species_well_high'+suffix)
        if suffix=='.pkl':
            payload=pickle.loads(p.read_bytes());payload['name']=ja;p.write_bytes(pickle.dumps(payload))
        p.rename(t.a/(ja+suffix))
    for root,opt in zip((t.a,t.b),opts):
        monkeypatch.chdir(root)
        opt.do_optimization()
        assert opt.shigh==1 and opt.shir==0
        assert opt.species.source_job==opt.log_name(1)
    # Every generated HIR request in the second worker must be a consumer.
    jobs_a=[j for j in opts[0].qc.l2.job_ids if j.startswith('hir/')]
    jobs_b=[j for j in opts[1].qc.l2.job_ids if j.startswith('hir/')]
    assert jobs_a and not jobs_b
    assert len(admissions)==1+len(jobs_a)
    assert len(list((t.b/'hir').glob('*.qc_reuse_consumer.json')))==len(jobs_a)
    # Repeated optimizer polls do not generate additional admissions.
    for root,opt in zip((t.a,t.b),opts):
        monkeypatch.chdir(root);opt.do_optimization()
    assert len(admissions)==1+len(jobs_a)
    # Complete native point handoffs; then exercise scan ingestion, projection,
    # final conventional publication and repeated finished-driver polling.
    monkeypatch.chdir(t.a)
    for i,job in enumerate(jobs_a):
        d=reuse.definition(opts[0].qc.l2,job)
        payload=dict(name=job,sym=d['symbols'],pos=np.asarray(d['positions']),
                     data=dict(status='normal',energy=-105.+.01*(1-np.cos(i*np.pi/2))))
        Path(job+'.pkl').write_bytes(pickle.dumps(payload))
        Path(job+'.log').write_text('Normal termination of Gaussian 16\n')
        Path(job+'_sella.log').write_text('converged\ndone\n')
        Path(job+'.exitcode').write_text(opts[0].qc.l2.job_ids[job]+' 0\n')
    for root,opt in zip((t.a,t.b),opts):
        monkeypatch.chdir(root);opt.do_optimization();opt.do_optimization()
        assert opt.shigh==opt.shir==1
        assert opt.species.source_job==opt.log_name(1)
    np.testing.assert_allclose(opts[0].species.reduced_freqs,opts[1].species.reduced_freqs,rtol=1.e-8)
    assert len(admissions)==1+len(jobs_a)


def test_existing_rotor_definitions_are_not_rewritten(tmp_path,monkeypatch):
    from kinbot.species_routing import routing_name
    from kinbot.stationary_pt import StationaryPoint
    t=pair(tmp_path,monkeypatch)
    points=[]
    for root,q,s in ((t.a,t.qa,t.sa),(t.b,t.qb,t.sb)):
        monkeypatch.chdir(root);Path('hir').mkdir()
        point=StationaryPoint('selected',0,2,atom=s.atom.tolist(),geom=s.geom)
        point.characterize();load_calculation_record(point,q,s.source_job)
        points.append(point)
    monkeypatch.chdir(t.a);minima.canonical_rotors(points[0],t.qa)
    monkeypatch.chdir(t.b)
    before=deepcopy(points[1].dihed)
    # Staging counts as existing work, even before sbatch/Slurm registration.
    staged=Path('hir')/(routing_name(points[1])+'_hir_0_00.py')
    staged.write_text('already staged')
    minima.canonical_rotors(points[1],t.qb)
    assert points[1].dihed==before and staged.read_text()=='already staged'
    staged.unlink()
    minima.canonical_rotors(points[1],t.qb)
    assert [P[i] for i in points[1].dihed[0]]==points[0].dihed[0]
    receipt=Path('hir')/(routing_name(points[1])+'.minimum_mapping.json')
    assert json.loads(receipt.read_text())['original_dihedrals']==before


def test_mapping_tampering_cannot_publish_a_result(tmp_path,monkeypatch):
    t=pair(tmp_path,monkeypatch);receipt=enroll_minima(t);finish_l2(t)
    receipt['mapping']['permutation']=[0,1,2,3,4,5]
    atomic_json(reuse._marker('species_well_high','consumer'),receipt)
    with pytest.raises(RuntimeError,match='atom mapping'):
        reuse.poll(t.qb,'species_well_high')
    assert not list(t.qb.db.select(name='species_well_high'))
