"""Local simulated-site tests. No test here claims a native Molpro calculation."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
import pytest

from kinbot import auto_radical_pair as auto
from kinbot import radical_pair_evidence as ev
from kinbot.execution import atomic_json,file_sha256
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot import molpro_reconnaissance as recon


@pytest.fixture
def case():
    p=Parameters('examples/anl/c2h6_anl1_f12/c2h6.json').par
    s=StationaryPoint('parent',0,1,smiles='CC');s.characterize()
    s.bond[0,1]=s.bond[1,0]=0;products,maps=s.start_multi_molecular();s.bond[0,1]=s.bond[1,0]=1
    r=SimpleNamespace(species=s,products=products,scan_coo=[0,1],instance_name='parent_hom_sci_1_2',family_name='hs')
    return p,r,auto.context(p,r)


def test_example_uses_real_opt_in_and_site_budget(case):
    p,r,c=case
    assert p['vrc_tst_active_space']=='auto_radical_pair' and not p['vrc_tst_surface_model_file']
    assert c['bases']==['avdz-f12','avtz-f12','avqz-f12']
    assert p['anl_f12_augmented']==1 and p['composite_method']=='ANL1-F12'
    assert p['queue_job_limit']==p['queue_live_policy']['max_jobs']==10
    assert p['ppn']==p['single_point_ppn']==p['rotdpy_processors']==8
    assert c['full_space']['active_electrons']==c['full_space']['active_orbitals']==14
    assert len(c['full_space']['orbital_selection']['centers'])==8
    assert c['spin']['gap_levels'][0]['basis']==c['spin']['reference']['basis']==['avtz-f12','avqz-f12']


@pytest.mark.parametrize('change',['spin','charge','products','bond','ring','orders','roots','counts'])
def test_applicability_is_chemical_not_a_molecule_name(case,change):
    p,r,c=case
    if change=='spin':r.species.mult=3
    if change=='charge':r.products[0].charge=1
    if change=='products':r.products[0].bond[:]=0
    if change=='bond':r.scan_coo=[2,3]
    if change=='ring':r.species.bond[2,3]=r.species.bond[3,2]=1
    if change=='orders':r.species.bond[0,1]=r.species.bond[1,0]=2
    if change=='roots':p['vrc_tst_selected_roots']={'1':2}
    if change=='counts':p['vrc_tst_scan_method']='caspt2-f12(4,4)'
    with pytest.raises(auto.ReviewRequired):auto.context(p,r)


def test_c_h_scission_supported_and_actual_products_verified(case):
    p,r,_=case;s=r.species
    h=next(i for i in range(len(s.atom)) if s.atom[i]=='H' and s.bond[0,i])
    s.bond[0,h]=s.bond[h,0]=0;r.products,_=s.start_multi_molecular();s.bond[0,h]=s.bond[h,0]=1
    r.scan_coo=[0,h];r.instance_name='a_hom_sci_1_'+str(h+1)
    c=auto.context(p,r)
    assert sorted(map(len,c['mappings']))==[1,7]


def seed_files(tmp_path,c):
    seeds={}
    for i,basis in enumerate(c['bases']):
        d=tmp_path/str(i);d.mkdir(parents=True)
        inp=d/'scf.inp';inp.write_text('{rhf;wf,18,1,2;orbital,2100.2}\n')
        out=d/'scf.out';out.write_text(' CANONICAL ORBITALS\n Orbital  Occupation    Energy    Cen Mu Typ  Coefficients\n'+
            ''.join(f'{j}.1 {2. if j<=8 else 1. if j<=10 else 0.:.5f} -0.1 1 1 s 1.0\n' for j in range(1,14)))
        wfu=d/'seed.wfu';wfu.write_bytes(b'SYNTHETIC seed, no scientific claim')
        seeds[basis]=dict(path=str(wfu),sha256=file_sha256(wfu),basis=basis,active_indices=[9,10],
            source_input=str(inp),source_input_sha256=file_sha256(inp),source_output=str(out),source_output_sha256=file_sha256(out))
    return seeds


def test_native_decks_use_identified_somos_single_spin_and_propagated_records(case,tmp_path):
    p,r,c=case;seeds=seed_files(tmp_path/'seeds',c)
    points={f'r{i}':dict(geometry_angstrom=auto.move(c,c['geometry'],d).tolist(),distance_angstrom=d)
            for i,d in enumerate(reversed(c['radial']))}
    anchors,angles,pairs=auto.anchors(c,points)
    for full in (False,True):
        for m in (1,3):
            name=('full' if full else 'compact')+str(m)
            source=auto.panel_source(c,p,name,m,anchors,seeds,full)
            file=atomic_json(tmp_path/(name+'.json'),source)
            path=recon.prepare(file,tmp_path/name);manifest=json.loads(path.read_text())
            assert len(manifest['jobs'])==len(anchors)*3
            first,second=manifest['jobs'][:2]
            a=(Path(first['directory'])/(first['name']+'.inp')).read_text()
            b=(Path(second['directory'])/(second['name']+'.inp')).read_text()
            assert f'wf,18,1,{m-1}' in a.lower()
            assert 'wf,18,1,2' in a.lower()  # high-spin initial ROHF, independent of MULTI spin
            assert first['root_states']==[dict(multiplicity=m,root=1)]
            assert second['orbital_predecessor']==first['name'] and not second['initial_record']
            assert 'start,2140.2,type=natural' in b
            if not full:
                assert 'start,2100.2,type=canonical' in a and 'config,csf' in a
                assert 'printci=' in a and 'Source active indices: 9,10' in a
                assert 'config,csf' in b and 'start,2100.2,type=canonical' not in b
                assert len(first['correlation_checks'])==(4 if m==3 else 3)
            else:
                assert first['probe_active_space']['active_orbitals']==14
                assert len(first['correlation_checks'])==1
            script=(Path(first['directory'])/(first['name']+'.sbatch')).read_text()
            assert '--exclusive' in script and '-n 8' in script
    assert len(angles)==12 and pairs
    for a,b in pairs:assert points[a]==points[b]
    for a,b in angles.items():
        assert points[a]['distance_angstrom']==points[b]['distance_angstrom']
        xyz=np.array(points[a]['geometry_angstrom']);assert np.linalg.norm(xyz[c['bond'][0]]-xyz[c['bond'][1]])==pytest.approx(points[a]['distance_angstrom'])


def test_read_native_xml_rejects_incomplete_and_associates_basis(tmp_path):
    p=tmp_path/'orbs.xml'
    p.write_text('''<molpro xmlns="http://www.molpro.net/schema/molpro-output" xmlns:c="http://www.xml-cml.org/schema" xmlns:x="http://www.w3.org/1999/xlink"><molecule id="m1"><c:atomArray><c:atom id="a1" elementType="H" x3="0" y3="0" z3="0"/></c:atomArray><basisSet length="1" cartesianLength="1"><basisGroup id="g1" minL="0" maxL="0"><basisExponents>1.0</basisExponents><basisContraction>1.0</basisContraction></basisGroup><association><bases x:href="#xpointer(//basisGroup[@id='g1'])"/><atoms x:href="#xpointer(//atom[@id='a1'])"/></association></basisSet><orbitals><orbital occupation="1">1</orbital></orbitals></molecule></molpro>''')
    o=ev.read_orbitals(p);assert o.occupations.tolist()==[1.] and o.atom_columns==[[0]]
    p.write_text(p.read_text().replace('length="1"','length="2"'))
    with pytest.raises(ValueError,match='Incomplete'):ev.read_orbitals(p)


def orbitals(matrix,occupations=None):
    n=len(matrix)
    return ev.Orbitals(np.asarray(matrix),np.ones(n) if occupations is None else np.array(occupations),
                       ['H']*n,np.zeros((n,3)),[[(0,(1.,),((1.,),))]]*n,[[i] for i in range(n)])


def test_subspace_and_ci_overlap_are_rotation_invariant_and_detect_wrong_orbitals():
    a=orbitals(np.eye(4));b=orbitals(np.eye(4));q=np.array([[1,1],[1,-1]])/np.sqrt(2)
    b.coefficients[:2,:2]=q
    assert ev.projection(a,b,[0,1],[0,1])==pytest.approx([1.,1.])
    b.positions[0,0]=1
    with pytest.raises(ValueError,match='geometry'):ev.transform(a,b)
    b.positions[0,0]=0;b.coefficients[[0,2]]=b.coefficients[[2,0]]
    assert min(ev.projection(a,b,[0,1],[0,1]))<.01
    covalent={'20':1/np.sqrt(2),'02':-1/np.sqrt(2)}
    assert ev.ionic_probability(covalent,q)==pytest.approx(0,abs=1.e-12)
    assert ev.ionic_probability({'20':1/np.sqrt(2),'02':1/np.sqrt(2)},q)==pytest.approx(1.)
    assert ev.ci_overlap({'ab':1.},covalent,q,1)==pytest.approx(1.)
    assert ev.ci_overlap({'ab':1.},{'20':1.},np.eye(2),1)==0.


def test_independent_fragment_embedding_rejects_wrong_basis_or_charge():
    d=orbitals(np.eye(2));f=orbitals(np.eye(1))
    s,p=ev.fragment_embedding(d,[f,f],[[0],[1]])
    assert p==[1.,1.] and np.array_equal(s,np.eye(2))
    g=deepcopy(f);g.atom_basis=[[]]
    with pytest.raises(ValueError,match='basis'):ev.fragment_embedding(d,[f,g],[[0],[1]])
    f.occupations[0]=2
    with pytest.raises(ValueError,match='SOMO'):ev.fragment_embedding(d,[f,f],[[0],[1]])


def test_ci_parser_requires_natural_transform_complete_vector_and_single_root():
    s=' PROGRAM * MULTI\nNatural orbital dump at molpro section 2140.2 (Orbital set 2)\nTransform the CI vector\n CI Coefficients of symmetry 1\n 20 .70710678\n 02 -.70710678\n Energy: -1\n PROGRAM * RSPT2\n'
    ci=ev.natural_ci(s,1);assert set(ci)=={'20','02'}
    with pytest.raises(ValueError,match='incomplete'):ev.natural_ci(s.replace('02 -.70710678','02 -.1'),1)
    with pytest.raises(ValueError,match='transformed'):ev.natural_ci(s.replace('Transform the CI vector',''),1)
    with pytest.raises(ValueError):ev.natural_ci(s,3)


def selected_jobs(tmp_path,c):
    jobs={}
    for mult in (1,3):
        group=[]
        for b in c['bases']:
            f=tmp_path/(b+str(mult)+'.wfu');f.write_bytes(b'fake native record for orchestration test')
            group.append(dict(basis=b,anchor='r0',wfu=str(f),wfu_sha256=file_sha256(f)))
        jobs[f'compact{mult}']=group
    return jobs


def test_single_spin_records_reach_sampling_gap_cc_and_immutable_levels(case,tmp_path):
    from kinbot.vrc_correction import correction_components,expected_levels
    from kinbot.vrc_multistate import state_parameters,resolve_state_model_paths
    p,r,c=case;jobs=selected_jobs(tmp_path,c);surface=auto.surface(c,jobs,tmp_path/'v.json')
    stage=surface['channels'][0]['bottlenecks'][0]
    leaf=state_parameters(p,stage,stage['states'][0],[2,2])
    components=correction_components(leaf,symbols=c['symbols'],charge=0)
    for item in components:
        m=item['multiplicity'] or 1
        assert item['state_model']['state_average']==[dict(multiplicity=m,roots=1,weights=[1.])]
        assert str(m)+'.wfu' in item['state_model']['orbital_selection']['path']
    levels=expected_levels(leaf,symbols=c['symbols'],charge=0)
    assert set(levels['state_models_by_spin'])=={'1','3'}
    assert resolve_state_model_paths(stage['parameters']['vrc_tst_state_model'],tmp_path)['schema']==3


def test_auto_mode_never_falls_back_to_frontier_selection(case):
    from kinbot.vrc_correction import resolved_state_model
    p,r,c=case
    with pytest.raises(ValueError,match='has not supplied'):resolved_state_model(p,c['symbols'],0,c['spin'],basis=c['bases'][0])


@pytest.fixture
def simulated_workflow(case,tmp_path,monkeypatch):
    p,r,c=case;monkeypatch.chdir(tmp_path);calls=[]
    artifact=tmp_path/'native.out';artifact.write_text('SYNTHETIC integration-test artifact')
    jobs=selected_jobs(tmp_path,c)
    def scf(root,c,p,far):
        calls.append(('scf',root.name));return {b:{} for b in c['bases']},{}
    def runner(root,sources,p):
        calls.append(('run',root.name));return jobs
    def relaxed(root,c,p,r,pilots):
        calls.append(('relax','scan'))
        points={f'r{i}':dict(geometry_angstrom=auto.move(c,c['geometry'],d).tolist(),distance_angstrom=d)
                for i,d in enumerate(reversed(c['radial']))}
        atomic_json(root/'relaxations.json',dict(points=points,receipts=[]))
        return points
    def assessment(c,*args):
        calls.append(('assess','native'))
        return dict(status='accepted',policy=c['policy'],contract_sha256=auto.digest(c),
            checks=[dict(test='synthetic_orchestration',passed=True)],files={str(artifact):file_sha256(artifact)})
    monkeypatch.setattr(auto,'scf_phase',scf);monkeypatch.setattr(auto,'_run_sources',runner)
    monkeypatch.setattr(auto,'relax',relaxed);monkeypatch.setattr(auto,'assess',assessment)
    monkeypatch.setattr(auto,'_files',lambda jobs:{})
    return p,r,c,calls,artifact


def test_complete_owner_orchestration_promotes_after_verification_and_reuses(simulated_workflow):
    from kinbot.vrc_multistate import channel_for_reaction,state_parameters
    from kinbot.vrc_reference_review import validate
    p,r,c,calls,artifact=simulated_workflow
    original=deepcopy(p)
    channel,provenance=channel_for_reaction(p,r)
    assert calls==[('scf','initial_scf'),('run','pilots'),('relax','scan'),('scf','relaxed_scf'),('run','panel'),('assess','native')]
    assert p==original and Path(provenance['verification']).is_file()
    stage=channel['bottlenecks'][0];leaf=state_parameters(p,stage,stage['states'][0],[2,2])
    validate(leaf,Path.cwd())
    before=list(calls)
    assert channel_for_reaction(p,r)==(channel,provenance) and calls==before
    artifact.write_text('tampered')
    from kinbot.vrc_multistate import SurfaceModelRequired
    with pytest.raises(SurfaceModelRequired,match='changed'):channel_for_reaction(p,r)
    assert calls==before


def test_failed_science_is_persistent_and_never_produces_surface(simulated_workflow,monkeypatch):
    from kinbot.vrc_multistate import channel_for_reaction,SurfaceModelRequired
    p,r,c,calls,artifact=simulated_workflow
    def failed(c,*args):
        return dict(status='needs_review',policy=c['policy'],contract_sha256=auto.digest(c),
                    checks=[dict(test='wrong_fragment',passed=False)],files={str(artifact):file_sha256(artifact)})
    monkeypatch.setattr(auto,'assess',failed)
    with pytest.raises(SurfaceModelRequired,match='^Native verification failed'):channel_for_reaction(p,r)
    assert not list(Path.cwd().glob('vrctst/auto_radical_pair/*/surface.json'))
    before=list(calls)
    with pytest.raises(SurfaceModelRequired,match='requires review'):channel_for_reaction(p,r)
    assert before==calls


def test_crash_after_durable_decision_finishes_publication_without_new_calcs(simulated_workflow,monkeypatch):
    p,r,c,calls,artifact=simulated_workflow;publish=auto._publish
    class Crash(BaseException):pass
    def crash(*args):raise Crash()
    monkeypatch.setattr(auto,'_publish',crash)
    with pytest.raises(Crash):auto.ensure(p,r)
    before=list(calls);monkeypatch.setattr(auto,'_publish',publish)
    channel,_=auto.ensure(p,r)
    assert calls==before and channel['bottlenecks'][0]['states'][0]['multiplicity']==1


def test_live_owner_preserved_without_starting_second_panel(simulated_workflow):
    from kinbot.vrc_multistate import SurfaceModelRequired
    p,r,c,calls,_=simulated_workflow
    root=(Path('vrctst/auto_radical_pair')/auto.digest(c)).resolve()
    with recon._manifest_lock(root):
        with pytest.raises(SurfaceModelRequired,match='live controller'):auto.ensure(p,r)
    assert not calls


def test_changed_contract_never_uses_an_old_approval(simulated_workflow):
    from kinbot.vrc_multistate import state_parameters
    from kinbot.vrc_reference_review import validate
    p,r,c,calls,_=simulated_workflow
    channel,_=auto.ensure(p,r);stage=channel['bottlenecks'][0]
    leaf=state_parameters(p,stage,stage['states'][0],[2,2])
    leaf['vrc_tst_scan_shift']=.2
    with pytest.raises(ValueError,match='another production'):validate(leaf,Path.cwd())


def test_run_sources_preserves_staged_files_shared_cap_and_waits_query_outage(case,tmp_path,monkeypatch):
    from kinbot import run_policy
    import subprocess
    p,r,c=case
    source=auto._base(c,p,'character',c['geometry'],mult=3)
    source.update(orbital_character_only=True,anchors=[dict(name='far',geometry_angstrom=c['geometry'])])
    shared=tmp_path/'.kinbot_resource_policy.json';shared.write_text('{}')
    monkeypatch.setattr(run_policy,'find',lambda path:shared)
    calls=[]
    def run(paths,**kwargs):
        calls.append(kwargs)
        if len(calls)==1:raise subprocess.CalledProcessError(1,['squeue','-h'])
    monkeypatch.setattr(recon,'run_many',run)
    def audit(path):
        s=json.loads(path.read_text())
        for j in s['jobs']:j['status']='complete'
        return s
    monkeypatch.setattr(recon,'audit',audit)
    waits=[];monkeypatch.setattr(auto.time,'sleep',waits.append)
    auto._run_sources(tmp_path/'phase',{'character':source},p)
    assert waits==[20.] and all(x['max_jobs']==10 and x['allocation_retries']==2 for x in calls)
    path=tmp_path/'phase/character/reconnaissance.json'
    saved=json.loads(path.read_text());assert saved['shared_policy_file']==str(shared)
    files={f:(f.read_bytes(),f.stat().st_mtime_ns) for f in path.parent.rglob('*') if f.is_file() and f!=path}
    auto._run_sources(tmp_path/'phase',{'character':source},p)
    assert all((f.read_bytes(),f.stat().st_mtime_ns)==v for f,v in files.items())
    source['charge']=2
    with pytest.raises(auto.ReviewRequired,match='Immutable'):auto._run_sources(tmp_path/'phase',{'character':source},p)


def test_partial_native_completion_cannot_be_promoted(case,tmp_path,monkeypatch):
    p,r,c=case
    source=auto._base(c,p,'character',c['geometry'],mult=3)
    source.update(orbital_character_only=True,anchors=[dict(name='far',geometry_angstrom=c['geometry'])])
    monkeypatch.setattr(recon,'run_many',lambda *a,**k:None)
    monkeypatch.setattr(recon,'audit',lambda path:dict(jobs=[dict(status='incomplete')]))
    with pytest.raises(auto.ReviewRequired,match='incomplete'):auto._run_sources(tmp_path/'phase',{'character':source},p)


@pytest.fixture
def scientific_panel(tmp_path,monkeypatch):
    # Analytic six-orbital toy Hamiltonian, two neutral three-electron fragments.
    # This exercises the assessor and native CI parser, not molecular accuracy.
    bases=['avdz-f12','avtz-f12','avqz-f12']
    c=dict(policy=deepcopy(auto.POLICY),symbols=['Li','Li'],electrons=6,mappings=[[0],[1]],bases=bases,
        full_space=dict(active_electrons=4,active_orbitals=4),
        sampling=dict(basis='avdz-f12'),spin=dict(reference=dict(basis=bases[1:],extrapolation_power=3.7),
            gap_levels=[dict(basis=bases[1:],extrapolation_power=3.7)]))
    points={k:dict(distance_angstrom=d,geometry_angstrom=[[0.,0.,0.],[d,0.,0.]])
            for k,d in [('r0',30.),('r1',20.),('r2',3.),('a0',3.),('back',3.)]}
    pairs=[('r2','back')];angles={'a0':'r2'}
    # Closed fragment orbitals precede the symmetric/antisymmetric SOMO pair.
    mat=np.eye(6)[[0,3,1,4,2,5]];mat[2:4]=np.array([[1,1],[1,-1]])@mat[2:4]/np.sqrt(2)
    atom_basis=[[(0,(1.,),((1.,),(1.,),(1.,)))]]*2
    def orbitals_at(point):
        return ev.Orbitals(mat.copy(),np.array([2.,2.,1.,1.,0.,0.]),c['symbols'],
            np.array(points[point]['geometry_angstrom']),deepcopy(atom_basis),[[0,1,2],[3,4,5]])
    jobs={};all_jobs={};scf={}
    for name in ('compact1','compact3','full1','full3'):
        jobs[name]=[];mult=int(name[-1]);nactive=4 if name.startswith('full') else 2
        for basis in bases:
            for point in points:
                ident=name+'_'+basis+'_'+point;d=tmp_path/ident;d.mkdir()
                ci=' 20 .7071067811865475\n 02 -.7071067811865475\n' if mult==1 else ' aa 1.0\n'
                text=(' PROGRAM * MULTI\n'+f'Number of closed-shell orbitals: {(6-nactive)//2} ({(6-nactive)//2})\n'
                      +f'Number of active orbitals: {nactive} ({nactive})\n'
                      +f'Number of active electrons: {nactive} Spin symmetry={"Singlet" if mult==1 else "Triplet"} Space symmetry=1\n'
                      +'Number of states: 1\nNatural orbital dump at molpro section 2140.2 (Orbital set 2)\n'
                      +'Transform the CI vector\nCI Coefficients of symmetry 1\n'+ci+' Energy: -15\n'
                      +' PROGRAM * RSPT2\n T1 diagnostic: 0.001\n')
                out=d/(ident+'.out');out.write_text(text)
                j=dict(name=ident,basis=basis,anchor=point,directory=str(d),output_sha256=file_sha256(out),
                       energy_hartree=-15.,probe_active_space=dict(active_electrons=nactive,active_orbitals=nactive),
                       correlated_energies_hartree=[dict(energy_hartree=-15.) for _ in range(4)],
                       orbitals=orbitals_at(point))
                jobs[name].append(j);all_jobs[(name,basis,point)]=j
    for i in range(2):
        scf[f'fragment{i}']=[]
        for basis in bases:
            o=ev.Orbitals(np.eye(3),np.array([2.,1.,0.]),['Li'],np.array([points['r0']['geometry_angstrom'][i]]),
                          [deepcopy(atom_basis[i])],[[0,1,2]])
            scf[f'fragment{i}'].append(dict(basis=basis,energy_hartree=-7.5,orbitals=o))
    monkeypatch.setattr(auto,'_xml',lambda job,kind='natural':job['orbitals'])
    monkeypatch.setattr(auto,'_files',lambda jobs:{'synthetic':'synthetic'})
    from kinbot import molpro_reference_identity
    monkeypatch.setattr(molpro_reference_identity,'report',lambda *a,**k:dict(
        rohf_somo_active_projection_eigenvalues=[1.,1.],closest_density_root=1,
        reference_states=[dict(root=1,density_distance=0.)]))
    return c,points,angles,pairs,jobs,scf,all_jobs


def test_scientific_assessor_accepts_consistent_analytic_evidence(scientific_panel):
    c,points,angles,pairs,jobs,scf,_=scientific_panel
    report=auto.assess(c,points,angles,pairs,jobs,scf)
    assert report['status']=='accepted', [x for x in report['checks'] if not x['passed']]
    assert 'CBS/full' in report['details']
    assert any(x['test'].endswith('/reference_return_energy') for x in report['checks'])
    assert 'sampling/corrected_angular_residual' in report['details']


@pytest.mark.parametrize('fault',['energy','angular','cbs','ionic','missing_ci','wrong_partition','t1','reference_density','plateau','missing_basis','wrong_orbitals'])
def test_scientific_assessor_rejects_concrete_failure_modes(scientific_panel,monkeypatch,fault):
    c,points,angles,pairs,jobs,scf,index=scientific_panel
    job=index[('compact1','avtz-f12','r2')]
    if fault=='energy':job['correlated_energies_hartree'][0]['energy_hartree']+=.1
    if fault=='angular':index[('compact1','avdz-f12','a0')]['correlated_energies_hartree'][0]['energy_hartree']+=.01
    if fault=='cbs':
        # Opposite sub-threshold errors at T/Q become super-threshold at CBS.
        for basis,sign in [('avtz-f12',-1),('avqz-f12',1)]:
            for name in ('compact1','compact3'):
                for k in ('r2','back','a0'):
                    index[(name,basis,k)]['correlated_energies_hartree'][0]['energy_hartree']+=sign*.4/auto.AUtoKCAL
    if fault in ('ionic','missing_ci','wrong_partition','t1'):
        job=index[('compact1','avtz-f12','r0')] if fault=='ionic' else job
        job=index[('compact3','avtz-f12','r2')] if fault=='t1' else job
        out=Path(job['directory'])/(job['name']+'.out');text=out.read_text()
        if fault=='ionic':text=text.replace('02 -.7071067811865475','02 .7071067811865475')
        if fault=='missing_ci':text=text.replace('Transform the CI vector','omitted')
        if fault=='wrong_partition':text=text.replace('Number of active orbitals: 2','Number of active orbitals: 3')
        if fault=='t1':text=text.replace('T1 diagnostic: 0.001','T1 diagnostic: 0.1')
        out.write_text(text);job['output_sha256']=file_sha256(out)
    if fault=='reference_density':
        from kinbot import molpro_reference_identity
        monkeypatch.setattr(molpro_reference_identity,'report',lambda *a,**k:dict(
            rohf_somo_active_projection_eigenvalues=[1.,1.],closest_density_root=2,
            reference_states=[dict(root=2,density_distance=0.)]))
    if fault=='plateau':
        for name in ('compact1','compact3'):index[(name,'avtz-f12','r1')]['correlated_energies_hartree'][0]['energy_hartree']+=.001
    if fault=='missing_basis':jobs['full1']=[j for j in jobs['full1'] if j['basis']!='avqz-f12']
    if fault=='wrong_orbitals':job['orbitals'].coefficients[[0,2]]=job['orbitals'].coefficients[[2,0]]
    if fault in ('missing_ci','wrong_partition','missing_basis'):
        with pytest.raises((ValueError,KeyError)):auto.assess(c,points,angles,pairs,jobs,scf)
    else:
        report=auto.assess(c,points,angles,pairs,jobs,scf)
        assert report['status']=='needs_review'
        failed=[x['test'] for x in report['checks'] if not x['passed']]
        assert failed
        if fault=='cbs':assert any(v.startswith('CBS/') for v in failed)


def test_complete_wavefunction_overlap_includes_inactive_space():
    ci={'20':2**-.5,'02':-2**-.5}
    t=np.eye(6)
    assert ev.wavefunction_overlap(ci,ci,t,2,1)==pytest.approx(1.)
    angle=.31
    t[:2,:2]=[[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
    assert ev.wavefunction_overlap(ci,ci,t,2,1)==pytest.approx(1.)
    # Same active pair is not sufficient when an inactive orbital is replaced.
    t=np.eye(6);t[[0,5]]=t[[5,0]]
    assert ev.ci_overlap(ci,ci,t[2:4,2:4],1)==pytest.approx(1.)
    assert ev.wavefunction_overlap(ci,ci,t,2,1)==pytest.approx(0.)
    # A rotation involving an active and closed orbital uses full determinants.
    t=np.eye(6);t[np.ix_([0,2],[0,2])]=[[0.,-1.],[1.,0.]]
    assert ev.wavefunction_overlap({'20':1.},{'20':1.},t,2,1)==pytest.approx(1.)


def test_pool_attachment_recovers_preparation_crash_but_never_reassigns_live_work(case,tmp_path,monkeypatch):
    from kinbot import run_policy
    p,r,c=case
    source=auto._base(c,p,'character',c['geometry'],mult=3)
    source.update(orbital_character_only=True,anchors=[dict(name='far',geometry_angstrom=c['geometry'])])
    root=tmp_path/'phase';root.mkdir();path=atomic_json(root/'character.json',source)
    manifest=recon.prepare(path,root/'character')
    shared=tmp_path/'pool.json';shared.write_text('{}')
    monkeypatch.setattr(run_policy,'find',lambda path:shared)
    monkeypatch.setattr(recon,'run_many',lambda *a,**k:None)
    monkeypatch.setattr(recon,'audit',lambda path:dict(jobs=[dict(status='complete')]))
    auto._run_sources(root,{'character':source},p)
    saved=json.loads(manifest.read_text());assert saved['shared_policy_file']==str(shared)
    saved.pop('shared_policy_file');saved['jobs'][0]['job_id']='1234'
    atomic_json(manifest,saved)
    with pytest.raises(auto.ReviewRequired,match='resource-pool'):
        auto._run_sources(root,{'character':source},p)
    assert json.loads(manifest.read_text())==saved


def test_relaxation_uses_production_constraints_and_reuses_complete_native_tasks(case,tmp_path,monkeypatch):
    from ase import Atoms
    from ase.io import write
    from kinbot import vrc_graph,vrc_tst_scan
    from kinbot.anl import dispatch
    p,r,c=case;monkeypatch.chdir(tmp_path)
    jobs=selected_jobs(tmp_path,c);pilots={'pilot1':jobs['compact1']}
    seen=[];created=[]
    def select(base,name,spec):
        seen.append(deepcopy(spec));d=(base/name).resolve();td=d/'tasks/relaxed_scan'
        if not td.exists():
            td.mkdir(parents=True);created.append(d)
            molecule=spec['molecule'];task=spec['tasks'][0]
            write(td/'final.xyz',Atoms(molecule['symbols'],positions=molecule['positions']))
            model=deepcopy(task['profile']['calculator_kwargs']['state_model'])
            checkpoint=td/'native.wfu';checkpoint.write_bytes(b'synthetic returned orbital checkpoint')
            model['orbital_selection'].update(path=str(checkpoint),sha256=file_sha256(checkpoint))
            atomic_json(td/'execution.json',dict(status='executed',details=dict(
                energy_ev=-100.,orbital_checkpoint_model=model)))
        return d
    monkeypatch.setattr(vrc_graph,'select_graph',select)
    monkeypatch.setattr(dispatch,'preflight_until_scheduler_ready',lambda *a,**k:None)
    monkeypatch.setattr(dispatch,'_load',lambda d:(None,None,dict(tasks=dict(relaxed_scan=dict(status='complete')))))
    def unexpected(*a,**k):raise AssertionError('A completed diagnostic task was resubmitted')
    monkeypatch.setattr(vrc_tst_scan,'_drive_vrc_dispatch',unexpected)
    root=tmp_path/'review';root.mkdir()
    first=auto.relax(root,c,p,r,pilots);count=len(created)
    assert len(first)==len(c['radial'])==count
    for spec in seen:
        task=spec['tasks'][0]
        assert task['optimizer']['constraints']['fixed_bonds']==[c['bond']]
        assert task['profile']['method']=='caspt2-f12(2,2)'
        assert task['profile']['calculator_kwargs']['state_model']['state_average']==[
            dict(multiplicity=1,roots=1,weights=[1.])]
        assert spec['molecule']['multiplicity']==1
    second=auto.relax(root,c,p,r,pilots)
    assert len(created)==count and second==first
