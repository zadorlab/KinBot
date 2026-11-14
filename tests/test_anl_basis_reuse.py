"""Accepted L1/L2 -> staged ANL -> parsed CBS; no licensed QC or Slurm needed."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from ase.db import connect
from ase.io import read
from ase.units import invcm, Hartree

from kinbot import constants
from kinbot.stationary_pt import StationaryPoint
from kinbot.calculation import load_calculation_record, publish_optimization_result
from kinbot.anl import automatic, anharmonic, basis_reuse as reuse, dispatch, protocol
from kinbot.anl.validation import automatic_composite_spec
from kinbot.anl.recipes import recipe
from kinbot.anl.extrapolation import extrapolation_coefficient
from kinbot.execution import atomic_json, file_sha256
from tests.test_anl_live_protocol import fixed_resources
from tests.test_anl_workflow import _complete_task, _vpt2_output
from tests.test_qc_reuse import worker


def accepted_pair(tmp_path, lower='aug-cc-pVTZ', upper='aug-cc-pVQZ'):
    atoms = molecule('CH4')
    atoms.positions[1] += .000123456789
    point = StationaryPoint.from_ase_atoms(atoms)
    point.characterize()
    point.wellorts = 0
    db = connect(tmp_path/'kinbot.db')
    def profile(basis):
        return SimpleNamespace(calculator='gaussian', method='B2PLYP', basis=basis,
                               calculator_kwargs={'EmpiricalDispersion': 'GD3BJ'})
    qc = SimpleNamespace(db=db, _fingerprint='4'*64,
        _routes={'minimum_tz': {'level':'l1'}, 'minimum_qz': {'level':'l2'}},
        profiles={'l1':profile(lower), 'l2':profile(upper)}, qc='gauss', use_sella=True)
    for name, scale in [('minimum_tz',1.),('minimum_qz',1.002)]:
        a=atoms.copy(); a.positions*=scale
        db.write(a,name=name,data={'status':'normal', 'energy':-40./constants.EVtoHARTREE,
                                   'zpe':.04, 'frequencies':[500.]*9})
    load_calculation_record(point,qc,'minimum_tz')
    reuse.seed(point)
    worker(tmp_path,job='minimum_qz',symbols=list(map(str,point.atom)),positions=point.geom.tolist(),
           extra_kwargs={'method':'B2PLYP','basis':upper,'mult':1,'charge':0})
    with pytest.MonkeyPatch.context() as mp:
        mp.chdir(tmp_path)
        reuse.bind_seed(point,qc,'minimum_qz')
    load_calculation_record(point,qc,'minimum_qz')
    return point,qc


def theory(basis='aug-cc-pVQZ'):
    return dict(calculator='gaussian', method='B2PLYP',basis=basis,dispersion='GD3BJ',
                composite_method='B2PLYP-D3BJ',label='B2PLYP-D3BJ/'+basis)


def spec_from(point,basis='aug-cc-pVQZ',tier='ANL1-F12'):
    record=automatic._l2_species_record(point,theory(basis))
    spec=automatic_composite_spec(record['molecule'],highest_tier=tier,
        vpt2_basis=basis,f12_augmented=True,augmented=True)
    fixed_resources(spec)
    return spec


@pytest.mark.parametrize('bases',[('cc-pVTZ','cc-pVQZ'),('aug-cc-pVTZ','aug-cc-pVQZ')])
@pytest.mark.parametrize('tier',['ANL1','ANL1-F12'])
def test_accepted_pair_no_duplicate_b2plyp_optimizations_or_harmonics(tmp_path,bases,tier):
    point,qc=accepted_pair(tmp_path,*bases)
    spec=spec_from(point,bases[1],tier)
    tasks={t['id']:t for t in spec['tasks']}
    assert not any(t.get('profile',{}).get('method')=='B2PLYP' for t in tasks.values())
    assert not any(t.get('result_parser',{}).get('kind')=='gaussian_harmonic' for t in tasks.values())
    pair=spec['intent']['anharmonic_basis_extension']
    assert pair['bases']==list(bases)
    assert pair['upper_task']=='gaussian_vpt2'
    assert tasks['gaussian_vpt2']['result_parser']['basis']==bases[1]
    tz=tasks[pair['lower_task']]
    assert '/'+bases[0] in tz['input_template']
    assert tz['geometry_record']['positions']==qc.db.get(name='minimum_tz').positions.tolist()
    assert anharmonic.extend(spec)==spec
    dispatch.validate_spec(spec)
    reduced=protocol.replacement(spec,{'anharmonic_qz':False})
    assert pair['lower_task'] in protocol.omitted(reduced)
    assert pair['upper_task'] not in protocol.omitted(reduced)


def test_published_l2_retains_pair_on_restart_and_reference_export(tmp_path):
    point,qc=accepted_pair(tmp_path)
    def publish(row,target):
        qc.db.write(row.toatoms(),name=target,data=dict(row.data))
    qc.publish_result=publish
    opt=SimpleNamespace(species=point,qc=qc,shigh=1,shir=1,name='minimum',
        par={'high_level':True},log_name=lambda high:'published_high')
    publish_optimization_result(opt)
    saved=qc.db.get(name='published_high')
    assert saved.data['anl_basis_reuse']['points'][0]['source_job']=='minimum_tz'
    clone=StationaryPoint.from_ase_atoms(saved.toatoms()); clone.characterize(); clone.wellorts=0
    qc._routes['published_high']={'level':'l2'}
    load_calculation_record(clone,qc,'published_high')
    assert reuse.exported(clone)==reuse.exported(point)
    from kinbot.mess import MESS
    # Exercise the real reference model writer while isolating unrelated HIR/symmetry work.
    writer=MESS.__new__(MESS);writer.par={}
    writer._stable_frequencies=lambda p:p.freq
    writer.make_rotors=lambda p,e:''
    writer._thermal_geometry=lambda p:p.geom
    writer._parent_symmetry=lambda p:1.
    clone.sigma_ext=clone.nopt=1
    model=writer.thermochemistry_model(clone)
    path=tmp_path/'models.json';path.write_text(json.dumps({'schema':1,'species':{clone.smiles:model}}))
    imported=automatic._model_record(path,{'smiles':clone.smiles,'charge':0,'multiplicity':1},theory())
    assert imported['molecule']['source']['anl_basis_reuse']==reuse.exported(point)


@pytest.mark.parametrize('problem',['hir_changed','missing','unverified_route','wrong_family'])
def test_no_speculative_reuse_of_unrelated_or_unverified_minimum(tmp_path,problem):
    point,qc=accepted_pair(tmp_path)
    if problem=='hir_changed':
        point.geom=point.geom.copy()
        point.geom[0,0]+=.1
    elif problem=='missing':
        del point.anl_basis_reuse
    elif problem=='unverified_route':
        qc._routes['minimum_tz']['source_profile_sha256']='5'*64
        load_calculation_record(point,qc,'minimum_tz');reuse.seed(point)
        load_calculation_record(point,qc,'minimum_qz')
    else:
        point,qc=accepted_pair(tmp_path, 'cc-pVTZ','aug-cc-pVQZ')
    spec=spec_from(point)
    tasks={t['id']:t for t in spec['tasks']}
    assert tasks['l3_vpt2_tz_geometry']['profile']['basis']=='aug-cc-pVTZ'
    assert 'l3_vpt2_qz_geometry' not in tasks
    assert 'geometry_record' not in tasks['gaussian_vpt2_tz']
    assert tasks['gaussian_vpt2']['geometry_from']=='initial'


@pytest.mark.parametrize('mutation',['hash','charge','multiplicity','symbols','masses','position'])
def test_invalid_record_is_not_silently_used(tmp_path,mutation):
    p,_=accepted_pair(tmp_path);s=spec_from(p)
    task=next(t for t in s['tasks'] if t['id']=='gaussian_vpt2_tz')
    point=task['geometry_record']
    if mutation=='hash':point['energy_hartree']+=1
    elif mutation=='charge':point['charge']=1;reuse.seal(point)
    elif mutation=='multiplicity':point['multiplicity']=3;reuse.seal(point)
    elif mutation=='symbols':point['symbols'][0]='N';reuse.seal(point)
    elif mutation=='masses':point['masses'][1]=2.;reuse.seal(point)
    else:point['positions'][0][0]+=.1;reuse.seal(point)
    with pytest.raises(ValueError):dispatch.validate_spec(s)


def test_atom_remapping_keeps_both_minima_and_does_not_invent_correspondence(tmp_path):
    p,qc=accepted_pair(tmp_path);record=reuse.exported(p)
    perm=[1,0,2,3,4]
    symbols=[str(p.atom[i]) for i in perm];positions=p.geom[perm]
    mapped=reuse.remap(record,perm,symbols,positions)
    mol=dict(symbols=symbols,positions=positions.tolist(),charge=0,multiplicity=1,
             source={'geometry_sha256':reuse.geometry_hash(symbols,positions),'anl_basis_reuse':mapped})
    selected=reuse.select(mol,'aug-cc-pVTZ')
    assert selected['positions']==qc.db.get(name='minimum_tz').positions[perm].tolist()
    assert selected['atom_permutation']==perm
    assert record==reuse.exported(p)


def test_full_staging_publication_race_resume_and_native_scalar_assembly(tmp_path,monkeypatch):
    monkeypatch.setattr(dispatch,'assign_partitions',fixed_resources)
    p,_=accepted_pair(tmp_path);spec=spec_from(p)
    pair=spec['intent']['anharmonic_basis_extension']
    spec['tasks']=[t for t in spec['tasks'] if t['id'] in pair['required_tasks']]
    graph=dispatch.prepare(spec,tmp_path/'graph')
    tasks={t['id']:t for t in spec['tasks']}
    for name in pair['required_tasks']:
        task=tasks[name]
        expected=task.get('geometry_record',spec['molecule'])['positions']
        assert np.allclose(read(graph/'tasks'/name/'geometry.xyz').positions,expected,rtol=0,atol=5.1e-9)
    state=dispatch._load(graph)[2]
    for n in tasks:state['tasks'][n].update(status='submitted',job_id='123')
    atomic_json(graph/'state.json',state)
    before={f:file_sha256(f) for f in (graph/'tasks').rglob('*') if f.is_file()}
    monkeypatch.setattr(dispatch,'_job_active',lambda job:True)
    for _ in range(2):
        current=dispatch.advance(graph,submit=False)
        assert all(e['status']=='submitted' for e in current['tasks'].values())
    assert before=={f:file_sha256(f) for f in before}
    assert not anharmonic.available(spec,current)
    _complete_task(graph,tasks[pair['lower_task']],_vpt2_output(pair['bases'][0],-137.))
    assert not anharmonic.available(spec,dispatch._load(graph)[2])
    _complete_task(graph,tasks[pair['upper_task']],_vpt2_output(pair['bases'][1],-136.))
    result=anharmonic.component(graph,state_id='test')
    assert result.value_hartree==pytest.approx((-136.+extrapolation_coefficient(4,power=3.7))*invcm/Hartree)
    assert result.basis=='CBS(aug-cc-pVTZ,aug-cc-pVQZ)'
    assert 'minimum_tz' in json.dumps(tasks[pair['lower_task']]['geometry_record'])
    for _ in range(2):
        assert anharmonic.component(graph,state_id='test')==result


def test_graph_owner_repeated_prepare_never_rewrites_live_work(tmp_path,monkeypatch):
    monkeypatch.setattr(dispatch,'assign_partitions',fixed_resources)
    p,_=accepted_pair(tmp_path);record=automatic._l2_species_record(p,theory())
    plan={'resources':dict(global_max_jobs=2,partition='test',scf_convergence=1e-8,cc_convergence=1e-7),
          'l2_theory':theory()}
    graph=tmp_path/'graph';entry={'highest_tier':'ANL1-F12'}
    automatic._prepare_stable_graph(plan,entry,record,graph)
    state=dispatch._load(graph)[2]
    state['tasks']['gaussian_vpt2'].update(status='submitted',job_id='live-job')
    atomic_json(graph/'state.json',state)
    before={f:file_sha256(f) for f in graph.rglob('*') if f.is_file()}
    for _ in range(2):automatic._prepare_stable_graph(plan,entry,record,graph)
    assert before=={f:file_sha256(f) for f in before}


def test_explicit_vpt2_refinement_starts_from_reused_tz_not_qz(tmp_path):
    p,_=accepted_pair(tmp_path);s=spec_from(p)
    task=next(t for t in s['tasks'] if t['id']=='gaussian_vpt2_tz')
    from kinbot.anl.vpt2_refinement import _tasks
    geometry,harmonic,vpt2=_tasks(s,task)
    assert geometry['geometry_record']==task['geometry_record']
    assert all('geometry_record' not in t for t in (harmonic,vpt2))
    assert geometry['profile']['basis']=='aug-cc-pVTZ'
    assert vpt2['geometry_from']==geometry['id']


@pytest.mark.parametrize('basis',reuse.BASES)
def test_basis_specific_recipe_labels_and_qz_fallback(basis):
    unpaired=recipe('ANL1-F12',vpt2_method='B2PLYP-D3BJ',vpt2_basis=basis)
    r=next(r for r in unpaired.requirements if r.key=='vpt2_correction')
    assert r.basis==basis
    paired=recipe('ANL1-F12',vpt2_method='B2PLYP-D3BJ',vpt2_cbs=True,vpt2_basis=basis)
    r=next(r for r in paired.requirements if r.key=='vpt2_correction')
    assert r.basis=='CBS('+','.join(reuse.basis_pair(basis))+')'


@pytest.mark.parametrize('problem',['none','different_seed','missing_worker','solvent'])
def test_no_conformer_search_and_existing_worker_lineage(tmp_path,monkeypatch,problem):
    point,qc=accepted_pair(tmp_path)
    monkeypatch.chdir(tmp_path)
    if problem=='solvent':qc.profiles['l1'].calculator_kwargs['scrf']='smd'
    load_calculation_record(point,qc,'minimum_tz')
    point.__dict__.pop('_anl_basis_point',None)
    reuse.seed(point,qc,'minimum_tz')
    if problem=='different_seed':
        worker(tmp_path,job='minimum_qz',symbols=list(map(str,point.atom)),
               positions=(point.geom*1.1).tolist(),extra_kwargs={'mult':1,'charge':0})
    elif problem=='missing_worker':
        (tmp_path/'minimum_qz.py').unlink()
    reuse.bind_seed(point,qc,'minimum_qz')
    load_calculation_record(point,qc,'minimum_qz')
    assert (reuse.exported(point) is not None)==(problem=='none')


@pytest.mark.parametrize('basis',reuse.BASES)
def test_ordinary_input_profile_resolution_preserves_actual_l2_basis(basis):
    resolved=automatic._resolved_l2_theory({'l2_profile':dict(calculator='gaussian',
        method='B2PLYP',basis=basis,calculator_kwargs={'EmpiricalDispersion':'GD3BJ'})})
    assert resolved['basis']==basis
    assert resolved['dispersion']=='GD3BJ'
