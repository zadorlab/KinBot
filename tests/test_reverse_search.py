"""Reverse-search staging, both fragment orders, and native worker routing."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
import pytest

from kinbot import reverse_search
from kinbot.execution import atomic_json, file_sha256
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.reactions.reac_abstraction import abstraction_align, Abstraction


def test_abstraction_alignment_preserves_atoms_in_both_fragment_orders():
    # OH abstracts H from methane. Exchange fragment storage order, including
    # an abstractor whose index is neither zero nor the first stored fragment.
    atoms=np.array(['O','H','C','H','H','H','H'])
    geom=np.array([[0,0,0],[0,0,.97],[8,0,0],[9.09,0,0],
                   [7.637,1.028,0],[7.637,-.514,.89],[7.637,-.514,-.89]])
    a=abstraction_align(geom,[0,3,2],atoms,2)
    order=np.array([2,3,4,5,6,0,1]);inverse=np.argsort(order)
    b=abstraction_align(geom[order],[5,1,0],atoms[order],5)
    for left,right in zip(a,b):
        assert np.allclose(left,right[inverse])
        for fragment in ([0,1],[2,4,5,6]):
            for i in fragment:
                for j in fragment:
                    assert np.linalg.norm(left[i]-left[j])==pytest.approx(np.linalg.norm(geom[i]-geom[j]))
    assert np.linalg.norm(a[1][0]-a[1][3])==pytest.approx(.97)


def test_bimolecular_atom_placement_never_overlaps_fragments():
    a=StationaryPoint('a',0,2,atom=['H'],geom=np.zeros((1,3)))
    b=StationaryPoint('b',0,2,atom=['H'],geom=np.zeros((1,3)))
    for p in (a,b):p.characterize()
    structure=[['H',0.,0.,0.],['H',0.,0.,0.]]
    original=deepcopy(structure)
    pair=StationaryPoint('pair',[0,0],1,structure=deepcopy(structure),fragA=a,fragB=b)
    assert np.linalg.norm(pair.geom[0]-pair.geom[1])>=4.
    assert structure==original


@pytest.fixture
def staged(tmp_path,monkeypatch):
    case=tmp_path/'case';case.mkdir()
    settings=dict(charge=[0,0],mult=[2,1,2],bimol=1)
    data=dict(schema=1,case=str(case),entries=[dict(id='pair',source=str(case),
        products=['p','q'],reactions=['rxn'],settings=settings)])
    monkeypatch.setattr(reverse_search,'plan',lambda p:deepcopy(data))
    copies=[]
    def seed(source,target,products):
        copies.append(str(target))
        (target/'kinbot.db').write_bytes(b'accepted native data')
        atomic_json(target/'fragment_seed.json',dict(products=products))
    monkeypatch.setattr(reverse_search,'_seed_fragment_qc',seed)
    return NS(case=case,data=data,copies=copies)


def test_staging_resumes_atomic_rename_and_never_overwrites_native_work(staged,monkeypatch):
    t=staged;original=reverse_search.atomic_json;interrupted=[]
    def write(path,data):
        if Path(path).name=='manifest.json' and data['entries'][0].get('input') and not interrupted:
            interrupted.append(True);raise OSError('interrupt after rename')
        return original(path,data)
    monkeypatch.setattr(reverse_search,'atomic_json',write)
    with pytest.raises(OSError):reverse_search.prepare(t.case)
    result=reverse_search.prepare(t.case)
    assert len(t.copies)==1
    inp=reverse_search.inputs(t.case)[0]
    assert inp.parent.joinpath('kinbot.db').read_bytes()==b'accepted native data'
    assert reverse_search.prepare(t.case)==result
    assert len(t.copies)==1
    atomic_json(inp,{'changed':True})
    with pytest.raises(ValueError,match='changed'):reverse_search.prepare(t.case)


def test_launch_is_idempotent_and_adopts_unreceipted_live_process(staged,monkeypatch):
    t=staged;reverse_search.prepare(t.case)
    inp=reverse_search.inputs(t.case)[0];calls=[];live=[dict(pid=123)]
    monkeypatch.setattr('kinbot.run_cleanup._processes',lambda p:live)
    def launch(*a,**kw):calls.append((a,kw));return NS(pid=321)
    monkeypatch.setattr(reverse_search.subprocess,'Popen',launch)
    reverse_search.launch(t.case,'/code')
    assert not calls
    live.clear()
    reverse_search.launch(t.case,'/code');reverse_search.launch(t.case,'/code')
    assert len(calls)==1
    assert calls[0][0][0][-2:]==['kinbot.kb',str(inp)]
    assert calls[0][1]['env']['PYTHONPATH']=='/code'
    with pytest.raises(RuntimeError,match='stale'):reverse_search.ensure_integrated(t.case)


def test_sella_abstraction_renders_saddle_worker_and_preserves_live_script(tmp_path,monkeypatch):
    from kinbot.qc import QuantumChemistry
    from kinbot.reac_family import carry_out_reaction
    monkeypatch.chdir(tmp_path)
    atomic_json('input.json',dict(bimol=1,charge=[0,0],mult=[2,1,2],
        structure=[['O',0.,0.,0.,'H',0.,0.,.97],
                   ['H',0.,0.,0.,'H',0.,0.,.74]],use_sella=1,
        method='B3LYP',basis='cc-pVDZ',barrier_threshold=150,
        queuing='local',calc_kwargs={'EmpiricalDispersion':'GD3BJ'}))
    par=Parameters('input.json').par
    a=StationaryPoint('a',0,2,structure=par['structure'][0]);a.characterize()
    b=StationaryPoint('b',0,1,structure=par['structure'][1]);b.characterize()
    pair=StationaryPoint('pair',[0,0],2,structure=deepcopy(par['structure']),fragA=a,fragB=b)
    pair.characterize()
    qc=QuantumChemistry(par)
    monkeypatch.setattr(qc,'check_qc',lambda n:'normal')
    monkeypatch.setattr(qc,'get_qc_geom',lambda *a,**kw:(0,pair.geom))
    monkeypatch.setattr(qc,'is_in_database',lambda n:False)
    submissions=[]
    monkeypatch.setattr(qc,'submit_qc',lambda *a,**kw:submissions.append((a,kw)) or 1)
    obj=Abstraction(pair,qc,par,[0,2,3],'pair_abstraction_1_3_4')
    step=carry_out_reaction(obj,1,'g16',bimol=1)
    script=Path(obj.instance_name+'.py').read_text()
    assert step==3 and 'order=1' in script and 'order=0' not in script
    assert 'B3LYP' in script and 'GD3BJ' in script
    compile(script,'worker','exec')
    saved=Path(obj.instance_name+'.py').read_bytes()
    monkeypatch.setattr(qc,'check_qc',lambda n:'running')
    assert carry_out_reaction(obj,2,'g16',bimol=1)==2
    assert Path(obj.instance_name+'.py').read_bytes()==saved and len(submissions)==1


def test_entrance_characterization_polls_same_jobs_and_requires_native_success(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    point=dict(energy_hartree=-10.,zpe_hartree=.01,frequencies=[20.])
    endpoint=NS(name='irc',atom=['H'],geom=[[0,0,0]],charge=0,mult=2)
    checks=iter([dict(status='waiting'),dict(status='complete',point=point)])
    monkeypatch.setattr('kinbot.l1_endpoint.ensure',lambda *a,**k:next(checks))
    monkeypatch.setattr(reverse_search.time,'sleep',lambda x:None)
    monkeypatch.setattr(reverse_search,'routing_name',lambda p:p)
    monkeypatch.setattr(reverse_search.l1_screening,'theory',lambda p:{})
    monkeypatch.setattr(reverse_search.l1_screening,'read_points',lambda *a,**k:[
        dict(energy_hartree=-4.,zpe_hartree=0.),dict(energy_hartree=-5.98,zpe_hartree=0.)])
    opt=NS(shigh=-1,species=endpoint)
    calls=[]
    def optimize():
        calls.append(1);opt.shigh=1 if len(calls)==2 else 0
    opt.do_optimization=optimize
    monkeypatch.setattr('kinbot.optimize.Optimize',lambda *a,**k:opt)
    result=reverse_search.characterize_entrance(endpoint,['p','q'],{'vdW_detection':.5},None)
    assert result['status']=='characterized' and len(calls)==2
    assert result['l2']['energy_hartree']==-10.


def test_terminal_entrance_failure_never_retries_or_invents_well(monkeypatch):
    monkeypatch.setattr('kinbot.l1_endpoint.ensure',lambda *a,**k:dict(status='unresolved',reason='failed native Hessian'))
    result=reverse_search.characterize_entrance(None,[],{},None)
    assert result['status']=='unresolved' and 'l2' not in result



def test_native_irc_assignment_and_stereopath_for_reverse_abstraction():
    from kinbot.irc import IRC
    from kinbot.reaction_path import prepare_stereopath, reaction_path_id
    # OH + H2 -> H2O + H, with common atom order O,H,H,H.
    a=StationaryPoint('a',0,2,atom=['O','H'],geom=np.array([[0.,0.,0.],[0.,.97,0.]]));a.characterize()
    b=StationaryPoint('b',0,1,atom=['H','H'],geom=np.array([[0.,0.,0.],[.74,0.,0.]]));b.characterize()
    incoming=np.array([[0.,0.,0.],[0.,.97,0.],[4.,0.,0.],[4.74,0.,0.]])
    outgoing=np.array([[0.,0.,0.],[0.,.97,0.],[.97,0.,0.],[5.,0.,0.]])
    pair=StationaryPoint('pair',0,2,atom=['O','H','H','H'],geom=incoming,fragA=a,fragB=b)
    pair.characterize()
    def geometry(name,*args,**kw):return 0,(incoming if '_R_' in name else outgoing).copy()
    obj=NS(instance_name='pair_abstraction_1_3_4',species=pair,qc=NS(get_qc_geom=geometry))
    par={'bimol':1,'electronic_state_models':[],'optical_population':'specified'}
    endpoint=IRC(obj,par).irc2stationary_pt()
    assert endpoint.name.endswith('_F_prod')
    assert obj.irc_reactant.name.endswith('_R_prod')
    ts=StationaryPoint('ts',0,2,atom=pair.atom,geom=(incoming+outgoing)/2,wellorts=1)
    ts.characterize()
    prepare_stereopath(ts,pair,endpoint)
    obj.ts=ts
    assert reaction_path_id(obj)
    # Return to the reactants at both ends is not an abstraction.
    obj.qc.get_qc_geom=lambda *args,**kw:(0,incoming.copy())
    assert IRC(obj,par).irc2stationary_pt()==0


@pytest.mark.parametrize('swapped', [False, True])
def test_export_keeps_native_jobs_through_single_fragment_renaming_and_l2(tmp_path, monkeypatch, swapped):
    from kinbot.irc import IRC
    monkeypatch.chdir(tmp_path)
    # Synthetic endpoint connectivity exercises the object-alias lifecycle;
    # it is not evidence of a physical reaction or a converged QC calculation.
    a=StationaryPoint('a',0,2,atom=['O','H'],geom=np.array([[0.,0.,0.],[0.,.97,0.]]))
    b=StationaryPoint('b',0,1,atom=['H','H'],geom=np.array([[0.,0.,0.],[.74,0.,0.]]))
    for p in (a,b): p.characterize(); p.energy=-1.; p.zpe=0.; p.freq=[500.]
    incoming=np.array([[0.,0.,0.],[0.,.97,0.],[4.,0.,0.],[4.74,0.,0.]])
    bound=np.array([[0.,0.,0.],[0.,.97,0.],[.97,0.,0.],[0.,0.,.97]])
    pair=StationaryPoint('pair',0,2,atom=['O','H','H','H'],geom=incoming,fragA=a,fragB=b)
    pair.characterize()
    direction='F' if swapped else 'R'
    def geometry(job,*args,**kwargs):
        return 0,(incoming if f'_IRC_{direction}_' in job else bound).copy()
    obj=NS(instance_name='pair_abstraction',species=pair,qc=NS(get_qc_geom=geometry))
    par=dict(bimol=1,electronic_state_models=[],optical_population='specified')
    obj.irc_prod=IRC(obj,par).irc2stationary_pt()
    jobs=dict(obj.irc_endpoint_jobs)
    assert jobs['reactant'].endswith(f'_IRC_{direction}_prod')
    assert jobs['product'] == obj.irc_prod.name
    fragments,_=obj.irc_prod.start_multi_molecular()
    assert fragments == [obj.irc_prod] and obj.irc_prod.name != jobs['product']
    # Actual product characterization reorders/optimizes this very object.
    obj.irc_prod.geom += 20.
    obj.irc_prod.energy=-8.; obj.irc_prod.zpe=.01; obj.irc_prod.freq=[300.]
    obj.ts_opt=NS(species=obj.irc_prod)
    obj.prod_opt=[NS(species=obj.irc_prod)]
    obj.do_vdW=False
    pair.reac_obj=[obj]; pair.reac_ts_done=[-1]
    queried=[]
    def read(directory,names,*args,**kwargs):
        queried.extend(names)
        return [dict(name=n,symbols=list(pair.atom),positions=geometry(n)[1].tolist(),
                     energy_hartree=-9.,zpe_hartree=None,frequencies=[]) for n in names]
    monkeypatch.setattr(reverse_search.l1_screening,'read_points',read)
    monkeypatch.setattr(reverse_search.l1_screening,'theory',lambda p:{})
    monkeypatch.setattr('kinbot.reaction_path.reaction_path_id',lambda o:'tested-separately')
    def characterize(endpoint,*args):
        assert endpoint.name == jobs['reactant'] and np.allclose(endpoint.geom,incoming)
        endpoint.geom += 30.  # Optimize's in-place mutation must not poison next export.
        return dict(status='characterized',l2=NS().__dict__)
    monkeypatch.setattr(reverse_search,'characterize_entrance',characterize)
    Path('kinbot.db').write_bytes(b'fixture')
    for _ in range(2):
        reverse_search.export_result(pair,par,None)
        r=json.loads(Path('reverse_search_result.json').read_text())
        row=r['reactions'][0]
        assert row['product_l1_job'] == jobs['product']
        assert row['product_endpoint_l1']['positions'] == bound.tolist()
        assert row['entrance_l1']['positions'] == incoming.tolist()
        assert not r['mess_ready']
    assert queried == [jobs['reactant'],jobs['product']]*2
