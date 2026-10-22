"""Exercise the reported r7 structures without scheduler or native QC calls."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot.execution import atomic_json
from test_requested_orbital_recon import helper, ROOT

PLAN=ROOT/'examples/anl/orbital_review_batches/r7_additional_20261021/batch.json'
CS='471430450000000000002-s06a7b4c221beb748'
SS='791671580700000000002-s36f1e2b924862f12'
CS_GEOM=[[.6324202323251805,-1.1752470042214924,-.46315675743475304],
 [.4882637466760829,.33629426627575826,.33379986596473493],
 [.37161066563934064,1.5610667569524943,.9796107569463828],
 [.10778984408129016,-.30390303514860806,1.5467389510656508],
 [.8633979585525545,1.0318289258603168,-.8499821246258238]]
SS_GEOM=[[1.2000484167532062,-.465541087630126,.101405732222727],
 [.9454445583129932,.7518815002969066,1.3204231600418774],
 [.2070747663828456,-.3354666964627526,2.8547906540460475],
 [-.13482612252117787,1.3245512569621833,.734096658492754],
 [.3346164228619777,-.9793332304359088,-.2861472967430523],
 [2.1495822734806818,-.9709787169152773,.1677215698015446]]


@pytest.fixture
def batch_tree(tmp_path):
    root=tmp_path/'production'; root.mkdir()
    atomic_json(root/'.kinbot_resource_policy.json',dict(schema=1,max_jobs=23,
        excluded_nodes=['badnode'],manage_pending=True,orbital_review_reserve=3))
    for case in ('ch3s','ch3ss'):
        (root/case).mkdir()
        atomic_json(root/case/'input.json',dict(queue_name='test-partition'))
    recipe=json.loads((PLAN.parent/'c_sh3_cs.json').read_text())
    request=root/'ch3s'/CS/'vrctst/orbital_review'/recipe['request_match']['reaction']/'request.json'
    request.parent.mkdir(parents=True)
    atomic_json(request,dict(recipe['request_match'],schema=1,status='needs_review',positions=CS_GEOM))
    well=root/'ch3ss'/SS; well.mkdir()
    structure=[]
    for atom,xyz in zip(['C','S','S','H','H','H'],SS_GEOM):
        structure.extend([atom,*xyz])
    atomic_json(well/(SS+'.json'),dict(structure=structure,charge=0,mult=2,smiles=''))
    (well/'kinbot_monitor.out').write_text('\n'.join('4 0 '+SS+'_hom_sci_'+b for b in ('1_2','2_3')))
    return root,tmp_path/'batch'


def test_three_screens_share_cap_skip_old_review_and_submit_nothing(batch_tree,monkeypatch):
    root,output=batch_tree
    before={p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:pytest.fail('No native/scheduler during preparation'))
    result=helper.batch(root,PLAN,output)
    assert (result['screens'],result['tasks'],result['max_jobs'])==(3,24,3)
    record=helper._saved(output)
    assert len(record['systems'])==3
    assert all(p.read_bytes()==b for p,b in before.items())
    for system in record['systems']:
        config=json.loads(Path(system['configuration']).read_text())
        manifest=json.loads(Path(system['manifest']).read_text())
        assert not config.get('correlation_checks')
        assert all(j['orbital_predecessor'] is None for j in manifest['jobs'])
        for j in manifest['jobs']:
            d=Path(j['directory']); text=(d/(j['name']+'.inp')).read_text().lower()
            assert '{multi,cas-ci;' in text and 'uccsd' not in text and 'rs2' not in text
            script=(d/(j['name']+'.sbatch')).read_text()
            assert '#SBATCH --exclusive' in script and 'partition_max' in script
        if 'ch2shs' in system['title']:
            assert manifest['jobs'][0]['probe_active_space']['active_electrons']==19
            assert manifest['jobs'][0]['probe_active_space']['active_orbitals']==12
            assert helper._saved(Path(system['manifest']).parent.parent)['provisional']
    with pytest.raises(FileExistsError,match='resume'):
        helper.batch(root,PLAN,output)


def test_reviewed_case_channel_is_not_screened_again(batch_tree):
    root,output=batch_tree
    atomic_json(root/'ch3s/input.json',dict(queue_name='test-partition',
        vrc_tst_surface_models={CS:dict(channels=[dict(bond=[1,2])])}))
    result=helper.batch(root,PLAN,output)
    assert result['screens']==2 and result['tasks']==16
    assert result['skipped'][0]['name']=='c_sh3_cs'


def test_provisional_rejected_search_blocks_whole_batch_before_publication(batch_tree):
    root,output=batch_tree
    (root/'ch3ss'/SS/'kinbot_monitor.out').write_text('-999 0 '+SS+'_hom_sci_1_2')
    with pytest.raises(ValueError,match='not an active/accepted search'):
        helper.batch(root,PLAN,output)
    assert not output.exists()


def test_batch_runs_all_manifests_under_current_reserve_and_isolates_failures(batch_tree,monkeypatch):
    root,output=batch_tree
    helper.batch(root,PLAN,output)
    policy=root/'.kinbot_resource_policy.json'; p=json.loads(policy.read_text())
    p.update(orbital_review_reserve=2,excluded_nodes=['badnode','newbad'])
    atomic_json(policy,p)
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:None)
    monkeypatch.setattr('kinbot.run_walltime.partition_walltime',lambda p:'31-00:00:00')
    calls=[]
    def run(paths,**kwargs):
        calls.append((paths,kwargs))
        manifests=[json.loads(p.read_text()) for p in paths]
        assert all(m['slurm_excluded_nodes']==['badnode','newbad'] for m in manifests)
        return manifests
    monkeypatch.setattr(helper.recon,'run_many',run)
    helper.run(output)
    assert len(calls)==1 and len(calls[0][0])==3
    assert calls[0][1]==dict(max_jobs=2,allocation_retries=2,keep_going=True)


def test_provisional_never_installs_model_and_can_bind_only_matching_live_request(batch_tree,monkeypatch):
    root,output=batch_tree
    helper.batch(root,PLAN,output)
    child=output/'ch2shs_ss'
    with pytest.raises(ValueError,match='Bind'):
        helper.review(child,root/'ch3ss/input.json',PLAN,apply=False)
    saved=helper._saved(child); m=json.loads(Path(saved['manifest']).read_text())
    for job in m['jobs']: job['status']='complete'
    atomic_json(saved['manifest'],m)
    original=json.loads((child/'request.json').read_text())
    live=deepcopy(original); live.pop('provisional'); live.pop('note')
    live['positions']=[[x+.01,y,z] for x,y,z in live['positions']]
    request=root/'ch3ss'/SS/'vrctst/orbital_review'/live['reaction']/'request.json'
    request.parent.mkdir(parents=True)
    wrong=deepcopy(live);wrong['products'][1]['multiplicity']=1
    atomic_json(request,wrong)
    with pytest.raises(ValueError,match='fragments differ'):
        helper.bind(child,request)
    assert helper._saved(child)['provisional']
    atomic_json(request,live)
    result=helper.bind(child,request)
    assert not result['production_model_accepted']
    assert helper._saved(child)['binding']['geometry_changed']
    assert not helper._saved(child)['provisional']
    assert helper._saved(output)['reviews']  # parent registry remains resumable
    assert json.loads((child/'provisional_request.json').read_text())==original
