import json
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms
from ase.db import connect
from kinbot.execution import atomic_json
from tests.test_requested_orbital_recon import helper, ROOT

PLAN=ROOT/'examples/anl/orbital_review_batches/r7_progress_20261024'
GEOMS={
'hcs2_h2_exit':[[1.3619841465,.4316796483,-.4363935818],[1.6286344595,2.003052028,.0598243878],[-.2219563325,-.0313960434,-.2589523896],[4.4281356519,1.0649159937,-.8718819614],[4.8025371513,.4911953276,-1.1543606403],[2.1430928239,-.2086023259,-.8299048761]],
'ch2s_sh_exit':[[1.3909722277,-.3188140614,.5235903149],[3.0190817109,2.2002042442,-1.7758151979],[.0520998342,.3731735087,-.0567885485],[1.7643340582,1.84348338,-1.4206797474],[2.3802995507,.0010891264,.2090944803],[1.3470283832,-1.1264959604,1.2471614427]]}


def test_actual_three_request_batch_stages_24_jobs_without_production_writes(tmp_path):
    production=tmp_path/'production';case=production/'ch3ss'
    atomic_json(case/'input.json',dict(queue_name='test-partition'))
    atomic_json(production/'.kinbot_resource_policy.json',dict(schema=1,max_jobs=22,
        excluded_nodes=['node1'],manage_pending=True,orbital_review_reserve=3))
    positions=[[1.8122818351,.0029772478,-.0010645061],[2.2322433887,1.7849627863,-.0036765215],
        [.0180913054,.0086000664,-.0026139346],[1.7933644276,2.0266677332,-1.2480033488],
        [2.2475352842,-.4907718979,-.8623412662],[2.2307468144,-.3955252319,.9183373719]]
    rows=[]
    for name in ('ch2_sh_s_cs','hcs2_h2_exit','ch2s_sh_exit'):
        recipe=json.loads((PLAN/(name+'.json')).read_text());match=recipe['request_match']
        worker=case/match['parent']
        atomic_json(worker/'vrctst/orbital_review'/match['reaction']/'request.json',
                    dict(match,schema=1,status='needs_review',positions=positions))
        if name in GEOMS:
            source=recipe['geometry_source']
            rows.append(dict(name=match['reaction'],product_complex=source['product_complex']))
            connect(worker/'kinbot.db').write(Atoms(match['symbols'],positions=GEOMS[name]),
                name=source['record_name'],data=dict(status='normal'))
    atomic_json(worker/'accepted_l2_network.json',dict(network=dict(reactions=rows)))
    before={p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    result=helper.batch(production,PLAN/'screens.json',tmp_path/'screens')
    assert result['tasks']==24 and result['screens']==result['max_jobs']==3
    saved=helper._saved(tmp_path/'screens')
    for system in saved['systems']:
        manifest=json.loads(Path(system['manifest']).read_text())
        assert len(manifest['jobs'])==8
        assert all(not j.get('job_id') for j in manifest['jobs'])
    assert all(p.read_bytes()==v for p,v in before.items())


@pytest.mark.parametrize('name',list(GEOMS))
def test_product_exit_screen_reads_database_without_writing_and_uses_correct_maps(tmp_path,name):
    recipe=PLAN/(name+'.json');control=json.loads(recipe.read_text());match=control['request_match']
    production=tmp_path/'production';case=production/'ch3ss';worker=case/match['parent']
    request=worker/'vrctst/orbital_review'/match['reaction']/'request.json'
    atomic_json(request,dict(match,schema=1,status='needs_review',positions=[[0,0,0]]*6))
    atomic_json(case/'input.json',dict(queue_name='test-partition'))
    atomic_json(production/'.kinbot_resource_policy.json',dict(schema=1,max_jobs=22,excluded_nodes=[],manage_pending=True,orbital_review_reserve=3))
    source=control['geometry_source']
    atomic_json(worker/'accepted_l2_network.json',dict(network=dict(reactions=[dict(name=match['reaction'],product_complex=source['product_complex'])])))
    dbpath=worker/'kinbot.db';connect(dbpath).write(Atoms(match['symbols'],positions=GEOMS[name]),
        name=source['record_name'],data=dict(status='normal'))
    before={p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    result=helper.prepare(request,case/'input.json',recipe,tmp_path/'screen')
    assert result['tasks']==8
    data=helper.recon._load(tmp_path/'screen/input.json')
    assert data['_probe_space']['active_electrons']==19 and data['_probe_space']['active_orbitals']==14
    assert np.allclose(data['geometry_angstrom'],GEOMS[name])
    for anchor in data['anchors']:
        xyz=helper.recon._anchor_geometry(data,anchor);i,j=[i-1 for i in anchor['bond']]
        assert np.linalg.norm(xyz[i]-xyz[j])==pytest.approx(anchor['distance_angstrom'])
        for fragment in source['fragments']:
            m=np.array(fragment['atoms'])-1
            dist=lambda a:np.linalg.norm(a[:,None]-a[None,:],axis=2)
            assert np.allclose(dist(xyz[m]),dist(np.array(GEOMS[name])[m]))
    assert all(p.read_bytes()==v for p,v in before.items())
    bad=dict(source,fragments=[dict(source['fragments'][0],atoms=[1,2,3]),dict(source['fragments'][1],atoms=[4,5,6])])
    with pytest.raises(ValueError,match='fragments'):
        helper.product_complex_geometry(request,json.loads(request.read_text()),bad)
    connect(dbpath).write(Atoms(match['symbols'],positions=GEOMS[name]),
        name=source['record_name'],data=dict(status='error'))
    with pytest.raises(ValueError,match='normal completion'):
        helper.product_complex_geometry(request,json.loads(request.read_text()),source)
