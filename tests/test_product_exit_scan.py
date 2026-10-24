"""Real mapped exit geometries through the production scan loop; native QC doubled."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect
from ase.io import write

from kinbot.execution import atomic_json, file_sha256
from kinbot.stationary_pt import StationaryPoint
from kinbot.vrc_tst_scan import VTS
from tests.test_product_exit_recon import GEOMS, PLAN


@pytest.mark.parametrize('name',list(GEOMS))
def test_outer_exit_keeps_product_maps_anchor_grid_and_auditable_geometries(tmp_path,monkeypatch,name):
    monkeypatch.chdir(tmp_path)
    recipe=json.loads((PLAN/(name+'.json')).read_text())
    xyz=np.array(GEOMS[name]);symbols=recipe['request_match']['symbols']
    maps=[[i-1 for i in f['atoms']] for f in recipe['geometry_source']['fragments']]
    anchors=[i-1 for i in recipe['request_match']['bond']]
    products=[];bonds=[]
    for m,fragment in zip(maps,recipe['geometry_source']['fragments']):
        point=StationaryPoint('fragment',fragment['charge'],fragment['multiplicity'],
            atom=[symbols[i] for i in m],geom=xyz[m])
        point.characterize();products.append(point)
        bonds.extend([[m[i],m[j]] for i in range(len(m)) for j in range(i+1,len(m)) if point.bond01[i,j]])
    # Deliberately different parent geometry: using it would scan the wrong path.
    parent=xyz.copy();parent[:,0]*=2
    reaction=SimpleNamespace(species=SimpleNamespace(atom=symbols,geom=parent,charge=0,mult=2),
        irc_prod=SimpleNamespace(geom=xyz,bondlist=bonds),products=products,maps=maps,
        scan_coo=anchors,equiv=[[a] for a in anchors],do_vdW=True)
    settings=dict(vrc_tst_fragment_geometry='asymptotic_scan',pp_orient='geometric',
        vrc_tst_scan_method='caspt2(7,6)',vrc_tst_scan_basis='avdz-f12',
        vrc_tst_scan_points=[2.3,3.7,8.1],vrc_tst_asymptote_distance=30.)
    vts=VTS(SimpleNamespace(),settings,None);vts.scan_reac['exit']=reaction
    calls=[]
    def native(reaction,geom,distance,logical_name,**kwargs):
        calls.append((distance,np.array(geom),logical_name))
        graph=tmp_path/'vrctst/molpro_scan'/logical_name;task=graph/'tasks/relaxed_scan'
        task.mkdir(parents=True,exist_ok=True)
        a=Atoms(symbols,positions=geom);write(task/'final.xyz',a)
        atomic_json(task/'task.json',{'task':{'id':'relaxed_scan'}})
        atomic_json(task/'execution.json',{'status':'executed'})
        atomic_json(graph/'workflow.json',{'name':logical_name})
        return a,-1.,graph,{'details':{'final_geometry_sha256':file_sha256(task/'final.xyz')}},{'seed':logical_name}
    monkeypatch.setattr(vts,'_molpro_multireference_scan_task',native)
    vts.do_scan(['exit'])
    assert [d for d,g,n in calls]==[30.,2.3,3.7,8.1]
    assert np.allclose(calls[0][1],VTS._shift_scan_geometry(reaction,xyz,30.))
    assert not np.allclose(calls[0][1],VTS._shift_scan_geometry(reaction,parent,30.))
    for d,geom,n in calls:
        assert np.linalg.norm(geom[anchors[0]]-geom[anchors[1]])==pytest.approx(d)
        for m in maps:
            distances=lambda p:np.linalg.norm(p[:,None]-p[None,:],axis=2)
            assert np.allclose(distances(geom[m]),distances(xyz[m]))
    rows=list(connect('kinbot.db').select())
    assert len(rows)==8  # relaxed/frozen at all three distances plus the asymptote
    assert {r.data.vrc_scan_role for r in rows}=={'relaxed','frozen'}
    assert all(r.data.vrc_scan_dispatch for r in rows)
    progress=list(Path('vrctst/scan_progress/exit').rglob('*.json'))
    assert progress
    text='\n'.join(p.read_text() for p in progress)
    assert 'coordinate' in text and 'separation_angstrom' in text
