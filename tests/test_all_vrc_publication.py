"""Full five-channel publication and diagnostic shutdown, with synthetic native files."""
from copy import deepcopy
import json
from pathlib import Path
import re

import pytest
from ase import Atoms
from ase.db import connect

from kinbot import input_refresh, run_policy, pool_drain
from kinbot.execution import atomic_json, file_sha256
from tests.test_requested_orbital_recon import completed_screen, helper
from tests.test_requested_orbital_batch import batch_tree, selected_batch, accepted_requests, live_cases, SS
from tests.test_product_exit_recon import PLAN, GEOMS
from tests.test_queue_progress_review import queue


@pytest.fixture
def all_reviews(batch_tree):
    root,_=batch_tree
    policy=helper._read(root/run_policy.NAME)
    atomic_json(root/run_policy.NAME,dict(policy,revision=1))
    old=root.parent/'r7_additional_orbital_recon_20261021'
    selected_batch((root,old));accepted_requests(root,old)
    review=root.parent/'r7_progress_review_20261024'
    xyz=[[1.8122818351,.0029772478,-.0010645061],[2.2322433887,1.7849627863,-.0036765215],
        [.0180913054,.0086000664,-.0026139346],[1.7933644276,2.0266677332,-1.2480033488],
        [2.2475352842,-.4907718979,-.8623412662],[2.2307468144,-.3955252319,.9183373719]]
    rows=[]
    for name in ('ch2_sh_s_cs','hcs2_h2_exit','ch2s_sh_exit'):
        recipe=helper._read(PLAN/(name+'.json'));match=recipe['request_match']
        worker=root/'ch3ss'/match['parent']
        request=worker/'vrctst/orbital_review'/match['reaction']/'request.json'
        atomic_json(request,dict(match,schema=1,status='needs_review',positions=xyz))
        if name in GEOMS:
            source=recipe['geometry_source']
            rows.append(dict(name=match['reaction'],product_complex=source['product_complex']))
            connect(worker/'kinbot.db').write(Atoms(match['symbols'],positions=GEOMS[name]),
                name=source['record_name'],data=dict(status='normal'))
            atomic_json(worker/'accepted_l2_network.json',dict(network=dict(reactions=rows)))
        completed_screen((request,root/'ch3ss/input.json',review/'screens'/name),PLAN/(name+'.json'))
    atomic_json(review/'screens/review_launch.json',dict(schema=1,files={},
        reviews=[str(review/'screens'/n) for n in ('ch2_sh_s_cs','hcs2_h2_exit','ch2s_sh_exit')],
        production_policy=str(root/run_policy.NAME)))
    candidate=review/'ch2_sh_s_followup/candidates'
    helper.candidates(review/'screens/ch2_sh_s_cs',PLAN/'ch2_sh_s_cs_candidates.json',candidate)
    for system in helper._read(candidate/'review_launch.json')['systems']:
        path=Path(system['manifest']);m=helper._read(path)
        for j in m['jobs']:
            d=Path(j['directory']);wfu=d/(j['name']+'.wfu');out=d/(j['name']+'.out')
            wfu.write_bytes(b'synthetic candidate');out.write_text('Synthetic comparison.\n')
            j.update(status='complete',wfu=str(wfu),wfu_sha256=file_sha256(wfu),output_sha256=file_sha256(out))
        atomic_json(path,m)
    original=(root/'ch3ss'/SS/(SS+'.json')).read_bytes()
    inp=root/'ch3ss/input.json'
    settings=helper._read(inp)
    settings.update(vrc_tst_scan_points=[2.3,3.7,8.1],vrc_tst_fragment_geometry='asymptotic_scan')
    atomic_json(inp,settings)
    live_cases(root);(root/'ch3ss'/SS/(SS+'.json')).write_bytes(original)
    return root,review


def test_all_five_exact_recipes_one_revision_and_repeat(all_reviews,monkeypatch):
    from kinbot.molpro import _render_vrc_state_model
    root,review=all_reviews
    before=(root/'ch3ss/input.json').read_bytes()
    production_ch3s={p:p.read_bytes() for p in (root/'ch3s').rglob('*') if p.is_file()}
    plan=PLAN/'all_selections.json'
    preview=helper.review_batch(root,root.parent,plan,require_all=True)
    assert len(preview['cases']['ch3ss']['channels'])==5
    assert (root/'ch3ss/input.json').read_bytes()==before
    original=input_refresh.submit;publications=[]
    def submit(*args,**kwargs):
        if kwargs.get('apply'):publications.append(kwargs['proposal'])
        return original(*args,**kwargs)
    monkeypatch.setattr(input_refresh,'submit',submit)
    result=helper.review_batch(root,root.parent,plan,require_all=True,apply=True)
    assert result['cases']['ch3ss']['status']=='queued' and len(publications)==1
    receipt=input_refresh.receipt(root/'ch3ss')['id']
    helper.review_batch(root,root.parent,plan,require_all=True,apply=True)
    assert input_refresh.receipt(root/'ch3ss')['id']==receipt
    assert all(p.read_bytes()==data for p,data in production_ch3s.items())
    data=helper._read(root/'ch3ss/input.json')
    assert sum(len(m['channels']) for m in data['vrc_tst_surface_models'].values())==5
    for parent,m in data['vrc_tst_surface_models'].items():
        for channel in m['channels']:
            b=channel['bottlenecks'][0];p=b['parameters'];a=p['vrc_tst_active_space']
            assert len(b['states'])==1 and b['states'][0]['root']==1
            assert p['vrc_tst_correction_model']=='direct'
            for model in p['vrc_tst_state_model']['basis_models'].values():
                ne,no=a['active_electrons'],a['active_orbitals'];closed=(41-ne)//2
                _,deck,_=_render_vrc_state_model(model,nelectron=41,active_electrons=ne,
                    active_orbitals=no,spin=1)
                # Simulate Molpro's swaps: noncontiguous H2 orbitals must end
                # up inside the active block, with no duplicate/missing source.
                order=list(range(100))
                for left,right in re.findall(r'rotate,(\d+)\.1,(\d+)\.1,0',deck):
                    i,j=int(left),int(right);order[i],order[j]=order[j],order[i]
                assert order[closed+1:closed+no+1]==model['orbital_selection']['active_indices']
                assert 'start,2140.2,type=natural' in deck
                assert f'occ,{closed+no}\nclosed,{closed}' in deck
    archive=next((root/'ch3ss/orbital_reviews/live_refresh').iterdir())
    selection=helper._read(archive/'selection.json')['selections'][2]
    assert len(selection['supporting_evidence'][0]['files'])==48
    summary=helper._read(archive/'review_summary.json')
    for row in summary['channels']:
        scan=row['scan_plan']
        assert scan['correction_distances_angstrom']==[2.3,3.7,8.1]
        assert scan['asymptote_distance_angstrom']==30.
        is_exit=Path(row['source_screen']).name in ('hcs2_h2_exit','ch2s_sh_exit')
        assert scan['sampling_coordinate']==('com' if is_exit else 'reactive')
        assert scan['sampling_distances_angstrom']==data['rotdpy_dist']
        if is_exit:
            assert scan['geometry_source'] and scan['anchor_atoms_one_based'] in ([6,4],[5,4])


@pytest.mark.parametrize('failure',['incomplete','changed','supporting_incomplete','supporting_changed'])
def test_batch_has_no_partial_write_or_binding(all_reviews,monkeypatch,failure):
    root,review=all_reviews
    if failure.startswith('supporting'):
        systems=helper._read(review/'ch2_sh_s_followup/candidates/review_launch.json')['systems']
        path=Path(systems[-1]['manifest'])
    else:path=review/'screens/ch2s_sh_exit/native/reconnaissance.json'
    m=helper._read(path);j=m['jobs'][-1]
    if failure.endswith('incomplete'):
        j['status']='submitted';atomic_json(path,m)
    else:Path(j['wfu']).write_bytes(b'changed')
    before={p:p.read_bytes() for base in (root,root.parent/'r7_additional_orbital_recon_20261021') for p in base.rglob('*.json')}
    if failure.endswith('changed'):
        with pytest.raises(ValueError,match='changed|unchanged'):
            helper.review_batch(root,root.parent,PLAN/'all_selections.json',require_all=True,apply=True)
    else:
        r=helper.review_batch(root,root.parent,PLAN/'all_selections.json',require_all=True,apply=True)
        assert r['cases']['ch3ss']['status']=='waiting'
    assert all(p.read_bytes()==data for p,data in before.items())


def test_retire_stops_submitters_before_jobs_preserves_evidence_and_restores_pool(all_reviews,monkeypatch):
    root,review=all_reviews
    evidence={p:p.read_bytes() for p in review.rglob('*') if p.is_file()}
    production={p:p.read_bytes() for p in root.rglob('*') if p.is_file() and p.name!=run_policy.NAME}
    events=[]
    monkeypatch.setattr(queue,'stop_owners',lambda *a:events.append('stop_owners'))
    monkeypatch.setattr(queue.run_cleanup,'cancel_allocations',lambda p: (
        events.append('cancel_allocations') or dict(cancelled_job_ids=['123'])))
    monkeypatch.setattr(pool_drain,'pulse',lambda p,force:events.append('drain') or dict(status='reconciled'))
    before=run_policy.read(root)
    for _ in range(2):
        result=queue.retire(root,review)
        assert result['phase']=='stopped' and result['cancelled_job_ids']==['123']
    assert events==['stop_owners','cancel_allocations','drain']*2
    after=run_policy.read(root)
    assert after['max_jobs']==23 and after['orbital_review_reserve']==0
    assert after['excluded_nodes']==before['excluded_nodes']
    assert all(p.read_bytes()==data for p,data in evidence.items())
    assert all(p.read_bytes()==data for p,data in production.items())
    with pytest.raises(ValueError,match='retired'):
        queue.perform(root,review,Path('/unused'),[],replace_owner=True)


@pytest.mark.parametrize('failure',[None,'incomplete','refresh'])
def test_full_publication_handoff_preflight_refresh_publish_export_repeat(all_reviews,monkeypatch,tmp_path,tmp_path_factory,failure):
    import importlib.util
    import sys
    monkeypatch.setitem(sys.modules,'queue_progress_review',queue)
    monkeypatch.syspath_prepend(str(PLAN.parents[1]/'tools'))
    path=PLAN.parents[1]/'tools/publish_progress_review.py'
    spec=importlib.util.spec_from_file_location('publish_progress_review_test',path)
    publish=importlib.util.module_from_spec(spec);spec.loader.exec_module(publish)
    root,review=all_reviews
    exports=tmp_path_factory.mktemp('publication_exports')
    events=[]
    monkeypatch.setattr(queue,'stop_owners',lambda *a:events.append('stop'))
    monkeypatch.setattr(queue.run_cleanup,'cancel_allocations',lambda p:
        events.append('cancel') or dict(cancelled_job_ids=['123']))
    monkeypatch.setattr(pool_drain,'pulse',lambda *a,**k:dict(status='reconciled'))
    # Test the actual archive writer with this fixture's finite evidence plan.
    plan=tmp_path/'plan';plan.mkdir()
    for p in PLAN.glob('*.json'):(plan/p.name).write_bytes(p.read_bytes())
    old=tmp_path/'r7_additional_20261021';old.mkdir()
    for p in (PLAN.parent/old.name).glob('*_selected.json'):(old/p.name).write_bytes(p.read_bytes())
    atomic_json(plan/'export.json',dict(sources=[dict(root='runs',path=review.name),
        dict(root='production',path='ch3ss/input.json')]))
    monkeypatch.setattr(publish,'PLAN',plan)
    original=publish.recon.review_batch
    def batch(*a,**kw):
        events.append('publish' if kw.get('apply') else 'preflight')
        return original(*a,**kw)
    monkeypatch.setattr(publish.recon,'review_batch',batch)
    before=(root/'ch3ss/input.json').read_bytes()
    original_run=publish.subprocess.run
    def refresh(args,**kw):
        if args[:1]!=['bash']:
            return original_run(args,**kw)
        assert args[-2:]==[str(root),'23'] and kw['check']
        assert run_policy.read(root)['orbital_review_reserve']==0
        events.append('refresh')
        if failure=='refresh':raise RuntimeError('simulated refresh failed before publication')
    monkeypatch.setattr(publish.subprocess,'run',refresh)
    if failure=='incomplete':
        manifest=review/'screens/ch2s_sh_exit/native/reconnaissance.json'
        record=helper._read(manifest);record['jobs'][-1]['status']='incomplete';atomic_json(manifest,record)
    if failure:
        with pytest.raises(RuntimeError,match='not all ready|refresh failed'):
            publish.perform(root,exports/'export')
        assert events==['stop','cancel','preflight']+(['refresh'] if failure=='refresh' else [])
        assert (root/'ch3ss/input.json').read_bytes()==before
        assert not (review/'publication_handoff.json').exists()
        return
    for i in range(2):
        result=publish.perform(root,exports/f'export{i}')
        assert result['cases']['ch3ss']['status']=='queued'
        archive=publish.exporter.verify(exports/f'export{i}/orbital-evidence.tar.gz')
        assert not archive['issues'] and any(f['member'].endswith('.wfu') for f in archive['files'])
        revision=input_refresh.receipt(root/'ch3ss')['id']
        if i:assert revision==prior
        prior=revision
    assert events==['stop','cancel','preflight','refresh','publish']*2
