"""Native-deck staging against supplied request and L2 product evidence."""
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot import requested_orbital_recon as review
from kinbot.execution import atomic_json

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_remaining_20261025'
PARENT = [[.5846511503115959,-.4473894706201628,-.2784601779532531],
 [.575829259100462,.08822847380033334,1.3556121135846924],
 [-.3934694339188291,-.4816026987731529,2.9599609456852365],
 [.4488465799514408,1.4374948831700918,1.2061461932292064],
 [-.12573832349725744,-1.2795037814491745,-.3408185512544808],
 [1.908852058798488,.1575620792914991,1.6336098268919377]]
PRODUCTS = {
 'csh2_sh_exit': [[.7219452853688184,.6876296231284607,-.8231308537479903],
 [.5751091898167605,.40464295171318393,.813094484471407],
 [-1.0487204147950917,-.46389456781748845,2.482718849038955],
 [.7686489724041617,1.4705786440049884,1.637702392437019],
 [-1.2537222988035466,-1.539500053118284,1.7056717485427364],
 [1.5632539313764098,-.38749478263098086,1.2977340922041016]],
 'hcs_h2s_exit': [[-.7156168281814114,-1.0245318258574119,.7932996495800226],
 [2.3480634840356265,1.0798024993655313,1.3372570594220932],
 [-1.9196034459884892,-.03322257557644806,.669038534624884],
 [1.3107377882972948,.24491280314285946,1.1433413869850342],
 [-.6619467942545391,-2.1089809104382637,.7507393629889204],
 [1.5579519810492546,2.161207982881922,1.2996443810404459]]}


@pytest.fixture
def tree(tmp_path):
    production = tmp_path/'production'; case = production/'ch3ss'
    atomic_json(case/'input.json', dict(queue_name='month-long-cpu'))
    atomic_json(production/'.kinbot_resource_policy.json', dict(schema=1, max_jobs=23,
        excluded_nodes=['node34'], revision=14, manage_pending=True))
    rows = []
    for entry in json.loads((PLAN/'screens.json').read_text())['screens']:
        recipe = json.loads((PLAN/entry['recipe']).read_text()); match = recipe['request_match']
        worker = case/match['parent']
        atomic_json(worker/'vrctst/orbital_review'/match['reaction']/'request.json',
                    dict(match, schema=1, status='needs_review', positions=PARENT))
        if entry['name'] in PRODUCTS:
            source = recipe['geometry_source']
            rows.append(dict(name=match['reaction'], product_complex=source['product_complex']))
            connect(worker/'kinbot.db').write(Atoms(match['symbols'], positions=PRODUCTS[entry['name']]),
                name=source['record_name'], data=dict(status='normal'))
    atomic_json(worker/'accepted_l2_network.json', dict(network=dict(reactions=rows)))
    return production, production/'orbital_recon'/'remaining'


def test_all_four_actual_requests_stage_and_keep_distinct_product_exits(tree, monkeypatch):
    production, output = tree
    before = {p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('Staging must not submit'))
    result = review.batch(production, PLAN/'screens.json', output, shared_pool=True)
    assert (result['screens'], result['tasks'], result['max_jobs']) == (4,32,3)
    assert all(p.read_bytes() == b for p,b in before.items())
    for child in review._saved(output)['reviews']:
        child = Path(child); recipe = json.loads((child/'recipe.json').read_text())
        loaded = review.recon._load(child/'input.json')
        manifest = json.loads((child/'native/reconnaissance.json').read_text())
        assert manifest['shared_policy_file'] == str(production/'.kinbot_resource_policy.json')
        assert all(not j.get('job_id') and j.get('orbital_predecessor') is None for j in manifest['jobs'])
        for job in manifest['jobs']:
            text = (Path(job['directory'])/(job['name']+'.inp')).read_text().lower()
            assert '{multi,cas-ci;' in text and 'rs2' not in text
        if child.name in PRODUCTS:
            assert loaded['_probe_space']['active_orbitals'] == 15
            assert np.allclose(loaded['geometry_angstrom'], PRODUCTS[child.name])
            for anchor in loaded['anchors']:
                xyz = review.recon._anchor_geometry(loaded, anchor)
                i,j = np.array(anchor['bond'])-1
                assert np.linalg.norm(xyz[i]-xyz[j]) == pytest.approx(anchor['distance_angstrom'])
                for fragment in recipe['geometry_source']['fragments']:
                    m = np.array(fragment['atoms'])-1
                    delta = lambda a: np.linalg.norm(a[:,None]-a[None,:], axis=2)
                    assert np.allclose(delta(xyz[m]), delta(np.array(PRODUCTS[child.name])[m]))
        else:
            assert np.allclose(loaded['geometry_angstrom'], PARENT)
            assert loaded['_probe_space']['active_orbitals'] == 12
    with pytest.raises(FileExistsError, match='resume'):
        review.batch(production, PLAN/'screens.json', output, shared_pool=True)


def test_bad_last_product_record_cannot_partially_publish_campaign(tree):
    production, output = tree
    worker = production/'ch3ss/791871710370000000002'
    source = json.loads((PLAN/'hcs_h2s_exit.json').read_text())['geometry_source']
    connect(worker/'kinbot.db').write(Atoms('CSSHHH', positions=PRODUCTS['hcs_h2s_exit']),
        name=source['record_name'], data=dict(status='error'))
    with pytest.raises(ValueError, match='normal completion'):
        review.batch(production, PLAN/'screens.json', output, shared_pool=True)
    assert not output.exists()


def test_campaign_launch_and_resume_share_retained_production_capacity(tree, monkeypatch):
    """Simulated site: terminal audits are stubbed, scientific acceptance is not tested."""
    from types import SimpleNamespace
    from kinbot import run_policy
    from kinbot import run_walltime
    production, output = tree
    policy_path = production/'.kinbot_resource_policy.json'
    policy = json.loads(policy_path.read_text()); policy.update(fair_submission=False)
    atomic_json(policy_path, policy)
    review.batch(production, PLAN/'screens.json', output, shared_pool=True)
    before = policy_path.read_bytes()
    running = {str(100+i): 1 for i in range(23)}
    done, launched = set(), []
    monkeypatch.setattr(run_policy, '_queue', lambda root: running.copy())
    monkeypatch.setattr(run_policy, 'slurm_job_record', lambda j:
        {'JobId':j, 'JobState':'COMPLETED'} if j in done else {})
    monkeypatch.setattr(run_walltime, 'partition_walltime', lambda p: '31-00:00:00')
    def scheduler(command, **kwargs):
        if command[0] == 'bash':
            return SimpleNamespace(returncode=0, stdout='/mock/molpro\n', stderr='')
        if command[0] == 'squeue':
            return SimpleNamespace(returncode=0, stdout='\n'.join(running), stderr='')
        assert command[0] == 'sbatch'  # Never scancel or a production worker command.
        assert len(running) < 23
        ident = str(1000+len(launched))
        launched.append((ident, str(kwargs['cwd'])))
        running[ident] = 1
        return SimpleNamespace(returncode=0, stdout=ident+'\n', stderr='')
    monkeypatch.setattr(review.subprocess, 'run', scheduler)
    def audit(path, **kwargs):
        state = json.loads(path.read_text())
        for job in state['jobs']:
            job['status'] = 'complete' if job.get('job_id') in done else 'incomplete'
        state['status'] = ('complete_review_required' if all(j['status']=='complete' for j in state['jobs'])
                           else 'incomplete')
        atomic_json(path, state)
        return state
    monkeypatch.setattr(review.recon, '_audit', audit)
    monkeypatch.setattr(review.recon, '_recover_slurm_startup_failures', lambda *a: None)
    monkeypatch.setattr(review.recon, '_retry_allocation_failures', lambda p,m,*a: m)
    def finish(_):
        if '100' in running:  # One production calculation naturally completes.
            running.pop('100'); done.add('100')
        for ident,_ in launched:
            if ident in running:
                running.pop(ident); done.add(ident)
    monkeypatch.setattr(review.recon.time, 'sleep', finish)
    result = review.run(output, report_failures=True)
    assert result['status'] == 'complete_review_required' and result['tasks'] == 32
    assert len(launched) == len({d for _,d in launched}) == 32
    assert running == {str(100+i):1 for i in range(1,23)}
    assert policy_path.read_bytes() == before
    assert json.loads((production/'.kinbot_policy_admissions.json').read_text()) == {}
    review.run(output, report_failures=True)
    assert len(launched) == 32
